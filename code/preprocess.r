library(glue)
library(jaffelab)
library(recount)
library(matrixStats)
library(BRETIGEA)
library(preprocessCore)
library(qsvaR)
library(RNOmni)
library(WGCNA)
library(tidyverse)

source("../code/qsva.r")

# Cleanup covariates like 'concordMapRate', 'overallMapRate', 'mitoRate' etc
# and recount RPKM (data from one sample may include multiple ‘lanes’)
merge_covariates <- function(rse) {
    rse_merged <- jaffelab::merge_rse_metrics(rse)
    assays(rse_merged)$RPKM     <- recount::getRPKM(
        rse_merged, length_var = "Length", mapped_var = "numMapped")
    assays(rse_merged)$logRPKM  <- log2(assays(rse_merged)$RPKM + 1)

    # Take mean of RIN numbers
    rse_merged$RIN <- sapply(rse_merged$RIN, mean)

    print(glue("Merged data from multi-lane samples."))
    return(rse_merged)
}

# Select samples
select_samples <- function(rse,
    minAge = 17, # Min SCZ diagnosis age is 17
    maxAge = 150,
    minRin = 6,
    region = c("DLPFC"), # First look only at DLPFC
    race = c("AA", "CAUC"),
    diagnosis = c("Control")
) {
    rse_subset <- rse %>% SummarizedExperiment::subset(
        subset = TRUE, # keep all rows (genes)
        select =
            (Age                >=   minAge       )   &
            (Age                <=   maxAge       )   &
            (Region            %in%  region       )   &
            (Race              %in%  race         )   &
            (Dx                %in%  diagnosis    )   &
            (sapply(RIN, mean)  >=   minRin       )   &
            (!is.na(snpPC1)                       )
    )

    print(glue("Selected {dim(rse_subset)[2]} of {dim(rse)[2]} samples."))
    return(rse_subset)
}

# Filter outlier samples
clean_outlier_samples <- function(rse, sd_threshold = 3) {
    exp <- assays(rse)$logRPKM
    # Sample-sample pairwise euclidean distance/disimilarity (based on all genes)
    IAC <- as.matrix(dist(t(exp)))
    # For each gene, get z-score of mean distance
    mean_IAC_zscore <- scale(matrixStats::rowMeans2(IAC))
    outlier_samples <- abs(mean_IAC_zscore) > sd_threshold
    names(outlier_samples) <- rownames(IAC)

    rse <- rse[, as.vector(!outlier_samples)]

    print(glue("Removed {sum(outlier_samples)} outliers leaving {dim(rse)[2]} samples."))
    return(rse)
}

# Filter genes for unexpressed and outliers, and mitochondrial genes
clean_genes <- function(rse) {
    exp <- assays(rse)$RPKM
    median_exp <- matrixStats::rowMedians(exp)
    names(median_exp) <- rownames(exp)

    # Gene filter 1: not more than 20% samples are zero
    gene_filter_1 <- as.vector(rowSums(exp == 0) <= ncol(exp) * 0.2)
    # Gene filter 2: median exp > 0.1
    gene_filter_2 <- as.vector(median_exp >= .1)
    median_exp_filtered_12 <- median_exp[gene_filter_1 & gene_filter_2]
    # Gene filter 3: z-score of median_exp < 3
    gene_filter_3 <- as.vector(abs(scale(median_exp_filtered_12)) <  3)
    genes_to_keep <- names(median_exp_filtered_12)[gene_filter_3]

    rse_clean <- rse %>% SummarizedExperiment::subset(
        subset = (
            (gencodeID %in% genes_to_keep) &
            # Also remove gene symbol 'MT-' (for mitochondrial genes)
            (!grepl("^MT-",Symbol))
        )
  )

  print(glue("Filtered {dim(rse)[1]-dim(rse_clean)[1]} of {dim(rse)[1]} genes."))
  return(rse_clean)
}

get_cell_proportions <- function(rse, assayname = 'RPKM', method = "SVD", nmarkers = 50) {
    exp <- assays(rse)[[assayname]]
    rownames(exp) <- rowData(rse)$Symbol   #Convert ENSEMBL names to GeneSymbol
    capture.output({
        #Scaled estimates
        zz1 <- BRETIGEA::brainCells(exp, nMarker=nmarkers, 
                species="human", method = method, scale = T)
        #Unscaled estimates
        zz2 <- BRETIGEA::brainCells(exp, nMarker=nmarkers, 
                species="human", method = method, scale = F)
        colnames(zz2) <- paste0(colnames(zz2),"_unscaled")
    })
    zz <- cbind(zz1,zz2)

    colData(rse) <- cbind(colData(rse), zz)

    print(glue("Added cell proportion estimates from BRETIGEA."))
    return(rse)
}

# Normalize samples (quantile norm)
normalise_samples <- function(rse) {
    # Quantile normalise samples. 
    # Ref: http://jtleek.com/genstats/inst/doc/02_05_normalization.html
    exp <- assays(rse, withDimnames = FALSE)$logRPKM
    exp_norm <- preprocessCore::normalize.quantiles(exp)
    assays(rse, withDimnames = FALSE)$logRPKM_quant_normalised <- exp_norm

    print(glue("Normalised logRPKM quantiles using preprocessCore."))
    return(rse)
}

# Regress out covariates
regress_covariates <- function(rse, rse_tx, 
                               assayname = 'logRPKM_quant_normalised',
                               formula = "~ Age + Sex + Race",
                               n_qsvs = NULL,
                               method = 'libd', P = 0) {
    source("../code/qsva.r")
    print(glue("Regressing covariates..."))

    rse_tx_matched <- match_rse_tx_samples(rse_tx, rse)
    print(glue("... matched gene & transcript samples"))
    qsvs <- make_qsvs(rse_tx_matched, formula = formula)
    print(glue("... made qSVs"))
    rse_clean <- get_residuals(rse, assayname, formula, qsvs, n_qsvs = n_qsvs)
    print(glue("... computed residuals"))

    print(glue("Regressed covariates."))
    return(rse_clean)
}

# Normalize genes (rank norm)
# https://www.rdocumentation.org/packages/RNOmni/versions/1.0.0/topics/RankNorm
normalise_genes <- function(rse,
                            assayname = 'residuals') {
    exp <- assays(rse, withDimnames = FALSE)[[assayname]]
    exp_ranknorm <- t(apply(exp, 1, RNOmni::RankNorm))
    assays(rse, withDimnames = FALSE)$ranknorm <- exp_ranknorm

    print(glue("Normalised genes by rank using RNOmni."))
    return(rse)
}