library(jaffelab)
library(recount)
library(matrixStats)
library(BRETIGEA)
library(preprocessCore)
library(qsvaR)
library(RNOmni)
library(WGCNA)
library(tidyverse)

# Cleanup covariates like 'concordMapRate', 'overallMapRate', 'mitoRate' etc
# and recount RPKM
# (data from one sample may include multiple ‘lanes’)
merge_covariates <- function(rse) {
    rse_merged <- jaffelab::merge_rse_metrics(rse)
    assays(rse_merged)$RPKM     <- recount::getRPKM(
        rse_merged, length_var = "Length", mapped_var = "numMapped")
    assays(rse_merged)$logRPKM  <- log2(assays(rse_merged)$RPKM + 1)

    # Take mean of RIN numbers
    rse_merged$RIN <- sapply(rse_merged$RIN, mean)
    return(rse_merged)
}

# Select samples
select_samples <- function(rse) {
    minAge = -1; maxAge = 150; minRin = 6;
    diagnosis = c("Control"); race = c("AA","CAUC")

    rse_subset <- rse %>% SummarizedExperiment::subset(
        subset = TRUE, # keep all rows (genes)
        select =
            (Age                >=   minAge       )   &
            (Age                <=   maxAge       )   &
            (Race              %in%  race         )   &
            (Dx                %in%  diagnosis    )   &
            (sapply(RIN, mean)  >=   minRin       )   &
            (!is.na(snpPC1)                       )
    )
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
    print(paste("Removed", sum(outlier_samples), "outlier samples"))

    rse <- rse[, as.vector(!outlier_samples)]
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

    rse <- rse %>% SummarizedExperiment::subset(
        subset = (
            (gencodeID %in% genes_to_keep) &
            # Also remove gene symbol 'MT-' (for mitochondrial genes)
            (!grepl("^MT-",Symbol))
        )
  )
  return(rse)
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
    return(rse)
}

# Normalize samples (quantile norm)
normalise_samples <- function(rse) {
    # Quantile normalise samples. 
    # Ref: http://jtleek.com/genstats/inst/doc/02_05_normalization.html
    exp <- assays(rse, withDimnames = FALSE)$logRPKM
    exp_norm <- preprocessCore::normalize.quantiles(exp)
    assays(rse, withDimnames = FALSE)$logRPKM_quant_normalised <- exp_norm
    return(rse)
}

# Regress out covariates
regress_covariates <- function(rse, rse_tx, 
                               assayname = 'logRPKM_quant_normalised',
                               formula = "~ Age + Sex + Race", 
                               n_qsvs = NULL,
                               method = 'libd', P = 0) {
    source("qsva.r")
    rse_tx_matched <- match_rse_tx_samples(rse_tx, rse)
    qsvs <- make_qsvs(rse_tx_matched, formula = formula)
    rse_clean <- get_residuals(rse, assayname, formula, qsvs, n_qsvs = n_qsvs)
    return(rse_clean)
}

# Normalize genes (rank norm)
# https://www.rdocumentation.org/packages/RNOmni/versions/1.0.0/topics/RankNorm
normalise_genes <- function(rse,
                            assayname = 'residuals') {
    exp <- assays(rse, withDimnames = FALSE)[[assayname]]
    exp_ranknorm <- t(apply(exp, 1, RNOmni::RankNorm))
    assays(rse, withDimnames = FALSE)$ranknorm <- exp_ranknorm
    return(rse)
}