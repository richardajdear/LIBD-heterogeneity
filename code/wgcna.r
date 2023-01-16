library(WGCNA)
library(SummarizedExperiment)
library(tidyverse)
library(comprehenr)
# library(BioNERO)
options(stringsAsFactors = FALSE)
options(max.print = 100)

outputDir <- "../outputs/"

# Pick soft threshold
get_soft_threshold <- function(rse, assayname = 'ranknorm', threads = NULL) {
    if (!is.null(threads)) {
        WGCNA::allowWGCNAThreads(nThreads = threads)
        WGCNAnThreads()
    }

    # WGCNA takes genes in columns, samples in rows
    exp <- t(assays(rse)[[assayname]])

    # Parameters for SFT
    args <- list(
        data = exp,
        powerVector = 1:25,
        RsquaredCut = 0.8,
        verbose = 5,
        blockSize = 9000,
        networkType = 'signed hybrid',
        moreNetworkConcepts = TRUE,
        corOpt = list(nThreads = threads)
    )

    fit <- do.call(WGCNA::pickSoftThreshold, args)
    return(fit)
}

# Fit WGCNA
fit_WGCNA <- function(rse, power, assayname = "ranknorm", threads = NULL, 
                      saveTOMs = FALSE, loadTOM = FALSE, fileBase = "test",
                      blocks = NULL) {
    if (!is.null(threads)) {
        WGCNA::allowWGCNAThreads(nThreads = threads)
        WGCNAnThreads()
    }

    # WGCNA takes genes in columns, samples in rows
    exp <- t(assays(rse)[[assayname]])

    # Get power estimate
    # power <- sft$powerEstimate

    # If using preloaded TOM, also try to read in blocks using same file name
    if (loadTOM && is.null(blocks)) {
        tryCatch(
            {
                blocksFile <- paste0(outputDir, fileBase, "-blocks.rds")
                print(paste("Reading blocks from", blocksFile))
                blocks <- readRDS(blocksFile)
            },
            error = function(condition) {
                message(condition)
                message("Will recalculate blocks.")
            }
        )
    }

    # Parameters for WGCNA
    args <- list(
        datExpr = exp,
        randomSeed = 123, verbose = 5, maxBlockSize = 9000,
        TOMType = "signed",
        loadTOM = loadTOM, blocks = blocks,
        saveTOMs = saveTOMs,
        saveTOMFileBase = paste0(outputDir, fileBase, "-TOM"),
        power = power,
        networkType = 'signed hybrid',
        minModuleSize = 40, ### changed
        pamStage = TRUE, pamRespectsDendro = TRUE,
        deepSplit = 4, ###
        reassignThreshold = 1e-06,
        mergeCutHeight = 0.15,
        numericLabels = FALSE,
        # minCoreKME = 0.5, ### changed
        # minCoreKMESize = 0, ### changed
        minKMEtoStay = 0.3
    )

    # Fit WGCNA by blocks
    net <- do.call("blockwiseModules", args)

    # Save blocks if it was null
    if (is.null(blocks) && saveTOMs) {
        blocksFile <- paste0(outputDir, fileBase, "-blocks.rds")
        print(paste("Saving blocks to", blocksFile))
        saveRDS(net$blocks, blocksFile)
    }

    # Fit intramodular connectivity
    args_IMC <- list(
        datExpr = exp,
        networkType = 'signed hybrid',
        colors = net$colors,
        power = power,
        getWholeNetworkConnectivity = TRUE
    )
    IMC <- do.call("intramodularConnectivity.fromExpr", args_IMC)
    IMC$modules <- net$colors
    rownames(IMC) <- names(net$colors)
    net$IMC <- IMC

    # Add KME to output
    net$KME <- WGCNA::signedKME(exp, net$MEs)

    return(net)
}


# Match the modules from one network to another and relabel
match_modules <- function(source_net, target_net) {

    # Match color of each gene from source to target
    matched_colors <- WGCNA::matchLabels(source_net$colors, target_net$colors)
    names(matched_colors) <- names(source_net$colors)

    # Get mapping from each old color to new
    old <- unique(source_net$colors)
    new <- unique(matched_colors)

    # Relabel the variables in source net
    source_net$colors <- matched_colors %>% c
    source_net$IMC$modules <- matched_colors %>% c
    ME_mapping <- setNames(paste0('ME', new), paste0('ME', old))
    source_net$MEs <- source_net$MEs %>% rename(ME_mapping)
    kME_mapping <- setNames(paste0('kME', new), paste0('kME', old))
    source_net$KME <- source_net$KME %>% rename(kME_mapping)

    return(source_net)
}


filter_kME_one_module <- function(kME, name, colors) {
    # Clean name
    name <- str_replace(name, "kME", "")
    # Set to NA if not this gene's module
    kME[colors != name] <- NA
    return(kME)
}
filter_kME <- function(net) {
    kME_filtered <- mapply(filter_kME_one_module,
                          net$KME, names(net$KME),
                          MoreArgs = list(colors=net$colors))
    order <- net$colors %>% table %>% sort
    order <- paste0("kME", names(order)) %>% rev
    rownames(kME_filtered) <- rownames(net$KME)
    return(kME_filtered[, order])
}


# Split rse data in two halves
split_samples <- function(rse) {
    n_samples <- dim(rse)[2]
    mask <- sample(c(TRUE, FALSE), size = n_samples, replace = TRUE)
    return(list(rse[, mask], rse[, !mask]))
}


# Fit WGCNA on a list of rse, using the max SFT power
make_nets <- function(rse_list, cor_method = "pearson", assayname = NULL) {
    # Select the right assay if specified (otherwise will use first assay)
    if (!is.null(assayname)) {
        for (i in seq_along(rse_list)) {
            assay(rse_list[[i]]) <- assays(rse_list[[i]])[[assayname]]
        }
    }
    # Find soft threshold power for each split 
    sft_list <- lapply(rse_list, function(rse) {
        rse %>% get_soft_threshold(threads = 9)
    })
    # Get max of powers
    power <- to_vec(for(sft in sft_list) sft$power) %>% max
    # Fit nets
    net_list <- lapply(rse_list, function(x) fit_WGCNA(x, power, threads=9))

    return(net_list)
}

# Fit WGCNA to multiple splits
fit_net_splits <- function(rse, n_splits = 2) {
    nets <- list()
    for (i in 1:n_splits) {
        print(paste0("..fitting split ", i, "/", n_splits))
        rse_split <- rse %>% split_samples
        invisible(capture.output(
            nets[[i]] <- make_nets(rse_split)
        ))
    }
    print("Done.")
    return(nets)
}













# Fit WGCNA with a given power
# Include the expression data and kME in the output
# fit_wgcna <- function(rse, power, cor_method = "pearson") {
#     net <- rse %>% exp2gcn(
#         SFTpower = power,
#         net_type = "signed hybrid",
#         cor_method = cor_method
#     )
#     # Don't save full expression mat to save data
#     # net$expression <- assay(rse) %>% `row.names<-`(rownames(rse))

#     # Add kME to the output
#     exp <- assay(rse) %>% `row.names<-`(rownames(rse))
#     net$kME <- get_kME(net, exp)
#     return(net)
# }

# Add kME matrix to each net object
# get_kME <- function(net, exp) {
#     # exp <- net$expression
#     MEs <- net$MEs
#     cor_method <- net$params$cor_method
#     # Calculate kME appropiately
#     if(cor_method == "spearman") {
#         MM <- WGCNA::signedKME(
#             t(exp), MEs, outputColumnName = "",
#             corOptions = "use = 'p', method = 'spearman'"
#         )
#     } else if(cor_method == "pearson") {
#         MM <- WGCNA::signedKME(t(exp), MEs, outputColumnName = "")
#     } else if(cor_method == "biweight") {
#         MM <- WGCNA::signedKME(
#             t(exp), MEs, corFnc = "bicor", outputColumnName = "",
#             corOptions = "maxPOutliers = 0.1"
#         )
#     }
#     # Set genes as rownames
#     rownames(MM) <- rownames(exp)
#     return(MM)
# }