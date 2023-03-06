library(biomaRt)
library(BiocParallel)

count_matches <- function(A, B) {
    modules <- names(A$counts)
    gene_list_A <- names(A$colors) %>% split(A$colors)
    gene_list_B <- names(B$colors) %>% split(B$colors)
    
    n_match <- setNames(rep(NA, length(modules)), modules)
    pct_match <- setNames(rep(NA, length(modules)), modules)
    for (module in modules) {
        if (exists(module, where=gene_list_B)) {
            matched <- sum(gene_list_A[[module]] %in% gene_list_B[[module]])
        } else {
            matched <- 0
        }
        total <- length(gene_list_A[[module]])
        n_match[module] <- matched
        pct_match[module] <- matched/total
    }
    return(data.frame(n_match, pct_match))
}


get_topN_matches_one_module <- function(source_weights, target_weights, 
                                        interval=1, how='topN') {
    # Sort the genes in source weights (dropping NAs)
    source_sorted <- order(-source_weights, na.last='NA') %>% suppressWarnings
    # Rank the genes in target weights (keeping NAs)
    target_ranked <- rank(-target_weights, na.last='keep')
    # Sort the target ranks by the source order
    target_ranked_sorted <- target_ranked[source_sorted]
    # Create an empty array to hold counts for top N, at intervals
    points <- seq(interval,length(source_sorted), interval)
    arr <- rep(NA, length(points))
    for (i in points) {
        if (how == 'topN') {
            # For each i, check how many genes with source_rank<i have target_rank<i
            arr[i/interval] <- sum(target_ranked_sorted[1:i] <= i, na.rm=TRUE)
        } else {
            # For each i, check how many genes with source_rank<i have target_rank not NA
            arr[i/interval] <- sum(!is.na(target_ranked_sorted[1:i]))
        }
    }
    return(arr)
}

# Count how many of the topN genes of source module are in topN genes of target module
get_topN_matches <- function(source_net, target_net, 
                             n_modules=20, interval=25, how='topN') {
    # Match target modules to source modules
    # NB: unmatched modules are moved to the back
    target_matched <- target_net %>% match_modules(source_net)
    # Get kME of only within-module genes, drop grey module
    kME_source <- source_net %>% filter_kME() %>% .[, colnames(.)!="kMEgrey"]
    kME_target <- target_matched %>% filter_kME() %>% .[, colnames(.)!="kMEgrey"]
    # Re-order target to match source, only keeping matched modules
    matched_modules <- colnames(kME_source)[colnames(kME_source) %in% colnames(kME_target)]
    kME_target <- kME_target[, matched_modules]
    # Move unmatched modules in source to back
    unmatched_modules <- colnames(kME_source)[!(colnames(kME_source) %in% matched_modules)]
    kME_source <- kME_source[, c(matched_modules, unmatched_modules)]

    # For each module in source, run count topN matching function
    # NB: must convert to dataframe to use with mapply
    topN_list <- mapply(get_topN_matches_one_module,
                        data.frame(kME_source[, 1:n_modules]),
                        data.frame(kME_target[, 1:n_modules]),
                        MoreArgs = list(interval=interval, how=how)
                    )
    # Combine list into dataframe, filling NAs
    topN_df <- data.frame(lapply(topN_list, `length<-`, max(lengths(topN_list))))
    rownames(topN_df) <- seq(1, nrow(topN_df))*interval
    return(topN_df)
}


### Compute enrichment

# Get GO annotations
get_GO_annotations <- function() {
    ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
    atts <- listAttributes(ensembl)
    GO <- getBM(mart=ensembl, attributes=c('ensembl_gene_id','name_1006'))
    GO <- GO[GO[2] != '', ]
    return(GO)
}

# Get enrichments for top n modules
get_enrichments <- function(net, annotation, p = 0.05,
                            module_names = NULL,
                            bp_param = BiocParallel::SerialParam()) {
    # Define background genes, clean names of version
    all_genes <- names(net$colors) %>% str_replace("\\..*", "")

    # Match annotations to background
    annotation_matched <- annotation[annotation[,1] %in% all_genes, ]
    annotation_list <- split(annotation_matched[, 1], annotation_matched[, 2])

    # Get modules gene lists as a list
    modules <- split(all_genes, net$colors)

    # Choose modules to test in order of size if not specified
    if (is.null(module_names)) {
        module_names <- net$colors %>% table %>% sort(decreasing=TRUE) %>% names
        module_names <- module_names[module_names != 'grey']
    }
    modules <- modules[module_names]

    # Run enrichment for all modules
    enrichment_allmodules <- BiocParallel::bplapply(
        seq_along(modules), function(x) {
            message("Enrichment analysis for module ",
                        names(modules)[x], "...")
            l <- par_enrich(
                genes = modules[[x]],
                reference = all_genes,
                genesets = annotation_list
            )

            # Filter for significant
            l <- l[l$padj < p, ]
            return(l)
    }, BPPARAM = bp_param)
    names(enrichment_allmodules) <- names(modules)
    
    # Remove modules with no enrichments
    enrichment_allmodules <- enrichment_allmodules[
        vapply(enrichment_allmodules, nrow, numeric(1)) > 0]

    # Combine into dataframe
    enrichment_df <- enrichment_allmodules %>% bind_rows(.id='module')
    return(enrichment_df)
}

# https://rdrr.io/github/almeidasilvaf/BioNERO/src/R/gcn_inference.R#sym-par_enrich
#' Helper function to perform Fisher's Exact Test with parallel computing
#'
#' @param genes Character vector containing genes on which enrichment will
#' be tested.
#' @param reference Character vector containing genes to be used as background.
#' @param genesets List of functional annotation categories.
#' (e.g., GO, pathway, etc.) with their associated genes.
#' @param adj Multiple testing correction method.
#' @param bp_param BiocParallel back-end to be used.
#' Default: BiocParallel::SerialParam()
#'
#' @return Results of Fisher's Exact Test in a data frame with TermID,
#' number of associated genes, number of genes in reference set,
#' P-value and adjusted P-value.
#'
#' @importFrom stats p.adjust fisher.test
#' @importFrom BiocParallel bplapply SerialParam
#' @noRd
par_enrich <- function(genes, reference, genesets, adj = "BH",
                       bp_param = BiocParallel::SerialParam()) {

    reference <- reference[!reference %in% genes]

    tab <- BiocParallel::bplapply(seq_along(genesets), function(i) {

        RinSet <- sum(reference %in% genesets[[i]])
        RninSet <- length(reference) - RinSet
        GinSet <- sum(genes %in% genesets[[i]])
        GninSet <- length(genes) - GinSet
        fmat <- matrix(
            c(GinSet, RinSet, GninSet, RninSet),
            nrow = 2, ncol = 2, byrow = FALSE
        )
        colnames(fmat) <- c("inSet", "ninSet")
        rownames(fmat) <- c("genes", "reference")
        fish <- stats::fisher.test(fmat, alternative = "greater")
        pval <- fish$p.value
        inSet <- RinSet + GinSet
        res <- c(GinSet, inSet, pval)
        return(res)
    }, BPPARAM = bp_param)

    rtab <- do.call(rbind, tab)
    rtab <- data.frame(as.vector(names(genesets)), rtab)
    rtab <- rtab[order(rtab[, 4]), ]
    colnames(rtab) <- c("TermID", "genes", "all", "pval")
    padj <- stats::p.adjust(rtab$pval, method = adj)
    tab.out <- data.frame(rtab, padj)

    return(tab.out)
}


count_enrichment_matches <- function(enrichments_pair) {
    A <- enrichments_pair[[1]]
    B <- enrichments_pair[[2]]
    modules <- unique(A$module)
    enrichment_list_A <- A$TermID %>% split(A$module)
    enrichment_list_B <- B$TermID %>% split(B$module)
    
    GO_A <- setNames(rep(NA, length(modules)), modules)
    GO_B <- setNames(rep(NA, length(modules)), modules)
    GO_match <- setNames(rep(NA, length(modules)), modules)
    GO_pct_match <- setNames(rep(NA, length(modules)), modules)
    for (module in modules) {
        if (exists(module, where=enrichment_list_B)) {
            matched <- sum(enrichment_list_A[[module]] %in% enrichment_list_B[[module]])
        } else {
            matched <- 0
        }
        
        GO_A[module] <- length(enrichment_list_A[[module]])
        GO_B[module] <- length(enrichment_list_B[[module]])
        GO_match[module] <- matched
        GO_pct_match[module] <- matched/GO_A[module]
    }
    return(data.frame(GO_A, GO_B, GO_match, GO_pct_match))
}



# View enrichments as table
pivot_enrichments_table <- function(enrichments) {
    l <- list()
    for (e in unique(enrichments$module)) {
        l[[e]] <- enrichments %>% filter(module==e) %>% .$TermID
    }
    df <- lapply(l, `length<-`, max(lengths(l))) %>% bind_cols
    return(df)
}




make_comparison_table <- function(source_net, target_net, n_modules=30, annotation=NULL) {
    # Match target net to source net
    target_net <- target_net %>% match_modules(source_net)
    print("Matched modules.")

    # Rename columns of count tables from source and target
    source_counts <- source_net$counts %>% 
        as_tibble %>% rename_with(~c('module', 'genes'))
    target_counts <- target_net$oldCounts %>% 
        as_tibble %>% rename_with(~c('module', 'genes')) %>% 
        mutate(new_name = target_net$counts %>% names, .before='genes')

    # Get match % of all genes in modules
    match_rate_all_genes <- count_matches(source_net, target_net) %>% 
        rownames_to_column('module')
    print("Got match % of all genes in modules.")

    # Get match % for top100 genes by kME for top 30 modules
    topN <- get_topN_matches(source_net, target_net, how='', n_modules=n_modules)
    match_rate_top100 <- topN %>% rownames_to_column('topN') %>% 
        pivot_longer(-topN, names_to='module', values_to='pct_kME100') %>% 
        filter(topN==100) %>% dplyr::select(-topN) %>% 
        mutate(module=str_replace(module, "kME","")) %>% 
        mutate(pct_kME100 = pct_kME100/100)
    print("Got match % of top 100 genes in modules.")

    summaryTable <- source_counts %>% 
        left_join(target_counts, by='module') %>% 
        rename_with(~ c('module', 'genes', 'module_B', 'genes_B')) %>% 
        left_join(match_rate_all_genes, by='module') %>% 
        left_join(match_rate_top100, by='module')

    # Get match % of GO enrichments
    if (!is.null(annotation)) {
        module_names <- source_net$colors %>% table %>% sort(decreasing=TRUE) %>% 
                    names %>% .[.!='grey'] %>% .[1:n_modules]
        enrichments_list <- list(source_net, target_net) %>% lapply(get_enrichments,
                annotation=annotation, module_names=module_names, bp_param=MulticoreParam(workers=8)
                )
        enrichment_matches <- count_enrichment_matches(enrichments_list) %>% 
            rownames_to_column('module')
        print("Got match % of module GO enrichments.")

        summaryTable <- summaryTable %>% 
        left_join(enrichment_matches, by='module')
    }
    
    return(summaryTable)
}