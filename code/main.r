options(stringsAsFactors = FALSE)
options(max.print = 100)

# setwd("code")

load("../data/rse_gene_unfiltered.Rdata")
load("../data/rse_tx_unfiltered.Rdata")

source("preprocess.r")
rse_merge_vars <- rse_gene %>% merge_covariates
rse_select_samples <- rse_merge_vars %>% select_samples
rse_clean_samples <- rse_select_samples %>% clean_outlier_samples
rse_clean_genes <- rse_clean_samples %>% clean_genes
rse_cell_proportions <- rse_clean_genes %>% get_cell_proportions
rse_norm_samples <- rse_cell_proportions %>% normalise_samples
saveRDS(rse_norm_samples, "data/rse_norm_samples.rds")

formula <- "~ Age + Sex + mitoRate + rRNA_rate + totalAssignedGene + RIN +
                         snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + 
                         snpPC6 + snpPC7 + snpPC8 + snpPC9 + snpPC10"
rse_regressed <- rse_norm_samples %>% 
    regress_covariates(rse_tx, formula = formula, n_qsvs = NULL)
rse_norm_genes <- rse_regressed %>% normalise_genes

saveRDS(rse_norm_genes, "../outputs/rse_processed_QSV.rds")

rse_QSV <- readRDS("../outputs/rse_processed_QSV.rds")
rse_noQSV <- readRDS("../outputs/rse_processed_noQSV.rds")

source("../code/wgcna.r")
sft <- get_soft_threshold(rse_noQSV)
power <- sft$powerEstimate

net_QSV <- fit_WGCNA(rse_QSV, power=3, threads = 9,
                loadTOM = TRUE, fileBase = "QSV")

net_noQSV <- fit_WGCNA(rse_noQSV, power=4, threads = 9,
                loadTOM = TRUE, fileBase = "noQSV")

plot_ngenes(net_noQSV)



source("../code/wgcna.r")

source_weights <- data.frame(kME_noQSV)[,"kMEturquoise"]
target_weights <- data.frame(kME_QSV)[,"kMEturquoise"]

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


source("../code/wgcna.r")
rse_noQSV_split <- rse_noQSV %>% split_samples
split_nets_noQSV <- make_nets(rse_noQSV_split)
saveRDS(split_nets_noQSV, "../outputs/split_nets_noQSV.rds")

topN <- get_topN_matches(split_nets_noQSV[[1]], split_nets_noQSV[[2]], how='')

source("../code/plots_wgcna.r")
plot_topN_lines(topN) +
ggtitle('% of top N genes by kME retained in matched module, random splits')
# ggtitle('% of top N genes by kME retained in *top N genes of* matched module, random splits')



split_nets_noQSV[[1]] %>% plot_kME_topN()



topN_QSV <- get_topN_matches(net_QSV, net_noQSV)
topN_QSV %>% plot_topN_lines(titlesub='QSV to no-QSV')

source_net <- net_noQSV
target_net <- net_QSV


net_QSV_matched <- net_QSV %>% match_modules(net_noQSV)
kME_noQSV <- net_noQSV %>% filter_kME()
kME_QSV <- net_QSV_matched %>% filter_kME()
matched_order <- colnames(kME_noQSV)[colnames(kME_noQSV) %in% colnames(kME_QSV)]
kME_QSV <- kME_QSV[, matched_order]

mapply(get_topN_matches_one_module, data.frame(kME_noQSV[,2:3]), data.frame(kME_QSV[,2:3]))



get_topN_matches_one_module(data.frame(kME_noQSV)[,"kMEturquoise"], data.frame(kME_QSV)[,"kM)
get_topN_matches_one_module(kME_noQSV[,5:5], kME_QSV[,5:5])



# Correlate all genes from noQSV data with QSVs
source("qsva.r")
rse_tx_matched <- match_rse_tx_samples(rse_tx, rse_select_samples)
# QSVs in columns, samples in rows
qsvs <- make_qsvs(rse_tx_matched, formula=formula)
# genes in columns, samples in rows
genes <- assays(rse_noQSV)$ranknorm %>% t
# correlate
gene_qsv_cor <- cor(genes, qsvs) %>% 
    data.frame %>% rownames_to_column('gene')
gene_qsv_cor

# Check which genes are in which module with/out QSV
gene_modules <- net_noQSV$colors %>% enframe %>% 
left_join(net_QSV$colors %>% enframe, by='name') %>% 
rename_all(~ c('gene','module_noQSV', 'module_QSV')) %>% 
mutate(
    module_noQSV_grey = module_noQSV=='grey',
    module_QSV_grey = module_QSV=='grey',
    module_category = case_when(
        module_noQSV_grey & module_QSV_grey ~ 'always grey',
        !module_noQSV_grey & !module_QSV_grey ~ 'never grey',
        !module_noQSV_grey & module_QSV_grey ~ 'only grey with qSVs',
        module_noQSV_grey & !module_QSV_grey ~ 'only grey without qSVs'
        )
)

# Join to correlations
gene_modules_cor <- gene_modules %>%
    left_join(gene_qsv_cor, by='gene') %>% 
    pivot_longer(cols=qSV1:qSV33, names_to='qSV', values_to='r') %>% 
    mutate(qSV = parse_number(qSV))
gene_modules_cor %>% 
    filter(qSV < 10) %>% 
    ggplot(aes(x=factor(qSV), y=abs(r), fill=module_category)) + 
    geom_boxplot(position='dodge') +
    xlab('qSV') + ylab('gene-qSV correlation (absolute)') +
    theme_classic()
