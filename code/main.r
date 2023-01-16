### Setup and load data
options(stringsAsFactors = FALSE)
options(max.print = 100)

setwd("code")

load("../data/rse_gene_unfiltered.Rdata")
load("../data/rse_tx_unfiltered.Rdata")

### Preprocess expression data
source("../code/preprocess.r")
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


### Fit WGCNA
# With/without qSVA
source("../code/wgcna.r")
sft <- get_soft_threshold(rse_noQSV)
power <- sft$powerEstimate

net_QSV <- fit_WGCNA(rse_QSV, power=3, threads = 9,
                loadTOM = TRUE, fileBase = "QSV")

net_noQSV <- fit_WGCNA(rse_noQSV, power=4, threads = 9,
                loadTOM = TRUE, fileBase = "noQSV")

plot_ngenes(net_noQSV)

# Splitting into random halves
source("../code/wgcna.r")
rse_noQSV_split <- rse_noQSV %>% split_samples
split_nets_noQSV <- make_nets(rse_noQSV_split)
saveRDS(split_nets_noQSV, "../outputs/split_nets_noQSV.rds")

split_nets_noQSV <- readRDS("../outputs/split_nets_noQSV.rds")


### Analyse stability
topN <- get_topN_matches(split_nets_noQSV[[1]], split_nets_noQSV[[2]], how='', n_modules=10)

source("../code/plots_wgcna.r")
plot_topN_lines(topN) +
ggtitle('% of top N genes by kME retained in matched module, random splits')
# ggtitle('% of top N genes by kME retained in *top N genes of* matched module, random splits')


split_nets_noQSV[[1]] %>% plot_kME_topN()

topN_QSV <- get_topN_matches(net_QSV, net_noQSV)
topN_QSV %>% plot_topN_lines(titlesub='QSV to no-QSV')

source("../code/plots_wgcna.r")
overlap <- get_overlap(split_nets_noQSV)
plot_overlap_counts(split_nets_noQSV)


### Compare enrichments
source("../code/analyse_coexp.r")
GO <- get_GO_annotations()
split_nets_noQSV_matched <- list(
    split_nets_noQSV[[1]],
    split_nets_noQSV[[2]] %>% match_modules(split_nets_noQSV[[1]])
)
# Check matched modules from each split
module_names <- split_nets_noQSV[[1]]$colors %>% table %>% sort(decreasing=TRUE) %>% 
                names %>% .[.!='grey'] %>% .[1:11]
enrichments_split_noQSV <- split_nets_noQSV_matched %>% lapply(get_enrichments,
            annotation=GO, module_names=module_names
            # , bp_param=MulticoreParam(workers=8)
            )

enrichments_split_noQSV[[2]] %>% pivot_enrichments_table %>% write_csv("../outputs/enrichments_noQSV_split2.csv")


net <- split_nets_noQSV_matched[[2]]
annotation <- GO
test <- par_enrich(
                genes = modules[['steelblue']],
                reference = all_genes,
                genesets = annotation_list
            )

count_enrichment_matches(enrichments_split_noQSV)


# Make table of comparisons
source("../code/wgcna.r")
# module A | n_genes | module B | n_genes | matched | pct_matched | top 100 kME pct_matched | enrichments A | enrichments B | pct enrichments
split_nets_noQSV_matched <- list(
    split_nets_noQSV[[1]],
    split_nets_noQSV[[2]] %>% match_modules(split_nets_noQSV[[1]])
)
A <- split_nets_noQSV_matched[[1]]$counts %>% as_tibble %>% rename_with(~c('module', 'genes'))
B <- split_nets_noQSV_matched[[2]]$oldCounts %>% as_tibble %>% rename_with(~c('module', 'genes')) %>% 
    mutate(new_name = split_nets_noQSV_matched[[2]]$counts %>% names, .before='genes')

topN_100 <- topN %>% rownames_to_column('topN') %>% 
    pivot_longer(-topN, names_to='module', values_to='pct_kME100') %>% 
    filter(topN==100) %>% dplyr::select(-topN) %>% 
    mutate(module=str_replace(module, "kME","")) %>% 
    mutate(pct_kME100 = pct_kME100/100)

A %>% 
    left_join(B, by='module') %>% rename_with(~ c('module', 'genes', 'module_B', 'genes_B')) %>% 
    left_join(count_matches(split_nets_noQSV_matched) %>% rownames_to_column('module'), by='module') %>% 
    left_join(topN_100, by='module') %>% 
    left_join(count_enrichment_matches(enrichments_split_noQSV) %>% rownames_to_column('module'), by='module') %>% 
    write_csv("../outputs/table_noQSV_splits.csv")

### Gene-qSV correlations

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
