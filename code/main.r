### Setup and load data
options(stringsAsFactors = FALSE)
options(max.print = 100)
# enableWGCNAThreads(nThreads = 9)
cor <- WGCNA::cor

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
saveRDS(rse_norm_samples, "../data/rse_norm_samples.rds")

formula <- "~ Dx + Age + Sex + Race + mitoRate + rRNA_rate + totalAssignedGene + RIN +
                         ast + end + mic + neu + oli + opc
                         "
rse_regressed <- rse_norm_samples %>% 
    regress_covariates(rse_tx, formula = formula, n_qsvs = 0, P=4)
rse_norm_genes <- rse_regressed %>% normalise_genes

saveRDS(rse_norm_genes, "../outputs/rse_processed_noQSV.rds")

rse_QSV <- readRDS("../outputs/rse_processed_QSV.rds")
rse_noQSV <- readRDS("../outputs/rse_processed_noQSV.rds")


### Fit WGCNA
# With/without qSVA
source("../code/wgcna.r")
sft <- get_soft_threshold(rse_QSV)
power <- sft$powerEstimate

net_QSV <- fit_WGCNA(rse_QSV, power=power, threads = 9)
                # loadTOM = FALSE, fileBase = "QSV")

net_noQSV <- fit_WGCNA(rse_noQSV, power=4, threads = 9,
                loadTOM = TRUE, fileBase = "noQSV")

plot_ngenes(net_noQSV)

# Splitting into random halves
source("../code/wgcna.r")
rse_QSV_split <- rse_QSV %>% split_samples
split_nets_QSV <- make_nets(rse_QSV_split)
saveRDS(split_nets_QSV, "../outputs/split_nets_QSV.rds")
print("Done QSV splits.\n")

rse_noQSV_split <- rse_noQSV %>% split_samples
split_nets_noQSV <- make_nets(rse_noQSV_split)
saveRDS(split_nets_noQSV, "../outputs/split_nets_noQSV.rds")
print("Done noQSV splits.\n")

QSV_noGrey_genes <- names(net_QSV$colors[net_QSV$colors!='grey'])
rse_noQSVnoGrey <- rse_noQSV[QSV_noGrey_genes,]
rse_noQSVnoGrey_split <- rse_noQSVnoGrey %>% split_samples
split_nets_noQSVnoGrey <- make_nets(rse_noQSVnoGrey_split)
saveRDS(split_nets_noQSVnoGrey, "../outputs/split_nets_noQSVnoGrey.rds")
print("Done noQSVnoGrey splits.\n")



### Module preservation statistics
multiExpr_QSV <- list(
    A = list(data=t(assays(rse_QSV_split[[1]])[['ranknorm']])),
    B = list(data=t(assays(rse_QSV_split[[2]])[['ranknorm']]))
)
multiExpr_noQSV <- list(
    A = list(data=t(assays(rse_noQSV_split[[1]])[['ranknorm']])),
    B = list(data=t(assays(rse_noQSV_split[[2]])[['ranknorm']]))
)

multiColor_QSV <- list(
    A = split_nets_QSV[[1]]$colors,
    B = split_nets_QSV[[2]]$colors
)
multiColor_noQSV <- list(
    A = split_nets_noQSV[[1]]$colors,
    B = split_nets_noQSV[[2]]$colors
)

mp_noQSV = modulePreservation(
    multiData = multiExpr_noQSV,
    multiColor = multiColor_noQSV,
    networkType = 'signed hybrid',
    nPermutations = 1000,
    randomSeed = 1,
    parallelCalculation = FALSE,
    verbose = 3
)
save(mp_noQSV, file='../outputs/modulePreservation_noQSV.RData')
mp_QSV = modulePreservation(
    multiData = multiExpr_QSV,
    multiColor = multiColor_QSV,
    networkType = 'signed hybrid',
    nPermutations = 1000,
    randomSeed = 1,
    parallelCalculation = FALSE,
    verbose = 3
)
save(mp_QSV, file='../outputs/modulePreservation_QSV.RData')


mp_noQSV$preservation$Z[[1]][[2]][WGCNA::standardColors(10),c('Zdensity.pres','Zconnectivity.pres')]
mp_noQSV$accuracy$Z[[1]][[2]][WGCNA::standardColors(10),c('moduleSize','Z.accuracy')]
mp_QSV$preservation$Z[[1]][[2]][WGCNA::standardColors(10),c('Zdensity.pres','Zconnectivity.pres')]
mp_QSV$accuracy$Z[[1]][[2]][WGCNA::standardColors(10),c('moduleSize','Z.accuracy')]



### Analyse stability
source("../code/analyse_coexp.r")
topN_QSV <- get_topN_matches(split_nets_QSV[[1]], split_nets_QSV[[2]], how='', n_modules=30)

source("../code/plots_wgcna.r")
plot_topN_lines(topN_QSV) +
ggtitle('% of top N genes by kME retained in matched module, random splits')
# ggtitle('% of top N genes by kME retained in *top N genes of* matched module, random splits')


# split_nets_noQSV[[1]] %>% plot_kME_topN()

# topN_QSV <- get_topN_matches(net_QSV, net_noQSV)
# topN_QSV %>% plot_topN_lines(titlesub='QSV to no-QSV')

# source("../code/plots_wgcna.r")
# overlap <- get_overlap(split_nets_noQSV)
# plot_overlap_counts(split_nets_noQSV)



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


split_nets_QSV_matched <- list(
    split_nets_QSV[[1]],
    split_nets_QSV[[2]] %>% match_modules(split_nets_QSV[[1]])
)
# Check matched modules from each split
module_names <- split_nets_QSV[[1]]$colors %>% table %>% sort(decreasing=TRUE) %>% 
                names %>% .[.!='grey'] %>% .[1:11]
enrichments_split_QSV <- split_nets_QSV_matched %>% lapply(get_enrichments,
            annotation=GO, module_names=module_names
            # , bp_param=MulticoreParam(workers=8)
            )


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
source("../code/analyse_coexp.r")
# module A | n_genes | module B | n_genes | matched | pct_matched | top 100 kME pct_matched | enrichments A | enrichments B | pct enrichments

source_net <- split_nets_QSV[[1]]
target_net <- split_nets_QSV[[2]]

GO <- get_GO_annotations()


split_nets_QSV <- readRDS("../outputs/split_nets_QSV.rds")
split_nets_noQSV <- readRDS("../outputs/split_nets_noQSV.rds")
split_nets_noQSVnoGrey <- readRDS("../outputs/split_nets_noQSVnoGrey.rds")



comparisonTable_QSV <- make_comparison_table(split_nets_QSV[[1]], split_nets_QSV[[2]], annotation=GO)
comparisonTable_QSV %>% write_csv("../outputs/comparisonTable_QSV.csv")

comparisonTable_noQSV <- make_comparison_table(split_nets_noQSV[[1]], split_nets_noQSV[[2]], annotation=GO)
comparisonTable_noQSV %>% write_csv("../outputs/comparisonTable_noQSV.csv")

comparisonTable_noQSVnoGrey <- make_comparison_table(
    split_nets_noQSVnoGrey[[1]], split_nets_noQSVnoGrey[[2]], n_modules=20, annotation=GO)
comparisonTable_noQSVnoGrey %>% write_csv("../outputs/comparisonTable_noQSVnoGrey.csv")




# split_nets_noQSV_matched <- list(
#     split_nets_noQSV[[1]],
#     split_nets_noQSV[[2]] %>% match_modules(split_nets_noQSV[[1]])
# )
# A <- split_nets_noQSV_matched[[1]]$counts %>% as_tibble %>% rename_with(~c('module', 'genes'))
# B <- split_nets_noQSV_matched[[2]]$oldCounts %>% as_tibble %>% rename_with(~c('module', 'genes')) %>% 
#     mutate(new_name = split_nets_noQSV_matched[[2]]$counts %>% names, .before='genes')

# topN_100 <- topN %>% rownames_to_column('topN') %>% 
#     pivot_longer(-topN, names_to='module', values_to='pct_kME100') %>% 
#     filter(topN==100) %>% dplyr::select(-topN) %>% 
#     mutate(module=str_replace(module, "kME","")) %>% 
#     mutate(pct_kME100 = pct_kME100/100)

# A %>% 
#     left_join(B, by='module') %>% rename_with(~ c('module', 'genes', 'module_B', 'genes_B')) %>% 
#     left_join(count_matches(split_nets_noQSV_matched) %>% rownames_to_column('module'), by='module') %>% 
#     left_join(topN_100, by='module') %>% 
#     left_join(count_enrichment_matches(enrichments_split_noQSV) %>% rownames_to_column('module'), by='module') %>% 
#     write_csv("../outputs/table_noQSV_splits.csv")



# split_nets_QSV_matched <- list(
#     split_nets_QSV[[1]],
#     split_nets_QSV[[2]] %>% match_modules(split_nets_QSV[[1]])
# )
# A <- split_nets_QSV_matched[[1]]$counts %>% as_tibble %>% rename_with(~c('module', 'genes'))
# B <- split_nets_QSV_matched[[2]]$oldCounts %>% as_tibble %>% rename_with(~c('module', 'genes')) %>% 
#     mutate(new_name = split_nets_QSV_matched[[2]]$counts %>% names, .before='genes')

# topN_QSV <- get_topN_matches(split_nets_QSV[[1]], split_nets_QSV[[2]], how='topN', n_modules=20)
# plot_topN_lines(topN_QSV) +
# ggtitle('% of top N genes by kME retained in *top N genes of* matched module, random splits, QSV')

# topN_100 <- topN_QSV %>% rownames_to_column('topN') %>% 
#     pivot_longer(-topN, names_to='module', values_to='pct_kME100') %>% 
#     filter(topN==100) %>% dplyr::select(-topN) %>% 
#     mutate(module=str_replace(module, "kME","")) %>% 
#     mutate(pct_kME100 = pct_kME100/100)

# A %>% 
#     left_join(B, by='module') %>% rename_with(~ c('module', 'genes', 'module_B', 'genes_B')) %>% 
#     left_join(count_matches(split_nets_QSV_matched) %>% rownames_to_column('module'), by='module') %>% 
#     left_join(topN_100, by='module') %>% 
#     left_join(count_enrichment_matches(enrichments_split_QSV) %>% rownames_to_column('module'), by='module') %>% 
#     write_csv("../outputs/table_QSV_splits.csv")









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


### Gene qSV correlations for curiosity

## Make qSVs to correlate with genes
rse_tx_matched <- match_rse_tx_samples(rse_tx, rse_norm_samples)
formula <- "~ Age + Sex + Race + mitoRate + rRNA_rate + totalAssignedGene + RIN +
                         ast + end + mic + neu + oli + opc"
qsvs <- make_qsvs(rse_tx_matched, formula)
gene_qsvs_correlations <- assays(rse_norm_samples)[['logRPKM_quant_normalised']] %>% t %>% cor(qsvs)
gene_qsvs_correlations %>% data.frame %>% rownames_to_column('ensembl_id') %>% 
write_csv("../outputs/gene_qsv_correlations.csv")
# plot gene qsvs corrs
gene_qsvs_correlations %>% data.frame %>% 
pivot_longer(everything(), names_to='qSV', values_to='r') %>% 
mutate(qSV = factor(qSV, ordered=T, levels=unique(.$qSV))) %>% 
ggplot(aes(x=r, color=qSV)) + geom_density() + theme_minimal() + 
xlim(-1,1) + ggtitle('Distributions of gene-qSV correlations')
# plot gene qsvs max corr
gene_qsvs_correlations %>% 
data.frame %>% 
mutate(gene=row_number()) %>% 
pivot_longer(-gene, names_to='qSV', values_to='r') %>%  
group_by(gene) %>% summarise(max_abs=max(abs(r))) %>% 
ggplot(aes(x=max_abs)) + 
stat_ecdf(geom = "step", color='blue') +
theme_minimal() + 
scale_x_continuous(breaks=seq(0,1,.2), limits=c(0,1), name='max absolute correlation of gene to any qSV') + 
ylab('% of genes') +
ggtitle('Cumulative distribution of max absolute gene-qSV correlation')
# plot gene qsvs corrs heatmap
library(pals)
gene_qsvs_correlations %>% data.frame %>% 
mutate_all(function(x) ifelse(abs(x)<.5,0,x)) %>% 
cor %>%
data.frame %>% rownames_to_column('qSV_x') %>% 
pivot_longer(-qSV_x, names_to='qSV_y', values_to='r') %>% 
mutate(
    qSV_x = factor(qSV_x, ordered=T, levels=unique(.$qSV_x)),
    qSV_y = factor(qSV_y, ordered=T, levels=unique(.$qSV_y))
) %>% 
ggplot(aes(x=qSV_x, y=qSV_y, fill=r)) + geom_raster() + 
scale_fill_gradientn(colors=rev(brewer.rdbu(200)), limits=c(-1,1)) + 
theme_classic() + theme(aspect.ratio=1) + xlab('') + ylab('') +
ggtitle('qSV-qSV correlations of gene-qSV correlations, filtered for ')