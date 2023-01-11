# library(usethis)
# usethis::edit_r_environ()

library(WGCNA)
library(BioNERO)
library(qsvaR)
library(limma)
library(biomaRt)
library(pals)
library(reshape2)
library(tidyverse)

# library(BiocFileCache)
# library(gsubfn)
#library(vscDebugger)
options(stringsAsFactors = FALSE)
options(max.print = 100)
# enableWGCNAThreads()

trace(exp2gcn, edit = TRUE)

source("qsva.r")
source("wgcna.r")
source("plots_wgcna.r")

load("data/rse_gene_BrainSeq_phases_1_2.n1627.Rdata")
# Set rownames to ensemblID (ie remove gene version)
#rownames(rse_gene) <- rowData(rse_gene)$ensemblID

rse_gene_controls <- rse_gene[,
  rse_gene$Dx == "Control" & rse_gene$Dataset == "BrainSeq_Phase2_DLPFC"]

rse_gene_controls_clean <- rse_gene_controls %>%
  .[rowRanges(rse_gene)$gene_type == "protein_coding", ] %>%
  exp_preprocess(remove_nonexpressed = TRUE,
                 min_exp = 10,
                 variance_filter = FALSE,
                 n = 2000,
                 Zk_filtering = FALSE,
                 cor_method = "pearson",
                 remove_confounders = FALSE,
                 vstransform = TRUE
                 )

rse_tx_controls <- get_rse_tx()
qsvs <- make_qsvs(rse_tx_controls)

rse_gene_controls_qsv <- rse_gene_controls_clean %>% 
  regress_qsvs(rse_tx_controls)

source('qsva.r')
rse_gene_controls_clean_methods <- list(
  # "limma" = regress_qsvs(rse_gene_controls_clean, qsvs, method='limma'),
  # "lm" = regress_qsvs(rse_gene_controls_clean, qsvs, method='lm'),
  "libd" = regress_qsvs(rse_gene_controls_clean, qsvs, method='libd'),
  "libd1" = regress_qsvs(rse_gene_controls_clean, qsvs, method='libd', P=1)
)
rse_gene_controls_clean_methods <- list(
  # "limma_6qsvs" = regress_qsvs(rse_gene_controls_clean, qsvs, n_qsvs=6, method='limma'),
  # "lm_6qsvs" = regress_qsvs(rse_gene_controls_clean, qsvs, n_qsvs=6, method='lm'),
  "libd_6qsvs" = regress_qsvs(rse_gene_controls_clean, qsvs, n_qsvs=6, method='libd'),
  "libd1_6qsvs" = regress_qsvs(rse_gene_controls_clean, qsvs, n_qsvs=6, method='libd', P=1)
) %>% c(rse_gene_controls_clean_methods)

source('wgcna.r')
methods_nets <- rse_gene_controls_clean_methods %>% make_nets(assayname='cleaned')
lapply(methods_nets, plot_ngenes_per_module)

plot_ngenes_per_module(methods_nets[['limma']])
plot_ngenes_per_module(methods_nets[['lm']])
plot_ngenes_per_module(methods_nets[['libd']])



rse_gene_controls_qsv_libd <- regress_qsvs_libd(rse_gene_controls_clean, rse_tx_controls)
qsv_nets['18 qsvs libd'] <- make_nets(list(rse_gene_controls_qsv_libd))




mod_gene <- model.matrix(~ Age + Sex + Race,
    data = colData(rse_gene_controls_clean)
)
# Add qSVs to design matrix
n_qsvs <- 1
if (!is.null(n_qsvs)) {
    qsvs <- qsvs[,1:n_qsvs]
}
mod <- cbind(mod_gene, qsvs)

exp <- assay(rse_gene_controls_clean) %>% as.matrix
library(jaffelab)
cleaned_exp <- cleaningY(exp, mod, P = 0)
assays(rse)$svacleaned = cleaned_matrix

lmFit(exp, mod) %>% eBayes %>% residuals(exp) %>% .[1:5,1:5]
lmFit(exp, mod) %>% residuals(exp) %>% .[1:5,1:5]

fit <- lm(t(exp) ~ mod)
fit %>% str
fit <- lmFit(exp, mod)

Hat <- solve(t(mod) %*% mod) %*% t(mod)
y <- exp
## For dealing with NAs
## https://stackoverflow.com/questions/16535084/matrix-multiplication-with-scattered-na-values
ty <- t(y)
ty[is.na(ty)] <- 0
beta <- (Hat %*% ty)
## Note that y might still have the NAs, and NA - a number = NA
## so there's no need to reset the NAs back on cleany
cleany <- y - t(as.matrix(mod[, -c(seq_len(P))]) %*% beta[-seq_len(P), ])
return(cleany)


rse_gene_controls_qsv1 <- rse_gene_controls_clean %>% 
  regress_qsvs(rse_tx_controls, n_qsvs=1)
rse_gene_controls_qsv6 <- rse_gene_controls_clean %>% 
  regress_qsvs(rse_tx_controls, n_qsvs=6)

qsv_nets <- list(
  '0 qsvs'=rse_gene_controls_clean,
  '1 qsv'=rse_gene_controls_qsv1,
  '6 qsvs'=rse_gene_controls_qsv6,
  '18 qsvs'=rse_gene_controls_qsv
) %>% make_nets

net_qsv18 <- list(rse_gene_controls_qsv) %>% make_nets
plot_ngenes_per_module(net_qsv18[[1]])
plot_dendro_and_colors(net_qsv18[[1]])

plot_ngenes_per_module(qsv_nets[['18 qsvs libd']])
plot_dendro_and_colors(qsv_nets[['18 qsvs libd']])




load("data/rse_gene_unfiltered.Rdata")
load("data/rse_tx_unfiltered.Rdata")

degradation_transcripts <- select_transcripts('cell_component')
sum(degradation_transcripts %in% rownames(rse_tx)) / length(degradation_transcripts)

sum(degradation_transcripts %in% rownames(rse_gene)) / length(degradation_transcripts)

# rse_gene_controls <- rse_gene[, rse_gene$Dx == "Control" & rse_tx$Region=="DLPFC"]
rse_tx_controls <- rse_tx[, rse_tx$Dx == "Control" & rse_tx$Region=="DLPFC"]

### qSVs using genes only
mod <- model.matrix(~ Age + Sex + Race,
    data = colData(rse_tx_controls)
)
# qSVs of matrix
qsvs <- qSVA(rse_tx = rse_tx_controls,
              type = "cell_component",
              mod = mod,
              assayname = "tpm"
)




# net_splits <- rse_gene_controls_qsv %>% fit_net_splits(n_splits = 10)
# save(net_splits, file = "spearman_splits_10.RData")

# load("pearson_splits_10.RData")

net_splits_overlap_pct <- net_splits %>% lapply(get_overlap_pct) %>% bind_rows(.id='split')
net_splits_overlap_pct %>% write_csv('pearson10_overlap_pct.csv')
source('plots_wgcna.r')
plot_overlap_pct(net_splits_overlap_pct, min_points = 8)

plot_overlap_counts(net_splits[[7]])

plot_dendro_and_colors(net_splits[[1]][[1]])
plot_ngenes_per_module(net_splits[[1]][[1]])


exp_base <- assay(rse_gene_controls_clean)
qsvs <- make_qsvs(rse_tx_controls)

gene_qSV_cors <- exp_base %>% as.matrix %>% melt %>% 
  rename(gene=Var1, sample=Var2, exp=value) %>% 
  mutate(sample=str_replace(sample, "_.*", "")) %>% 
  left_join(qsvs %>% data.frame %>% rownames_to_column('sample'), by='sample') %>% 
  pivot_longer(cols=qSV1:qSV18, names_to='qSV', values_to='qSV_value') %>% 
  group_by(gene, qSV) %>%
  summarise(
    r=cor.test(exp, qSV_value)$estimate,
    p=cor.test(exp, qSV_value)$p.value
  ) %>% 
  mutate(star = case_when(p < 0.001 ~ '***',p < 0.01 ~ '**',p < 0.05 ~ '*', T ~ ''))

net <- net_pair[[2]]
colours <- net$genes_and_modules %>% group_by(Modules) %>% summarise(n=n()) %>% arrange(-n) %>% .$Modules

gene_qSV_cors %>% 
  left_join(net$genes_and_modules %>% rename(gene=Genes, module=Modules), by='gene') %>%
  mutate(module = factor(module, ordered=T, levels=colours)) %>% 
  mutate(qSV = factor(qSV, ordered=T, levels=paste0('qSV',seq(1:18)))) %>% 
  filter(module %in% c('grey', 'turquoise', 'blue', 'brown', 'yellow', 'green', 'red')) %>% 
  ggplot(aes(x=qSV, y=r, group=interaction(qSV,module))) + 
  geom_boxplot(aes(fill=module), alpha=.5, position=position_dodge(), width=.7) +
  scale_fill_identity() +
  ylab('correlation') +
  scale_x_discrete('qSV') +
  theme_classic() + 
  theme(axis.text.x=element_text(angle=0, hjust=1)) +
  ggtitle('Pre-qSV gene expression vs qSV, distribution of correlations grouped by post-qSV module')

gene_qSV_cors %>% 
  left_join(net$genes_and_modules %>% rename(gene=Genes, module=Modules), by='gene') %>%
  mutate(module = factor(module, ordered=T, levels=colours)) %>% 
  mutate(qSV = factor(qSV, ordered=T, levels=paste0('qSV',seq(1:18)))) %>% 
  group_by(qSV, module) %>% summarise(r=mean(abs(r))) %>% 
  ggplot(aes(x=qSV, y=module)) + 
  geom_tile(aes(fill=r)) +
  # geom_text(aes(label=star)) +
  scale_fill_gradientn(colors=brewer.rdbu(100), limits=c(-1,1)) +
  scale_y_discrete(name='module', limits=rev) +
  xlab('qSV') +
  coord_fixed(ratio=.3) +
  theme_classic() + 
  theme(
    axis.text.y=element_text(angle=0, hjust=1, colour=colours %>% rev),
    axis.text.x=element_text(angle=30, hjust=1)
  ) +
  ggtitle('Pre-qSV gene expression vs qSV, absolute correlation, mean by post-qSV module')


  mutate(Modules = factor(Modules, ordered=T, levels=unique(.$Modules))) %>% 
  ggplot(aes(x=Modules, y=log10(variances), color=Modules)) + 
  geom_jitter(size=.1, alpha=.5) + 
  geom_boxplot(width=.3, outlier.shape=NA) +
  scale_color_identity() +
  theme_classic() + theme(axis.text.x=element_text(angle=30, hjust=1))

rse_split <- rse_gene_controls_qsv %>% split_samples
net_pair <- rse_split %>% make_nets

exp <- rse_split[[2]] %>% assay


t(exp) %>% cbind(qsvs)

vars_ = exp %>% as.matrix %>% rowVars
df_module_vars = net$genes_and_modules %>% 
  mutate(variances=vars_) %>% 
  group_by(Modules) %>% mutate(n = n()) %>% arrange(-n)
df_module_vars %>% 
  mutate(Modules = factor(Modules, ordered=T, levels=unique(.$Modules))) %>% 
  ggplot(aes(x=Modules, y=log10(variances), color=Modules)) + 
  geom_jitter(size=.1, alpha=.5) + 
  geom_boxplot(width=.3, outlier.shape=NA) +
  scale_color_identity() +
  theme_classic() + theme(axis.text.x=element_text(angle=30, hjust=1))


plot_overlap_counts(net_pair)
plot_dendro_and_colors(net_pair[[2]])
plot_ngenes_per_module(net_pair[[2]])

# load("test_pair.RData")
# save(net_pair, file = "test_pair.RData")


overlap <- get_overlap(net_splits[[1]])

melt_table <- overlap$countTable %>% melt

melt_table %>% group_by(Var1) %>% summarise(total=sum(value)) %>% 
  mutate(Var2 = Var1) %>% 
  left_join(melt_table, by = c('Var1', 'Var2')) %>% select(-Var2) %>% 
  mutate(pct=value/total)

data.frame(total=totals) %>% left_join(match=matches)
totals <- rowSums(overlap$countTable)
matches <- diag(overlap$countTable)
totals[as.numeric(names(matches))] <- matches
totals

merge(
  rowSums(overlap$countTable), 
  diag(overlap$countTable))

plot_overlap_pct(net_pair) %>% data.frame

modules1 <- net_pair[[1]]$genes_and_modules$Modules
modules2 <- net_pair[[2]]$genes_and_modules$Modules
relabel=T
if (relabel) {
    modules2 <- WGCNA::matchLabels(modules2, modules1)
} 
overlap <- WGCNA::overlapTable(modules1, modules2)

all_colours <- WGCNA::standardColors(200)
colours1 <- c(all_colours[all_colours %in% unique(modules1)], 'grey')
colours2 <- c(all_colours[all_colours %in% unique(modules2)], 'grey')
overlap$countTable %>% rownames
overlap$countTable[colours1, colours2]

n_modules <- dim(overlap$countTable)
colours <- list(
    c(WGCNA::standardColors(n_modules[1]-1), 'grey'),
    c(WGCNA::standardColors(n_modules[2]-1), 'grey')
)


overlap$countTable %>% melt() %>% 
left_join(overlap$pTable %>% melt, by=c("Var1", "Var2")) %>% 
rename(net1=Var1, net2=Var2, n=value.x, p=value.y) %>% 
mutate(net1=factor(net1, ordered=T, levels=WGCNA::standardColors(50))) %>% 
mutate(net2=factor(net2, ordered=T, levels=WGCNA::standardColors(50))) %>% 
ggplot(aes(x=net2, y=net1, fill=n)) + 
geom_tile() + 
scale_fill_gradientn(colours=brewer.blues(100)) +
scale_y_discrete(limits=rev) +
#geom_text(aes(label=paste0(value.x, '\n', round(value.y,3)))) +
geom_text(aes(label=n, color=ifelse(n>2000,'white','black'))) +
scale_color_identity(guide='none') + 
theme_classic() + 
theme(
  axis.text.x=element_text(angle=30,hjust=1,color=WGCNA::standardColors(50))
)

colours1 <- net_pair[[1]]$genes_and_modules %>% 
group_by(Modules) %>% summarise(n=n()) %>% arrange(-n)

colours2 <- net_pair[[2]]$genes_and_modules %>% 
group_by(Modules) %>% summarise(n=n()) %>% arrange(-n)



merge(
  net_pair[[1]]$genes_and_modules,
  net_pair[[2]]$genes_and_modules,
  by='Genes'
) %>% group_by(Modules.x, Modules.y) %>% 
summarise(n=n())

cbind(net_pair[[1]]$genes_and_modules, net_pair[[1]]$kME) %>% 
  pivot_longer(!c(Genes, Modules), names_to = 'Module', values_to = 'kME') %>% 
  mutate(kME = ifelse(Module==Modules, kME, NA)) %>% 
  pivot_wider(id_cols=c(Genes, Modules), 
      names_from = 'Module', values_from = 'kME') %>% head

a <- net_pair[[1]]$kME
b <- net_pair[[2]]$kME
r <- cor(a,b)
r %>% melt %>% 
  ggplot(aes(Var1, Var2, fill=value)) + 
  geom_tile() + scale_fill_gradientn(colours=rev(brewer.rdbu(100))) +
  theme_classic() + theme(axis.text.x=element_text(angle=30, hjust=1))

get_hubs_gcn(rse_split[[1]] %>% assay, net) %>% head





    #exp <- handleSE(exp)
    exp <- rse_split[[1]] %>% assay
    genes <- rownames(rse_split[[1]])
    net <- net_pair[[1]]
    cor_method <- net$params$cor_method
    genes_modules <- net$genes_and_modules
    MEs <- net$MEs
    kIN <- net$kIN[, 2, drop=FALSE]

    # Add kIN info
    genes_modulesk <- merge(
        genes_modules, kIN, by.x = "Genes", by.y = "row.names"
    )

    # Calculate kME
    if(cor_method == "spearman") {
        MM <- WGCNA::signedKME(
            t(exp), MEs, outputColumnName = "",
            corOptions = "use = 'p', method = 'spearman'"
        )
    } else if(cor_method == "pearson") {
        MM <- WGCNA::signedKME(t(exp), MEs, outputColumnName = "")
    } else if(cor_method == "biweight") {
        MM <- WGCNA::signedKME(
            t(exp), MEs, corFnc = "bicor", outputColumnName = "",
            corOptions = "maxPOutliers = 0.1"
        )
    }

    # Add kME info
    kme_info <- merge(genes_modulesk, MM, by.x="Genes", by.y="row.names")

    # Keep top 10% genes with the highest degree for each module
    genes_mod_list <- split(kme_info, kme_info[,2])
    genes_mod_list$grey <- NULL
    top_10 <- lapply(genes_mod_list, function(x) {
        indices <- round(nrow(x) * 0.1)
        return(x[order(x$kWithin, decreasing=TRUE), ][seq_len(indices), ])
    })

    # Pick genes from the top 10% degree with kME above 0.8
    hubs <- Reduce(rbind, lapply(top_10, function(x) {
        cols <- c("Genes", "Modules", "kWithin", unique(x$Modules))
        h <- x[, cols]
        h <- h[h[,4] > 0.8, c(1,2,3)]
    }))
    rownames(hubs) <- seq_len(nrow(hubs))
    colnames(hubs) <- c("Gene", "Module", "kWithin")
    return(hubs)

    lapply(top_10, function(x) {
      name <- x$Modules %>% unique
      x[,c('Genes','Modules','kWithin',name)] %>% 
        rename(kME=all_of(name)) %>% 
        mutate(kIN_rank = rank(-kWithin),
               kME_rank = rank(-kME)
        ) %>% arrange(kME_rank) %>% head(1)
    }) %>% bind_rows %>% View()
