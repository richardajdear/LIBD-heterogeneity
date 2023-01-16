library(reshape2)
library(ggplot2)
library(pals)

plot_ngenes <- function(net) {
    counts <- net$colors %>% table %>% 
        as_tibble %>% rename_all(~ c('module','n')) %>% 
        arrange(-n)

    counts %>% 
    mutate(module = factor(module, ordered=TRUE, levels=.$module)) %>% 
    ggplot(aes(y=module, x=n, fill=module)) +
    geom_col() + scale_fill_identity() +
    geom_text(aes(label=n), hjust=0, angle=0) +
    scale_y_discrete('module', limits=rev) +
    coord_cartesian(clip='off') +
    xlab('n genes') +
    theme_classic() + theme(plot.margin=margin(r=20))
}

plot_topN_lines <- function(topN, xlim=1000) {
    topN %>% 
    rownames_to_column('n') %>% 
    pivot_longer(-n, names_to='module', values_to='matches') %>% 
    filter(!is.na(matches)) %>%
    mutate(n = as.numeric(n)) %>% 
    mutate(pct = matches/n) %>% 
    mutate(module = str_replace(module, "kME", "")) %>% 
    ggplot(aes(x=n, y=pct, color=module)) +
    geom_line() +
    scale_color_identity() +
    scale_y_continuous(limits=c(0, 1), name='% of top N genes matched',
                        breaks=seq(0,1,.2)) + 
    scale_x_continuous(limits=c(0, xlim), name='top N genes by kME', 
                        breaks=c(100,seq(0,xlim,250))) +
    theme_classic() +
    theme(panel.grid.major=element_line(size=.5))
}

plot_kME_distribution <- function(net) {
    net %>% filter_kME() %>% 
    data.frame %>% .[, 1:10] %>% 
    pivot_longer(everything(), names_to='module', values_to='kME') %>%
    filter(!is.na(kME), module!='kMEgrey') %>%
    mutate(module = str_replace(module, "kME", "")) %>%
    ggplot(aes(x=kME, color=module)) +
    geom_density() +
    scale_color_identity() +
    theme_classic() +
    ggtitle("kME distribution, largest 10 modules")
}

plot_kME_topN <- function(net) {
    net %>% filter_kME() %>%
    data.frame %>% .[, 1:10] %>%
    pivot_longer(everything(), names_to='module', values_to='kME') %>%
    filter(!is.na(kME), module!='kMEgrey') %>%
    mutate(module = str_replace(module, "kME", "")) %>%
    arrange(module, -kME) %>%
    group_by(module) %>%
    mutate(
        n = seq_along(kME),
        avg = cumsum(kME)/n) %>%
    ggplot(aes(x=n, y=kME, group=module, color=module)) +
    geom_line() +
    scale_color_identity() +
    scale_x_continuous(limits=c(0,1000), breaks=c(100, seq(0,1000,250))) +
    xlab('top Nth gene by kME') +
    ylab('kME mean') +
    theme_classic() +
    ggtitle("kME of top Nth gene, largest 10 modules")
}


get_overlap <- function(net_pair, relabel=TRUE) {
    # Get overlaps, optionally match, and reorder by size
    modules1 <- net_pair[[1]]$colors
    modules2 <- net_pair[[2]]$colors

    if (relabel) {
        modules2 <- WGCNA::matchLabels(modules2, modules1)
    }
    overlap <- WGCNA::overlapTable(modules1, modules2)

    all_colours <- WGCNA::standardColors(200)
    colours1 <- c(all_colours[all_colours %in% unique(modules1)], 'grey')
    colours2 <- c(all_colours[all_colours %in% unique(modules2)], 'grey')

    overlap$countTable <- overlap$countTable[colours1, colours2]
    overlap$pTable <- overlap$pTable[colours1, colours2]

    return(overlap)
}

plot_overlap_counts <- function(net_pair, relabel=TRUE) {
    overlap <- get_overlap(net_pair, relabel=relabel)
    colours_y <- rownames(overlap$countTable)
    colours_x <- colnames(overlap$countTable)

    overlap$countTable %>% melt() %>%
    left_join(overlap$pTable %>% melt, by=c("Var1", "Var2")) %>%
    dplyr::rename(y=Var1, x=Var2, n=value.x, p=value.y) %>%
    mutate(y=factor(y, ordered=T, levels=colours_y)) %>%
    mutate(x=factor(x, ordered=T, levels=colours_x)) %>%
    ggplot(aes(y=y, x=x, fill=n)) +
    geom_tile() +
    scale_fill_gradientn(colours=brewer.blues(100)) +
    scale_y_discrete(limits=rev) +
    #geom_text(aes(label=paste0(value.x, '\n', round(value.y,3)))) +
    geom_text(aes(label=n, color=ifelse(n>2000,'white','black'))) +
    scale_color_identity(guide='none') +
    theme_classic() +
    theme(
        axis.text.y=element_text(color=colours_y %>% rev),
        axis.text.x=element_text(angle=30,hjust=1,color=colours_x)
    )
}

get_overlap_pct <- function(net_pair, relabel=TRUE) {
    overlap <- get_overlap(net_pair, relabel=relabel)
    colours <- rownames(overlap$countTable)

    melt_table <- overlap$countTable %>% melt

    pct_table <- melt_table %>% 
    group_by(Var1) %>% summarise(total=sum(value)) %>% 
    mutate(Var2 = Var1) %>% 
    left_join(melt_table, by = c('Var1', 'Var2')) %>% 
    select(module=Var1, total, n=value) %>% 
    mutate(pct=n/total)

    return(pct_table)
}

plot_overlap_pct <- function(pct_table, min_points = 5, boxplot = TRUE) {
    colours <- unique(pct_table$module) %>% as.character()

    p <- pct_table %>%
    tidyr::replace_na(list(pct=0)) %>%
    group_by(module) %>% mutate(n_points=n()) %>%
    filter(n_points >= min_points) %>%
    ggplot(aes(x=module, y=pct)) +
    geom_point(aes(colour=module)) + scale_color_identity() +
    ylim(c(0,1)) +
    theme_classic() +
    theme(axis.text.x=element_text(angle=30,hjust=1,color=colours))

    if (boxplot) {
        p <- p + geom_boxplot(aes(fill=module), alpha=.3) + scale_fill_identity()
    } else {
        # p <- p + geom_point(aes(colour=module)) + scale_color_identity()
    }
    return(p)
}