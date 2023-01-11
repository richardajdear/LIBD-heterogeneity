library(ggplot2)

get_overlap <- function(net_pair, relabel=TRUE) {
    # Get overlaps, optionally match, and reorder by size
    modules1 <- net_pair[[1]]$genes_and_modules$Modules
    modules2 <- net_pair[[2]]$genes_and_modules$Modules

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
    rename(y=Var1, x=Var2, n=value.x, p=value.y) %>%
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