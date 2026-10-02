IsoseqsSwitchList$isoformFeatures %>% 
    ggplot(aes(x=dIF, y=-log10(isoform_switch_q_value))) + 
    geom_point(
        aes( color=abs(dIF) > 0.5 & isoform_switch_q_value < 0.05 ), # default cutoff
        size=1
    ) +
    ggrepel::geom_text_repel(data = IsoseqsSwitchList$isoformFeatures %>% filter(abs(dIF) > 0.5 & isoform_switch_q_value < 0.05), aes(label=gene_name),size=3, max.overlaps = 20) +
    geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
    geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
    facet_wrap( ~ condition_2) +
    #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
    labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
    theme_bw()
ggsave('figures/Isoform_Switching_DexSeq_Isoseqs.png', width=8, height=6)

AstroDEGs <- read_tsv('./from_collaborator/AstroDEGs_stim_vs_unstim.txt')

AstroDEGs %>% 
    ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
    geom_point(
        aes( color=abs(log2FoldChange) > 1 & padj < 0.05 ), # custom cutoff
        size=1
    ) +
    geom_hline(yintercept = -log10(0.05), linetype='dashed') +
    geom_vline(xintercept = c(-1, 1), linetype='dashed') +
    scale_color_manual('Signficant\nDEG', values = c('black','red')) +
    labs(x='Log2 Fold Change', y='-Log10 ( Adjusted P Value )') +
    coord_cartesian(ylim=c(0, 10), xlim=c(-6,6)) +
    theme_bw()
ggsave('figures/Astrocyte_DEGs_stim_vs_unstim.png', width=8, height=6)