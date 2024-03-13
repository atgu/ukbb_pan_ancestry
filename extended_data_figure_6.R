source('figure3.R')

haploinsufficiency_plot = function(novel=T, sig_threshold=5e-8) {
  sig_genes = all_variants %>%
    filter(Pvalue < sig_threshold) %>%
    select(gene=nearest_genes) %>%
    separate_rows(gene) %>% distinct
  
  gene_lists %>% filter(gene_list == 'universe' & gene %in% all_genes$gene) %>%
    mutate(novel_gwas_hit = gene %in% all_novel_genes$gene,
           gwas_hit = gene %in% sig_genes$gene) %>% 
    select(gene, novel_gwas_hit, gwas_hit) -> annotated_genes
  
  genes_in_lists = annotated_genes %>%
    left_join(gene_lists, multiple = 'all') %>%
    mutate(gene_list = if_else(grepl('haploinsufficiency', gene_list), 'Haploinsufficient', gene_list),
           presence = 1)
  if (novel) {
    genes_in_lists %<>% mutate(hit=novel_gwas_hit)
  } else {
    genes_in_lists %<>% mutate(hit=gwas_hit)
  }
  
  gene_list_order = c('Haploinsufficient', 'Autosomal Dominant', 'Clinvar P/LP', 'Autosomal Recessive', 'All genes')
  
  gene_list_gwas_data = genes_in_lists %>%
    complete(gene_list, hit, fill = list(presence = 0)) %>%
    count(gene_list, hit, wt = presence) %>%
    group_by(gene_list) %>%
    mutate(prop_in_bin = n / sum(n)) %>%
    ungroup %>%
    mutate(
      gene_list = fct_recode(gene_list, 'Autosomal Recessive' = "all_ar", 'Autosomal Dominant' = "all_ad",
                             'Clinvar P/LP' = "clinvar_path_likelypath", 'All genes' = 'universe'),
      gene_list = fct_relevel(gene_list, gene_list_order))
  
  gene_list_gwas_data %>%
    count(gene_list, wt = n) %>%
    mutate(label=paste0(gene_list, '\n(n = ', comma(n), ')')) -> gene_list_labels_df
  
  gene_list_labels = gene_list_labels_df$label
  names(gene_list_labels) = gene_list_labels_df$gene_list
  
  gene_list_colors['Clinvar P/LP'] = colour_ramp(c(color_dominant, color_recessive))(0.5)
  gene_list_colors['All genes'] = 'lightgray'
  
  gene_list_gwas_data %>%
    filter(gene_list %in% gene_list_order) %>%
    ggplot + aes(x = gene_list, y = prop_in_bin, alpha = hit, fill=gene_list) +
    geom_bar(stat='identity') + 
    scale_fill_manual(values=gene_list_colors) +
    scale_y_continuous(labels=percent, name='Percentage of gene list') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    scale_x_discrete(labels=gene_list_labels, name=NULL) +
    scale_alpha_discrete(name=paste0(if_else(novel, 'Novel a', 'A'), 'ssociation\nnear gene')) +
    # theme(legend.position = 'bottom') +
    # guides(alpha = guide_legend(title.position="top", title.hjust = 0.5)) +
    guides(fill=F) %>% 
    return
}

ed6a = gene_list_gwas_data %>%
  filter(gene_list %in% gene_list_order) %>%
  ggplot + aes(x = gene_list, y = prop_in_bin, alpha = hit, fill=gene_list) +
  geom_bar(stat='identity') + 
  scale_fill_manual(values=gene_list_colors) +
  scale_y_continuous(labels=percent, name='Percentage of gene list') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_x_discrete(labels=gene_list_labels, name=NULL) +
  scale_alpha_discrete(name=paste0(if_else(novel, 'Novel a', 'A'), 'ssociation\nnear gene')) +
  # theme(legend.position = 'bottom') +
  # guides(alpha = guide_legend(title.position="top", title.hjust = 0.5)) +
  guides(fill=F) 

ed6a = haploinsufficiency_plot()
dmd = image_read_pdf("DMD_pan_ukb_X_31854782_META_500000_1Mb_continuous-23100-both_sexes--irnt.pdf")
ed6c = ggplot_pdf(dmd)


extended_data_figure6 = function(output_format = 'png', create_subplots=F) {
  # image_widths = c(image_info(pitx)$width, image_info(gdf5)$width)/300
  # image_heights = c(image_info(pitx)$height, image_info(gdf5)$height)/300

  # output_type(output_format, paste0('extended_data_figure6.', output_format), height=5.27, width=3.44*2)
  # # print(ggarrange(p1a, p1b, p1c, p1d, nrow = 2, ncol = 2))
  # right = ggarrange(ed6b + plot_annotation(tag_levels = list(c('b', ''))),
  #                   ed6c + plot_annotation(tag_levels = list(c('c', ''))), nrow=2
  # )
  # p = ggarrange(ed6a + plot_annotation(tag_levels = list(c('a', ''))), right, ncol=2)
  # print(p)
  # dev.off()
  scale = 1.2
  output_type(output_format, paste0('extended_data_figure6.', output_format), height=scale * 5.9 / 2, width=scale * 3.5 * 2)
  print(ggarrange(ed6a, ed6c, nrow = 1, ncol = 2))
  dev.off()
}
extended_data_figure6()
extended_data_figure6('pdf')

s_figure_32 = function(output_format = 'png', create_subplots=F) {
  p1 = haploinsufficiency_plot(F)
  output_type(output_format, paste0('s_figure_32.', output_format), height=4, width=4)
  print(p1)
  dev.off()
}
s_figure_32()