source('~/ukbb_pan_ancestry/constants.R')

efo = read_tsv(gzfile('data/known_ukbb_loci.meta_hq_annotated.txt.bgz'),
               col_types = cols(locus=col_character(),
                                lead_locus=col_character(),
                                locus_otg=col_character()))

known_novel_efo = function() {
  # 85,960
  efo %>% select(idx) %>% distinct %>% nrow
  
  # 35,950 (terms)
  efo %>% 
    filter(trait_efo_term == trait_efo_term_otg) %>%
    distinct(idx) %>%
    nrow
  
  # 14,708 (no term)
  efo %>% 
    filter(is.na(trait_efo_term)) %>%
    distinct(idx) %>%
    nrow
  
  # novel = 35302
  85960-35950-14708
  
  # 65,821
  efo %>% 
    filter(trait_efo_category == trait_efo_category_otg) %>%
    distinct(trait_efo_category, idx) %>%
    nrow
  
  # novel = 35302
  85960-65821-14708
  
  efo %>%
    group_by(trait_efo_category, idx) %>%
    summarize(any_same_category=any(trait_efo_category == trait_efo_category_otg)) %>%
    group_by(trait_efo_category) %>%
    summarize(total_count=n(),
              known_count=sum(any_same_category, na.rm=T),
              novel_count=total_count-known_count,
              known_frac=known_count/total_count,
              novel_frac=1-known_frac) %>%
    mutate(trait_efo_category = fct_reorder(trait_efo_category, novel_frac)) %>%
    pivot_longer(-trait_efo_category) %>%
    separate(name, c('novelty', 'scale'), ('_')) %>%
    filter(!is.na(trait_efo_category)) -> data
  
  efo %>%
    filter(trait_efo_category == 'Metabolic disorder' & trait_efo_category_otg != 'Metabolic disorder') %>%
    select(lead_locus, lead_alleles, nearest_genes, description) %>% distinct
  
  efo %>%
    filter(trait_efo_category == 'Cancer' & trait_efo_category_otg != 'Cancer') %>%
    count(description, description_more)
  
  # data %>%
  #   mutate(value=if_else(scale == 'frac', value, -value)) %>%
  #   ggplot() +
  #   aes(y = value, x = trait_efo_category, fill = novelty) +
  #   geom_bar(stat='identity') +
  #   facet_share(~scale, dir = "h", scales = "free", reverse_num = T,
  #               labeller = labeller(scale = c('count'='Number of hits', 'frac'='Fraction novel'))) + 
  #   coord_flip() +
  #   theme_classic() +
  #   ylab(NULL) + xlab(NULL) +
  #   scale_fill_discrete(labels=c('known' = 'Known', 'novel' = 'Novel'), name=NULL)
  
  fill_colors = c(
    'known' = muted('lightblue'),
    'novel' = muted('orange', l=70, c=100)
  )
  data %>%
    filter(scale == 'count' & novelty == 'novel') %>% 
    arrange(desc(value))
  p_num = data %>%
    filter(scale == 'count' & novelty != 'total') %>%
    ggplot +
    aes(x = value, y = trait_efo_category, fill = novelty, label=comma(value)) +
    geom_bar(stat='identity', position='stack') +
    geom_text(data=data %>% filter(novelty == "total"), hjust=1, nudge_x = -200, size=2.5) +
    theme_classic() + scale_x_reverse(labels=comma, name='Number of significant loci',
                                      expand = expansion(c(0.2, 0))) +
    scale_y_discrete(position='right') +
    theme(axis.text.y=element_blank(),
          legend.position=c(0.01, 0.99),
          legend.justification=c(0, 1)) +
    ylab(NULL) +
    scale_fill_manual(labels=c('known' = 'Known', 'novel' = 'Novel'), name=NULL,
                      values=fill_colors)
  
  data2 = data %>%
    filter(novelty != 'total' & scale == 'count') %>%
    arrange(desc(novelty)) %>%
    group_by(trait_efo_category) %>%
    mutate(frac = value / sum(value),
           x =  cumsum(frac) - frac / 2) %>% ungroup
  
  p_frac = data2 %>%
    ggplot +
    aes(x = frac, y = trait_efo_category, fill = novelty) +
    geom_col() +
    geom_text(aes(x=x, label=percent(frac, accuracy=1), color=novelty),
              data=data2 %>% filter(frac > 0.1),
              hjust=0.5, size=2.5) +
    theme_classic() + scale_x_continuous(labels=percent, name='Fraction novel') +
    theme(axis.text.y=element_blank()) +
    scale_color_manual(values=c('known' = 'white', 'novel' = 'black'), guide=F) +
    ylab(NULL) + guides(fill=F) +
    scale_fill_manual(labels=c('known' = 'Known', 'novel' = 'Novel'), name=NULL,
                      values=fill_colors)
  # p_frac
  p_labels = data %>%
    filter(novelty == 'known' & scale == 'count') %>%
    ggplot + aes(y = trait_efo_category, label = trait_efo_category) +
    geom_text(x = 0.5, size = 3, hjust = 0.5) + 
    theme(axis.text.y=element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank()) + xlab(NULL) + ylab(NULL)
  
  p = egg::ggarrange(p_num + theme(plot.margin = margin(0,0,0,8, "pt")),
                     p_labels + theme(plot.margin = margin(0,0,0,0, "pt")),
                     p_frac + theme(plot.margin = margin(0,8,0,0, "pt")),
                     nrow=1, ncol=3)
  # return(p)
  return(list(p_num + theme(plot.margin = margin(0,0,0,8, "pt")),
              p_labels + theme(plot.margin = margin(0,0,0,0, "pt")),
              p_frac + theme(plot.margin = margin(0,8,0,0, "pt"))))
}



p4a_subs = known_novel_efo()
p4b = haploinsufficiency_plot()
# pitx = readPNG("PITX2_pan_ukb_4_112282681_META_800000_continuous-5099-both_sexes--irnt.png")
pitx = image_read_pdf("PITX2_pan_ukb_4_112282681_META_650000_1Mb_continuous-5099-both_sexes--irnt_test.pdf")
p4c = ggplot() +
  annotation_raster(pitx, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_blank() + theme_nothing()
gdf5 = image_read_pdf("GDF5_pan_ukb_20_34025756_META_500000_1Mb_phecode-716.2-both_sexes--.pdf")
p4d = ggplot() +
  annotation_raster(gdf5, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_blank() + theme_nothing()


figure4 = function(output_format = 'png', create_subplots=F) {
  # image_widths = c(image_info(pitx)$width, image_info(gdf5)$width)/300
  # image_heights = c(image_info(pitx)$height, image_info(gdf5)$height)/300
  
  
  relative_heights = c(0.8, 1)
  factor = 0.9
  width = 3.44 * 2 * factor
  height = 5.27 * sum(relative_heights) * factor
  if (create_subplots) {
    fig4_plots = list(p4a, p4b, p4c)
    for (i in seq_along(fig4_plots)) {
      output_type(output_format, paste0('figure4_panel', i, '.', output_format),
                  height=4, width=6)
      print(fig4_plots[i])
      dev.off()
    }
  } else {
    output_type(output_format, paste0('figure4.', output_format), height=height, width=width)
    # print(ggarrange(p1a, p1b, p1c, p1d, nrow = 2, ncol = 2))
    top = (p4a_subs[[1]] | p4a_subs[[2]] | p4a_subs[[3]]) + plot_annotation(tag_levels = list(c('a', '')))
    right = ggarrange(p4c + plot_annotation(tag_levels = list(c('c', ''))),
                      p4d + plot_annotation(tag_levels = list(c('d', ''))), nrow=2
    )
    bottom = ggarrange(p4b + plot_annotation(tag_levels = list(c('b', ''))), right, ncol=2)
    print(ggarrange(top, bottom, nrow=2, ncol=1, align='none', heights = relative_heights))
    dev.off()
  }
}
figure4()
figure4('pdf')

all_variants = read_tsv(gzfile('data/top_meta_with_top_pop_by_variant_full.txt.bgz'))

all_genes = all_variants %>%
  select(gene=nearest_genes) %>%
  separate_rows(gene) %>% distinct

gene_lists = load_all_gene_list_data()

all_variants %>%
  mutate(Pvalue=if_else(Pvalue < 1e-10, 1e-10, Pvalue)) %>%
  ggplot + aes(x = Pvalue) + geom_histogram(bins=100) + scale_x_log10()

sig_genes = all_variants %>%
  filter(Pvalue < 5e-8) %>%
  select(gene=nearest_genes) %>%
  separate_rows(gene) %>% distinct

novel_hits = efo %>%
  group_by(trait_efo_category, idx) %>% # , trait_type, phenocode, pheno_sex, coding, modifier, lead_locus, nearest_genes) %>%
  mutate(any_same_category=any(trait_efo_category == trait_efo_category_otg, na.rm=T)) %>%
  filter(!any_same_category) %>%
  select(matches("lead_locus") | (!contains("otg") & !matches("locus"))) %>%
  distinct %>% ungroup

novel_hits %>%
  select(gene=nearest_genes) %>%
  separate_rows(gene) %>%
  distinct -> all_novel_genes

gene_lists %>% filter(gene_list == 'universe' & gene %in% all_genes$gene) %>%
  mutate(novel_gwas_hit = gene %in% all_novel_genes$gene,
         gwas_hit = gene %in% sig_genes$gene) %>% 
  select(gene, novel_gwas_hit, gwas_hit) -> annotated_genes

genes_in_lists = annotated_genes %>%
  left_join(gene_lists, multiple = 'all') %>%
  mutate(gene_list = if_else(grepl('haploinsufficiency', gene_list), 'Haploinsufficient', gene_list),
         presence = 1)

genes_in_lists %>% 
  filter(gene == 'ATP8A2')

gene_list_order = c('Haploinsufficient', 'Autosomal Dominant', 'Clinvar P/LP', 'Autosomal Recessive', 'All genes')

gene_list_gwas_data = genes_in_lists %>%
  complete(gene_list, novel_gwas_hit, fill = list(presence = 0)) %>%
  count(gene_list, novel_gwas_hit, wt = presence) %>%
  group_by(gene_list) %>%
  mutate(prop_in_bin = n / sum(n)) %>%
  ungroup %>%
  mutate(
    gene_list = fct_recode(gene_list, 'Autosomal Recessive' = "all_ar", 'Autosomal Dominant' = "all_ad",
                           'Clinvar P/LP' = "clinvar_path_likelypath", 'All genes' = 'universe'),
    gene_list = fct_relevel(gene_list, gene_list_order))

gene_list_gwas_data %>%
  count(gene_list, wt = n) %>%
  mutate(label=paste0(gene_list, '\n(n = ', n, ')')) -> gene_list_labels_df

gene_list_labels = gene_list_labels_df$label
names(gene_list_labels) = gene_list_labels_df$gene_list

gene_list_colors['Clinvar P/LP'] = colour_ramp(c(color_dominant, color_recessive))(0.5)
gene_list_colors['All genes'] = 'lightgray'

p4b = gene_list_gwas_data %>%
  filter(gene_list %in% gene_list_order) %>%
  ggplot + aes(x = gene_list, y = prop_in_bin, alpha = novel_gwas_hit, fill=gene_list) +
  geom_bar(stat='identity') + 
  scale_fill_manual(values=gene_list_colors) +
  scale_y_continuous(labels=percent, name='Percentage of gene list') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_x_discrete(labels=gene_list_labels, name=NULL) +
  guides(fill=F) + scale_alpha_discrete(name='Novel assocation near gene') +
  theme(legend.position = 'bottom') +
  guides(alpha = guide_legend(title.position="top", title.hjust = 0.5))

novel_hits %>%
  rename(gene=nearest_genes) %>%
  separate_rows(gene) %>%
  left_join(gene_lists, multiple = 'all') %>%
  filter(grepl('haploinsufficiency', gene_list)) %>%
  View


novel_hits %>% count(trait_efo_category)

novel_hits %>% 
  mutate(freq_max = pmax(AFR_af, AMR_af, CSA_af, EAS_af, MID_af),
         fc = freq_max / EUR_af) %>% 
  filter(fc > 2) %>%
  filter(!is.na(trait_efo_category) & trait_efo_category != 'Other measurement') %>%
  filter(trait_type %in% c('icd10', 'phecode')) %>%
  arrange(desc(freq_max)) %>% View

novel_hits %>% 
  filter(EUR_af < 0.05 & AFR_af > 0.05) %>%
  filter(!is.na(trait_efo_category) & trait_efo_category != 'Other measurement') %>%
  count(trait_type, phenocode, pheno_sex, coding, modifier) %>% count
