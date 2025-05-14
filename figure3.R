source('~/ukbb_pan_ancestry/constants.R')

reticulate::use_miniconda('r-reticulate')

pheno_summary = load_ukb_file('pheno_summary_full.txt.bgz', 'phenos/')

efo = load_ukb_file('known_ukbb_loci.meta_hq_annotated.txt.bgz', subfolder = 'known_novel/',
                    force_cols = cols(locus=col_character(),
                                lead_locus=col_character(),
                                locus_otg=col_character()))

efo_no_exclusion = load_ukb_file('known_ukbb_loci_no_exclusion.meta_hq_annotated.txt.bgz',subfolder = 'known_novel/',
                                 force_cols = cols(locus=col_character(),
                                         lead_locus=col_character(),
                                         locus_otg=col_character()))

meta_eur_comparison_all = load_ukb_file('meta_eur_comparison.tsv.bgz', subfolder = 'known_novel/', force_cols = cols(coding=col_character()))
meta_eur_comparison_all %>% filter(N_pops > 1) %>% mutate(
  # nlog10p(meta) - nlog10p(EUR) = log10(P[EUR]/P[Meta])
  ratio = Pvalue_meta - Pvalue_EUR,
  freq_max = pmax(freq_AFR, freq_AMR, freq_CSA, freq_EAS, freq_MID),
  sig_in=case_when(Pvalue_meta > log_threshold & Pvalue_EUR < log_threshold ~ 'meta_only',
                   Pvalue_meta < log_threshold & Pvalue_EUR > log_threshold ~ 'EUR_only',
                   Pvalue_meta > log_threshold & Pvalue_EUR > log_threshold ~ 'both',
                   TRUE ~ 'neither'),
  orig_sig_in=case_when(Pvalue_meta > log_orig_threshold & Pvalue_EUR < log_orig_threshold ~ 'meta_only',
                        Pvalue_meta < log_orig_threshold & Pvalue_EUR > log_orig_threshold ~ 'EUR_only',
                        TRUE ~ 'both')
) %>% filter(Pvalue_meta > log_orig_threshold | Pvalue_EUR > log_orig_threshold) -> meta_eur_comparison

meta_eur_comparison %>%
  filter(sig_in == 'meta_only') %>%
  count(phenocode, nearest_genes, description) %>%
  count(nearest_genes) %>%
  arrange(desc(n)) %>%
  print(n=20)

# 11.5 million significant associations, 237,360 LD-independent
meta_eur_comparison %>%
  filter(Pvalue_meta > log_orig_threshold & high_quality) %>%
  summarize(clump=sum(!is.na(in_clump)), n=n())

# spanning 431 phenotypes
meta_eur_comparison %>%
  filter(Pvalue_meta > log_orig_threshold & high_quality) %>%
  select_at(key_fields) %>% distinct %>% nrow

# 511,841 associations were not significant in EUR (14,676 LD-independent loci)
meta_eur_comparison %>%
  filter(Pvalue_meta > log_orig_threshold & high_quality & Pvalue_EUR < log_orig_threshold) %>%
  summarize(clumped=sum(!is.na(in_clump)), n=n())

# 19,493 have at least one genome-wide significant association
meta_eur_comparison %>%
  filter(Pvalue_meta > log_orig_threshold & high_quality & freq_max >= 0.01 & freq_EUR < 0.01) %>%
  select(chrom, pos, ref, alt) %>% distinct %>% nrow

# Not included stat: 43,230 total associations among 3,485 LD-independent associations
meta_eur_comparison %>%
  filter(Pvalue_meta > log_orig_threshold & high_quality & freq_max >= 0.01 & freq_EUR < 0.01) %>%
  summarize(clumped=sum(!is.na(in_clump)), n=n())

# 5,259 associations that are not significant in EUR (559 LD-independent associations)
meta_eur_comparison %>%
  filter(Pvalue_meta > log_orig_threshold & high_quality & freq_max >= 0.01 & freq_EUR < 0.01 & Pvalue_EUR < log_orig_threshold) %>%
  summarize(clumped=sum(!is.na(in_clump)), n=n())

plot_n_sig_per_pheno = function(indep_only=T) {
  plot_data = meta_eur_comparison %>%
    filter(Pvalue_meta > -log10(5e-8) & high_quality)
  
  if (indep_only) plot_data %<>% filter(!is.na(in_clump))
  axis_name = paste0('Number of significant', if_else(indep_only, ' independent', ''), ' associations')
  
  plot_data %>%
    group_by_at(key_fields) %>%
    summarize(n=n()) %>%
    ggplot + aes(x = n) +
    geom_histogram(bins=50) + 
    scale_x_log10(labels=comma, name=axis_name) +
    ylab('Number of phenotypes') %>%
    return
}

meta_eur_comparison %>%
  lm(Pvalue_meta ~ Pvalue_EUR, data=.) %>% 
  summary -> res
print(res$coefficients)

meta_eur_freq_fc = function(pop = 'CSA', threshold=1e-10, clump_only=F, suggestive_line=T) {
  plot_data = meta_eur_comparison %>%
    filter(Pvalue_meta > -log10(threshold) | Pvalue_EUR > -log10(threshold)) 
  if(pop == 'max') {
    high_col = '#008080'
    low_col = 'black'
    label_name = expression(paste(log[10], '(', frac(max~AF, EUR~AF), ')'))
  } else {
    high_col = ukb_pop_colors[[pop]]
    low_col = color_eur
    label_name = paste0('log10(', pop, ' AF / EUR AF)')
  }
  if (clump_only) {
    plot_data %<>% filter(!is.na(in_clump))
  }
  
  ggplot(data=plot_data %>% head) +
    scattermore::geom_scattermore(
      aes(Pvalue_EUR, Pvalue_meta, color = log10(.data[[paste0('freq_', pop)]]/.data[['freq_EUR']])),
      pointsize = 2, pixels = c(1024, 1024), alpha = 0.5,
      data = plot_data # %>% head(1000)
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    scale_x_continuous(trans = locusviz::trans_loglog_p(), labels=comma) +
    scale_y_continuous(trans = locusviz::trans_loglog_p(), labels=comma) +
    scale_color_gradient2(high = high_col, low = low_col, mid = 'gray90',
                          name=label_name) +
    geom_segment(x = 0, xend = -log10(threshold), y = -log10(threshold), yend = -log10(threshold), color = "black", linetype='dashed') +
    geom_segment(x = -log10(threshold), xend = -log10(threshold), y = 0, yend = -log10(threshold), color = "black", linetype='dashed') +
    labs(x = expression(paste(-log[10], "P (EUR)")), y = expression(paste(-log[10], "P (Meta-analysis)"))) +
    guides(color=guide_colorbar(title.position='top', title.hjust=1, direction = 'horizontal')) +
    theme(legend.position = c(1, 0), legend.justification = c(1, 0),
          legend.background=element_rect(fill = alpha("white", 0)),
          plot.margin = margin(8, 16, 0, 0, "pt")) -> plt
  if (suggestive_line) {
    suggestive_threshold = 1e-6
    color = 'grey50'
    linetype = 'dotted'
    ymax = 15
    plt = plt + geom_segment(x = -log10(suggestive_threshold), xend=-log10(suggestive_threshold),
                             y=-log10(threshold), yend=ymax, linetype=linetype, color=color) +
      geom_segment(x = -log10(threshold), xend=-log10(threshold),
                   y=-log10(threshold), yend=ymax, linetype=linetype, color=color)
  }
  # plt
  return(plt)
}
p3a = meta_eur_freq_fc('max', threshold = 5e-8)
p3a

get_efo_stats = function(dat, category_only=F) {
  if (!category_only) {
    print('Total:')
    dat %>% select(idx) %>% distinct %>% nrow %>% print
    
    print('Missing:')
    dat %>% 
      group_by(idx) %>%
      summarize(all_missing=all(is.na(trait_efo_term))) %>% 
      count(all_missing) %>% print
    
    print('Novelty:')
    dat %>%
      filter(!is.na(trait_efo_term)) %>%
      mutate(trait_efo_term_otg = if_else(is.na(trait_efo_term_otg), 'none', trait_efo_term_otg)) %>%
      group_by(idx) %>%
      summarize(known=any(trait_efo_term == trait_efo_term_otg)) %>%
      count(known) %>% print
  }
  print('Category:')
  dat %>%
    filter(!is.na(trait_efo_category)) %>%
    mutate(trait_efo_category_otg = if_else(is.na(trait_efo_category_otg), 'none', trait_efo_category_otg)) %>%
    group_by(idx) %>%
    summarize(known=any(trait_efo_category == trait_efo_category_otg)) %>%
    count(known) %>% print
  
  # print('Category without X:')
  # dat %>%
  #   filter(!is.na(trait_efo_category)) %>%
  #   filter(!str_starts(lead_locus, 'X')) %>%
  #   mutate(trait_efo_category_otg = if_else(is.na(trait_efo_category_otg), 'none', trait_efo_category_otg)) %>%
  #   group_by(idx) %>%
  #   summarize(known=any(trait_efo_category == trait_efo_category_otg)) %>%
  #   count(known) %>% print
  
  if (!category_only) {
    print('X total:')
    dat %>% filter(str_starts(lead_locus, 'X')) %>%
      select(idx) %>% distinct %>% nrow %>% print
    
    print('X missing:')
    dat %>% 
      filter(str_starts(lead_locus, 'X')) %>%
      group_by(idx) %>%
      summarize(all_missing=all(is.na(trait_efo_term))) %>% 
      count(all_missing) %>% print
    
    print('X novelty:')
    dat %>%
      filter(!is.na(trait_efo_term)) %>%
      filter(str_starts(lead_locus, 'X')) %>%
      mutate(trait_efo_term_otg = if_else(is.na(trait_efo_term_otg), 'none', trait_efo_term_otg)) %>%
      group_by(idx) %>%
      summarize(known=any(trait_efo_term == trait_efo_term_otg)) %>%
      count(known) %>% print
  
    print('X categories matching:')
    dat %>%
      filter(!is.na(trait_efo_category)) %>%
      filter(str_starts(lead_locus, 'X')) %>%
      mutate(trait_efo_category_otg = if_else(is.na(trait_efo_category_otg), 'none', trait_efo_category_otg)) %>%
      group_by(idx) %>%
      summarize(known=any(trait_efo_category == trait_efo_category_otg)) %>%
      count(known) %>% print
  }
}

explorations = function() {
  
  efo_no_exclusion %>%
    filter(!is.na(trait_efo_category)) %>%
    mutate(trait_efo_category_otg = if_else(is.na(trait_efo_category_otg), 'none', trait_efo_category_otg)) %>%
    group_by(idx) %>%
    filter(!any(trait_efo_category == trait_efo_category_otg)) %>% View
  
  efo_no_exclusion %>%
    filter(!is.na(trait_efo_category)) %>%
    filter(!str_starts(lead_locus, 'X') & trait_type != 'biomarkers') %>%
    mutate(trait_efo_category_otg = if_else(is.na(trait_efo_category_otg), 'none', trait_efo_category_otg)) %>%
    group_by(idx) %>%
    filter(!any(trait_efo_category == trait_efo_category_otg)) %>% ungroup %>%
    distinct(trait_type, phenocode, coding, modifier, lead_locus, lead_alleles, rsid, description, description_more, idx, trait_efo, trait_efo_term, trait_efo_category) %>%
    # mutate(X=str_starts(lead_locus, 'X')) %>%
    count(trait_type, phenocode, modifier, coding, description) %>% arrange(desc(n)) %>% View
  
  efo_no_exclusion %>%
    filter(!is.na(trait_efo_category)) %>%
    filter(!str_starts(lead_locus, 'X') & trait_type != 'biomarkers') %>%
    mutate(trait_efo_category_otg = if_else(is.na(trait_efo_category_otg), 'none', trait_efo_category_otg)) %>%
    group_by(idx) %>%
    filter(!any(trait_efo_category == trait_efo_category_otg)) %>% ungroup %>%
    distinct(trait_type, phenocode, coding, modifier, lead_pvalue, lead_locus, lead_alleles, rsid, description, description_more, idx, trait_efo, trait_efo_term, trait_efo_category) %>%
    # filter(phenocode == "21001") %>%
    filter(trait_type %in% c('icd10', 'phecode')) %>%
    View
  
  efo_no_exclusion %>%
    filter(!is.na(trait_efo_category)) %>%
    mutate(trait_efo_category_otg = if_else(is.na(trait_efo_category_otg), 'none', trait_efo_category_otg)) %>%
    group_by(idx) %>%
    filter(!any(trait_efo_category == trait_efo_category_otg)) %>% ungroup %>%
    filter(phenocode == "23104" & rsid == "rs2582842" 
           # & grepl('NEALE', study_id_otg)
    ) %>% View
  
  efo_no_exclusion %>%
    filter(phenocode == "30100" & rsid == "rs73633938")
  
  
  efo %>%
    filter(!is.na(trait_efo_term)) %>%
    mutate(trait_efo_term_otg = if_else(is.na(trait_efo_term_otg), 'none', trait_efo_term_otg)) %>%
    group_by(idx, trait_type, phenocode, pheno_sex, coding, modifier, lead_locus, lead_alleles) %>%
    summarize(known=any(trait_efo_term == trait_efo_term_otg), .groups = 'drop') %>%
    filter(!known) %>%
    left_join(meta_eur_comparison %>%
                mutate(lead_locus=paste(chrom, pos, sep=':'),
                       lead_alleles=paste0('[\"', ref, '\",\"', alt, '\"]'),
                       coding=as.character(coding)) %>%
                select(trait_type, phenocode, pheno_sex, coding, modifier, lead_locus, lead_alleles,
                       Pvalue_meta, Pvalue_EUR, high_quality, description, nearest_genes)) %>%
    filter(Pvalue_meta > log_orig_threshold & high_quality & Pvalue_EUR < log_orig_threshold) %>% View
  
  efo %>%
    filter(!is.na(trait_efo_category)) %>%
    mutate(trait_efo_category_otg = if_else(is.na(trait_efo_category_otg), 'none', trait_efo_category_otg)) %>%
    group_by(idx, trait_type, phenocode, pheno_sex, coding, modifier, lead_locus, lead_alleles) %>%
    summarize(known=any(trait_efo_category == trait_efo_category_otg), .groups = 'drop') %>%
    filter(!known) %>%
    left_join(meta_eur_comparison %>%
                mutate(lead_locus=paste(chrom, pos, sep=':'),
                       lead_alleles=paste0('[\"', ref, '\",\"', alt, '\"]'),
                       coding=as.character(coding)) %>%
                select(trait_type, phenocode, pheno_sex, coding, modifier, lead_locus, lead_alleles,
                       Pvalue_meta, Pvalue_EUR, freq_max, freq_EUR, high_quality, description, nearest_genes)) %>%
    filter(Pvalue_meta > log_orig_threshold & high_quality & Pvalue_EUR < log_orig_threshold) %>% View
  # 429 remain
  # filter(high_quality & freq_max > freq_EUR * 2 & Pvalue_meta > Pvalue_EUR) %>% View
  
  efo %>%
    filter(trait_efo_category == 'Metabolic disorder' & trait_efo_category_otg != 'Metabolic disorder') %>%
    select(lead_locus, lead_alleles, nearest_genes, description) %>% distinct
  
  efo %>%
    filter(trait_efo_category == 'Cancer' & trait_efo_category_otg != 'Cancer') %>%
    count(description, description_more)
}

all_variants = load_ukb_file('top_meta_with_top_pop_by_variant_full.txt.bgz', subfolder='sig_hits/')

all_genes = all_variants %>%
  select(gene=nearest_genes) %>%
  separate_rows(gene) %>% distinct

gene_lists = load_all_gene_list_data()

efo %>%
  group_by(trait_efo_category, idx) %>% # , trait_type, phenocode, pheno_sex, coding, modifier, lead_locus, nearest_genes) %>%
  mutate(any_same_category=any(trait_efo_category == trait_efo_category_otg, na.rm=T)) %>%
  filter(!any_same_category) %>%
  select(matches("lead_locus") | (!contains("otg") & !matches("locus"))) %>%
  distinct %>% ungroup %>%
  select(gene=nearest_genes) %>%
  separate_rows(gene) %>%
  distinct -> all_novel_genes

# Of 22,776,573 variants 
nrow(all_variants)

# 1,589,664 had at least one significant association
all_variants %>%
  filter(Pvalue < 5e-8) %>%
  nrow

# near 17,285 genes
all_variants %>%
  filter(Pvalue < 5e-8) %>%
  select(gene=nearest_genes) %>%
  separate_rows(gene) %>% distinct %>%
  nrow

# of 19,842 genes
all_genes %>% nrow

# We filtered the 19,842 available genes to 17,428 that are in the gene list "universe"
gene_lists %>% filter(gene_list == 'universe' & gene %in% all_genes$gene) %>% nrow

# CAMK2D:
meta_eur_comparison %>% filter(pos == 114604622 & phenocode == 30870) %>% t
pheno_summary %>% filter(phenocode == 30870) %>%
  filter(pheno_data.pop %in% c('AFR', 'CSA', 'EUR', 'MID')) %>%
  summarize(n=sum(pheno_data.n_cases))

# PITX2:
meta_eur_comparison %>% filter(pos == 112282681 & phenocode == 5099) %>% t
pheno_summary %>% filter(phenocode == 5099) %>%
  filter(pheno_data.pop %in% c('CSA', 'EUR')) %>%
  summarize(n=sum(pheno_data.n_cases))

# DMD:
meta_eur_comparison %>% filter(pos == 31854782 & phenocode == 23100) %>% t
pheno_summary %>% filter(phenocode == 23100) %>%
  filter(pheno_data.pop %in% c("AFR", "CSA", "EAS", "EUR", "MID")) %>%
  summarize(n=sum(pheno_data.n_cases))

# G6PD:
meta_eur_comparison %>% filter(pos == 153764217 & phenocode == 30750) %>% t
pheno_summary %>% filter(phenocode == 30750) %>%
  filter(pheno_data.pop %in% c("AFR", "AMR", "EUR", "MID")) %>%
  summarize(n=sum(pheno_data.n_cases))

funnel_plot = function() {
  total_hits = efo %>% 
    group_by(idx) %>%
    summarize(all_missing=all(is.na(trait_efo_term))) %>% 
    filter(!all_missing) %>% nrow
  
  novel_hits_term = efo %>%
    filter(!is.na(trait_efo_term)) %>%
    mutate(trait_efo_term_otg = if_else(is.na(trait_efo_term_otg), 'none', trait_efo_term_otg)) %>%
    group_by(idx) %>%
    summarize(known=any(trait_efo_term == trait_efo_term_otg)) %>%
    filter(!known) %>% nrow
  
  novel_hits_category = efo %>%
    filter(!is.na(trait_efo_category)) %>%
    mutate(trait_efo_category_otg = if_else(is.na(trait_efo_category_otg), 'none', trait_efo_category_otg)) %>%
    group_by(idx) %>%
    summarize(known=any(trait_efo_category == trait_efo_category_otg)) %>%
    filter(!known) %>% nrow
  
  # was 4915
  minimal = efo_no_exclusion %>% 
    select(idx, trait_efo_category_otg, trait_efo_category,
           AFR_af, AMR_af, CSA_af, EAS_af, MID_af, EUR_af,
           chrom, pos, ref, alt, trait_type, phenocode, pheno_sex, coding, modifier) %>%
    filter(!is.na(trait_efo_category))
  
  minimal = minimal %>%
    mutate(freq_max = pmax(AFR_af, AMR_af, CSA_af, EAS_af, MID_af),
           coding=as.character(coding),
           match=trait_efo_category_otg == trait_efo_category) %>%
    select(idx, match, freq_max,
           chrom, pos, ref, alt, trait_type, phenocode, pheno_sex, coding, modifier,
           EUR_af)
  
  novel_hits_category_no_exclusion = minimal %>%
    select(idx, match) %>%
    group_by(idx) %>%
    summarize(known=any(match, na.rm=T)) %>%
    filter(!known) %>% nrow
  
  minimal %>%
    # mutate(trait_efo_category_otg = if_else(is.na(trait_efo_category_otg), 'none', trait_efo_category_otg)) %>%
    left_join(meta_eur_comparison_all %>%
                select(chrom, pos, ref, alt, trait_type, phenocode, pheno_sex, coding, modifier, Pvalue_meta, Pvalue_EUR),
              by=c("chrom", "pos", "ref", "alt", key_fields)) -> minimal2
  
  minimal2 %>%
    mutate(enriched = freq_max / EUR_af >= 10,
           eur_suggestive=Pvalue_EUR >= -log10(1e-6)) -> test
  
  # test %>%
  #   filter(!is.na(trait_efo_category)) %>%
  #   mutate(trait_efo_category_otg = if_else(is.na(trait_efo_category_otg), 'none', trait_efo_category_otg)) %>%
  #   group_by(idx) %>%
  #   filter(!any(trait_efo_category == trait_efo_category_otg) & !eur_suggestive) %>% ungroup %>%
  #   distinct(trait_type, phenocode, coding, modifier, lead_pvalue, chrom, pos, ref, alt, rsid, description, description_more, idx, trait_efo, trait_efo_term, trait_efo_category) %>%
  #   View
  
  # X chromsome novelty (573 / 2448)
  minimal %>%
    filter(chrom == 'X') %>%
    select(idx, match) %>%
    group_by(idx) %>%
    summarize(known=any(match, na.rm=T)) %>%
    count(!known) %>% print
  
  res = test %>%
    group_by(idx) %>%
    summarize(x=any(chrom=='X'),
              known=any(match, na.rm=T),
              enriched=any(enriched),
              eur_suggestive=any(eur_suggestive, na.rm=T)) 
  res %>% filter(!known) %>% count(x, enriched)
  res %>% filter(!known) %>% count(enriched, eur_suggestive)
  enriched_n = res %>% filter(!known & enriched) %>% nrow 
  suggestive_n = res %>% filter(!known & eur_suggestive) %>% nrow 
  
  test %>%
    filter(chrom == 'X') %>%
    group_by(idx) %>%
    summarize(known=any(match, na.rm=T),
              enriched=any(enriched),
              eur_suggestive=any(eur_suggestive, na.rm=T)) %>%
    count(enriched, eur_suggestive)
  
  # overall novelty: 5249, 3023 w/o X, 2544 w/o X or biomarkers
  # no exclusion:    4915, 2689 w/o X, 2213 w/o X or biomarkers
  # no UMich-SAIGE:  4921, 2695 w/o X
  
  start = 60
  end = 80
  colors = colorRampPalette(c('#2E8BC0', '#B1D4E0'))(4)
  colors = c(colors, colors[4])
  m <- list(
    l = 0,
    r = 0,
    b = 25,
    t = 25,
    pad = 0
  )
  
  p <- plot_ly(
    type = "funnelarea",
    values = c(4.5, 3.5, 3, 7.5, 0),
    textinfo='text',
    aspectratio=1/2,
    baseratio=1/2,
    text = c(paste0("Independent genome-wide significant associations:\n", comma(total_hits)),
             paste0("Novel associations (EFO term), excluding previous UKB:\n", comma(novel_hits_term)),
             paste0("Novel associations (EFO category), excluding previous UKB:\n", comma(novel_hits_category)),
             paste0("Novel associations (EFO category)\nincluding previous UKB:\n", comma(novel_hits_category_no_exclusion),
                    "\n• Suggestive in EUR: ", comma(suggestive_n),
                    "\n• Non-EUR ancestry-enriched: ", comma(enriched_n)), 
             ''),
    marker = list(colors = colors, line = list(color = colors)),
    textfont = list(family = "Arial", size = 40, color = "black"),
    opacity = 1,
    width = 500, height = 400) %>%
    layout(autosize = F, margin = m)
  p
  return(p)
}

p3a = meta_eur_freq_fc('max', threshold = 5e-8)
p3b = funnel_plot()
camk2d = image_read_pdf("figure3c_CAMK2D_4_114604622_META_500000_biomarkers-30870-both_sexes--irnt.pdf")
p3c = ggplot_pdf(camk2d)
pitx = image_read_pdf("figure3d_PITX2_4_112282681_META_650000_continuous-5099-both_sexes--irnt.pdf")
p3d = ggplot_pdf(pitx)

figure3 = function(output_format = 'png', create_subplots=F) {
  nrows = 4
  height = 5.5 * nrows / 2
  width = 6.5 / 2
  scale = 1.1
  pdf('figure3_panel1.pdf', height=scale * height / nrows, width=scale * width)
  print(p3a)
  dev.off()
  save_image(p3b, 'figure3_panel2.pdf') #, height=100*height/2, width=100*width/2)
  p3a1 = ggplot_pdf(image_read_pdf('figure3_panel1.pdf'))
  p3b1 = ggplot_pdf(image_read_pdf('figure3_panel2.pdf'))
  
  output_type(output_format, paste0('figure3.', output_format), height=height, width=width)
  print(ggarrange(p3a1, p3b1, p3c, p3d, nrow = nrows, ncol = 1, labels='auto'))
  # output_type(output_format, paste0('figure3.', output_format), height=5.27, width=3.44*2)
  # right = ggarrange(p3b + plot_annotation(tag_levels = list(c('b', ''))),
  #                   p3c + plot_annotation(tag_levels = list(c('c', ''))), nrow=2
  # )
  # p = ggarrange(p3a + plot_annotation(tag_levels = list(c('a', ''))), right, ncol=2)
  # print(p)
  dev.off()
}
figure3()
figure3('pdf')


# novel_hits %>%
#   rename(gene=nearest_genes) %>%
#   separate_rows(gene) %>%
#   left_join(gene_lists, multiple = 'all') %>%
#   filter(grepl('haploinsufficiency', gene_list)) %>%
#   View

# novel_hits %>% count(trait_efo_category)

# novel_hits %>% 
#   mutate(freq_max = pmax(AFR_af, AMR_af, CSA_af, EAS_af, MID_af),
#          fc = freq_max / EUR_af) %>% 
#   filter(fc > 2) %>%
#   filter(!is.na(trait_efo_category) & trait_efo_category != 'Other measurement') %>%
#   filter(trait_type %in% c('icd10', 'phecode')) %>%
#   arrange(desc(freq_max)) %>% View

# novel_hits %>% 
#   filter(EUR_af < 0.05 & AFR_af > 0.05) %>%
#   filter(!is.na(trait_efo_category) & trait_efo_category != 'Other measurement') %>%
#   count(trait_type, phenocode, pheno_sex, coding, modifier) %>% count
