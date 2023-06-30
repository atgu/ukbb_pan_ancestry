source('~/ukbb_pan_ancestry/constants.R')
library(gganimate)
library(gifski)

theme_set(theme_classic())

threshold = 1e-10
log_threshold = -log10(threshold)
meta_only_color = muted('green')
both_color = muted('purple')
eur_only_color = muted(color_eur)

meta_eur_comparison = read_tsv(gzfile('data/meta_eur_comparison.tsv.bgz'))
meta_eur_comparison %>% filter(N_pops > 1) %>% mutate(
  # nlog10p(meta) - nlog10p(EUR) = log10(P[EUR]/P[Meta])
  ratio = Pvalue_meta - Pvalue_EUR,
  freq_max = pmax(freq_AFR, freq_AMR, freq_CSA, freq_EAS, freq_MID),
  sig_in=case_when(Pvalue_meta > log_threshold & Pvalue_EUR < log_threshold ~ 'meta_only',
                   Pvalue_meta < log_threshold & Pvalue_EUR > log_threshold ~ 'EUR_only',
                   TRUE ~ 'both')
) %>% filter(Pvalue_meta > log_threshold | Pvalue_EUR > log_threshold) -> meta_eur_comparison

meta_eur_comparison %>%
  filter(sig_in == 'meta_only') %>%
  count(phenocode, nearest_genes, description) %>%
  count(nearest_genes) %>%
  arrange(desc(n)) %>%
  print(n=20)

meta_eur_density_plot = function(annotate_zones=F) {
  p = meta_eur_comparison %>% head %>%
    ggplot
  right = max(meta_eur_comparison$Pvalue_EUR, na.rm=T)
  top = max(meta_eur_comparison$Pvalue_meta, na.rm=T)
  if (annotate_zones) {
    box_alpha = 0.2
    p = p +
      geom_rect(aes(xmin=0, xmax=log_threshold, ymin=log_threshold, ymax=100), fill=meta_only_color, alpha=box_alpha) +
      geom_rect(aes(xmin=log_threshold, xmax=right, ymin=log_threshold, ymax=top), fill=both_color, alpha=box_alpha) +
      geom_rect(aes(xmin=log_threshold, xmax=100, ymin=0, ymax=log_threshold), fill=eur_only_color, alpha=box_alpha)
  }
  p + geom_hex(aes(x = Pvalue_EUR, y = Pvalue_meta),
             meta_eur_comparison %>% filter(Pvalue_meta > log_threshold | Pvalue_EUR > log_threshold),
             bins=100) +
    scale_x_continuous(trans = locusviz::trans_loglog_p()) +
    scale_y_continuous(trans = locusviz::trans_loglog_p()) +
    geom_segment(x = 0, xend = log_threshold, y = log_threshold, yend = log_threshold, color = "black", linetype='dashed') +
    geom_segment(x = log_threshold, xend = log_threshold, y = 0, yend = log_threshold, color = "black", linetype='dashed') +
    labs(x = "-log10 P (EUR)", y = "-log10 P (Meta-analysis)") +
    scale_fill_gradient(low = low_color, high = high_color, labels=comma, name='Number of\ntrait-loci pairs') +
    theme(legend.position = c(0.01, 1), legend.justification = c(0, 1),
          legend.background=element_rect(fill = alpha("white", 0))) -> p
  return(p)
}

meta_eur_freq_fc = function(pop = 'CSA') {
  plot_data = meta_eur_comparison 
  if(pop == 'max') {
    high_col = '#008080'
    low_col = 'black'
    label_name = expression(paste(log[10], '(', frac(max~AF, EUR~AF), ')'))
  } else {
    high_col = ukb_pop_colors[[pop]]
    low_col = color_eur
    label_name = paste0('log10(', pop, ' AF / EUR AF)')
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
    geom_segment(x = 0, xend = log_threshold, y = log_threshold, yend = log_threshold, color = "black", linetype='dashed') +
    geom_segment(x = log_threshold, xend = log_threshold, y = 0, yend = log_threshold, color = "black", linetype='dashed') +
    labs(x = expression(paste(-log[10], "P (EUR)")), y = expression(paste(-log[10], "P (Meta-analysis)"))) +
    guides(color=guide_colorbar(title.position='top', title.hjust=1, direction = 'horizontal')) +
    theme(legend.position = c(1, 0), legend.justification = c(1, 0),
          legend.background=element_rect(fill = alpha("white", 0))) -> plt
  plt
  return(plt)
}
meta_eur_freq_fc('max')

meta_eur_heterogeneity = function() {
  plot_data = meta_eur_comparison
  
  high_col = muted('purple')
  ggplot(data=plot_data %>% head) +
    scattermore::geom_scattermore(
      aes(Pvalue_het, ratio, color = log10(.data[['freq_max']]/.data[['freq_EUR']])),
      pointsize = 2, pixels = c(1024, 1024), alpha = 0.5,
      data = plot_data
    ) +
    scale_x_continuous(trans = locusviz::trans_loglog_p()) +
    scale_color_gradient2(high = high_col, low = color_eur, mid = 'gray90', name=paste0('log10(max AF / EUR AF)')) +
    labs(x = "-log10 P (het)", y = "log-ratio meta / EUR") +
    theme(legend.position = c(0.01, 0.01), legend.justification = c(0, 0),
          legend.background=element_rect(fill = alpha("white", 0))) -> plt
  # plt
  return(plt)
}

prepare_table = function() {
  data_proc = meta_eur_comparison %>%
    mutate(mean_ldscore = (ldscore_EUR + ldscore_AFR + ldscore_AMR +
                             ldscore_CSA + ldscore_EAS + ldscore_MID)/6,
           mean_freq = (freq_EUR + freq_AFR + freq_AMR +
                          freq_CSA + freq_EAS + freq_MID)/6,
           common_non_EUR=freq_AFR > 0.01 | freq_AMR > 0.01 | freq_CSA > 0.01 | freq_EAS > 0.01 | freq_MID > 0.01,
           common_AFR = freq_AFR > 0.01,
           common_EUR = freq_EUR > 0.01)
  
  data_proc %>%
    group_by(sig_in) %>%
    summarize(
      `Number of variants (raw)`=n(),
      `Percent low INFO score` = 100 * sum(info < 0.9)/n(),
      `Percent low quality`=100 * sum(!high_quality)/n(),
      `Percent heterogeneous`=100 * sum(Pvalue_het > 2)/n()
      ) %>% gather(var, value, -sig_in) %>% 
    spread(sig_in, value) %>% arrange(var) -> pre_filter
  
  pre_filter %>% union_all(
    data_proc %>%
      filter(Pvalue_het < 2 & high_quality) %>%
      group_by(sig_in) %>%
      summarize(
        `Number of variants (filtered)`=n(),
        # `Percent common in EUR`=100 * sum(common_EUR)/n(),
        # `Percent common in AFR`=100 * sum(common_AFR)/n(),
        # `Percent common in AFR not EUR`=100 * sum(common_AFR & !common_EUR)/n(),
        `Percent common outside EUR (rare in EUR)`=100 * sum(common_non_EUR & !common_EUR)/n(),
        # mean_ldscore_EUR=mean(ldscore_EUR, na.rm=T),
        # mean_ldscore_AFR=mean(ldscore_AFR, na.rm=T),
        # mean_ldscore_resid_EUR=mean(ldscore_resid_EUR, na.rm=T),
        # mean_ldscore_resid_AFR=mean(ldscore_resid_AFR, na.rm=T),
        # mean_ldscore_diff_EUR_AFR=mean(ldscore_resid_EUR - ldscore_resid_AFR, na.rm=T)
        ) %>% gather(var, value, -sig_in) %>% 
      spread(sig_in, value) %>% arrange(var)
    ) %>% select(var, meta_only, both, EUR_only) -> summaries
  
  summaries = summaries %>% mutate(across(where(is.numeric), ~ gsub('\\.0$', '', sprintf("%.1f", .x))))
  summaries = union_all(data.frame(var='', meta_only='Meta-only', both='Both', EUR_only='EUR-only'), summaries) %>% tibble
  return(summaries)
}

generate_table = function() {
  summaries = prepare_table()
  header_cell_color = 'white'
  thm = ttheme(
    colnames.style = colnames_style(color = "white", fill = "#8cc257"),
    tbody.style = tbody_style(color = "black", fill = 'white',
                              hjust = as.vector(matrix(c(0, 1, 1, 1), ncol = 4, nrow = nrow(summaries), byrow = TRUE)),
                              x = as.vector(matrix(c(0, 0.9, 0.9, 0.9), ncol = 4, nrow = nrow(summaries), byrow = TRUE)))
  )
  ggtab <- ggtexttable(summaries, rows=NULL, cols=NULL, theme = thm)
  ggtab <- table_cell_bg(ggtab, row = 1, column = 2, fill = meta_only_color, color=header_cell_color)
  ggtab <- table_cell_bg(ggtab, row = 1, column = 3, fill = both_color, color=header_cell_color)
  ggtab <- table_cell_bg(ggtab, row = 1, column = 4, fill = eur_only_color, color=header_cell_color)
  for (i in 1:3) {
    ggtab <- table_cell_font(ggtab, row = 1, column = i+1, color=header_cell_color)
  }
  return(ggtab)
}

power_curve = function(return_plots=F, point_size = 1) {
  staging = meta_eur_comparison %>%
    filter(modifier == 'irnt') %>%
    mutate(beta_EUR = abs(beta_EUR),
           beta_meta = abs(beta_meta),
           MAF_EUR = 0.5 - abs(0.5 - freq_EUR),
           MAF_max = 0.5 - abs(0.5 - freq_max))
  
  staging %>% count(trait_type, phenocode, pheno_sex, coding, modifier) %>% count
  
  set.seed(663)
  plot_data = staging %>%
    filter(sig_in != 'meta_only') %>%
    sample_frac(0.01) %>%
    union_all(staging %>% filter(sig_in == 'meta_only') %>% sample_frac(0.5))
  
  settings = list(ylim(0, 1.15), xlab('MAF'), ylab('Beta'),
                  theme(legend.position = c(1, 1), legend.justification = c(1, 1)))
  p1 = plot_data %>%
    filter(sig_in != 'meta_only') %>%
    ggplot + aes(x = MAF_EUR, y = beta_EUR) +
    geom_point(color=color_eur, size = point_size) +
    settings + xlab('EUR MAF') + ylab('EUR beta')

  p2 = plot_data %>%
    filter(sig_in != 'meta_only') %>%
    ggplot + aes(x = MAF_EUR, y = beta_EUR) +
    geom_point(color=color_eur, size = point_size) +
    geom_point(aes(x = MAF_EUR, y = beta_EUR, color=max_maf_pop), size = point_size,
               data = plot_data %>% filter(sig_in == 'meta_only' & Pvalue_meta > 12)) +
    settings + xlab('EUR MAF') + ylab('EUR beta')
  
  p2_black = plot_data %>%
    filter(sig_in != 'meta_only') %>%
    ggplot + aes(x = MAF_EUR, y = beta_EUR, color=color_eur) +
    geom_point(size = point_size) + 
    geom_point(aes(x = MAF_EUR, y = beta_EUR, color='black'), size = point_size,
               data = plot_data %>% filter(sig_in == 'meta_only' & Pvalue_meta > 12)) +
    settings + xlab('EUR MAF') + ylab('EUR beta') + scale_color_identity()
  
  p3 = plot_data %>%
    filter(sig_in != 'meta_only') %>%
    ggplot + aes(x = MAF_EUR, y = beta_EUR) +
    geom_point(color=color_eur, size = point_size) +
    geom_point(aes(x = MAF_max, y = beta_meta, color=max_maf_pop), size = point_size,
               data = plot_data %>% filter(sig_in == 'meta_only' & Pvalue_meta > 12)) +
    settings + xlab('Ancestry-specific MAF') + ylab('Meta-analysis beta') + pop_color_scale
  p4 = plot_data %>%
    filter(sig_in != 'meta_only') %>%
    ggplot + aes(x = MAF_EUR, y = beta_EUR) +
    geom_point(color=color_eur, size = point_size) +
    geom_segment(aes(x = MAF_EUR, y = beta_EUR, xend=MAF_max, yend = beta_meta, color=max_maf_pop),
                 arrow=arrow(length=unit(0.05, "inches")),
                 data = plot_data %>% filter(sig_in == 'meta_only' & Pvalue_meta > 12)) +
    settings + xlab('Ancestry-specific MAF') + ylab('Meta-analysis beta') + pop_color_scale
  
  if (return_plots) {
    return(list(p1, p2, p2_black, p3, p4))
  }
  width = 4.5
  height = 3.5
  png('EUR_discovery.png', width = width, height = height, res = 300, units='in')
  print(p1)
  dev.off()
  png('meta_discovery_in_EUR.png', width = width, height = height, res = 300, units='in')
  print(p2)
  dev.off()
  png('meta_discovery_in_EUR_anon.png', width = width, height = height, res = 300, units='in')
  print(p2_black)
  dev.off()
  png('meta_discovery_in_meta.png', width = width, height = height, res = 300, units='in')
  print(p3)
  dev.off()
  
  before_state = 1
  after_state = 2
  test = plot_data %>% filter(sig_in != 'meta_only') %>% mutate(state=before_state, max_maf_pop = 'EUR') %>%
    union_all(plot_data %>% filter(sig_in != 'meta_only') %>% mutate(state=after_state, max_maf_pop = 'EUR')) %>%
    union_all(plot_data %>% filter(sig_in == 'meta_only' & Pvalue_meta > 12) %>% mutate(state=before_state)) %>%
    union_all(plot_data %>% filter(sig_in == 'meta_only' & Pvalue_meta > 12) %>% mutate(state=after_state) %>%
                mutate(MAF_EUR = MAF_max, beta_EUR = beta_meta))

  # possible to change from EUR color to other color by creating more states and interpolating manually
  xlabel = c('1' = 'EUR MAF', '2' = 'Ancestry-specific MAF')
  ylabel = c('1' = 'EUR beta', '2' = 'Meta-analysis beta')
  p = test %>%
    ggplot + aes(x = MAF_EUR, y = beta_EUR, color=max_maf_pop) +
    geom_point() +
    transition_states(state, wrap = F) +
    settings + 
    ease_aes('quadratic-out') +
    labs(x = '{xlabel[[closest_state]]}', y = '{ylabel[[closest_state]]}')
  
  ren = animate(p, duration = 3, fps = 10, width = width, height = height, res = 300, units='in',
                renderer = gifski_renderer(loop = FALSE), rewind=FALSE)
  anim_save("meta_discovery.gif", ren)
  ren_loop = animate(p, duration = 3, fps = 10, width = width, height = height, res = 300, units='in',
                renderer = gifski_renderer(loop = TRUE), rewind=FALSE)
  anim_save("meta_discovery_loop.gif", ren_loop)
}

efo = read_tsv(gzfile('data/known_ukbb_loci.meta_hq_annotated.txt.bgz'),
               col_types = cols(locus=col_character(),
                                lead_locus=col_character(),
                                locus_otg=col_character()))
known_novel_efo = function() {
  n_known = efo %>% 
    filter(trait_efo_category == trait_efo_category_otg) %>%
    distinct(trait_efo_category, idx) 
  # 66,280
  
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
  return(p)
}

power_plots = power_curve(return_plots=T, point_size = 0.8)
p3a = power_plots[[3]] +
  scale_color_identity(breaks = c("black", color_eur),
                       labels = c("Discovered only in meta-analysis", "Discovered only in EUR"),
                       guide='legend') +
  theme(legend.key.size = unit(0.05, 'in'), legend.position=c(1, 1.25)) +
  guides(color=guide_legend(nrow=2, byrow=TRUE, title=NULL))
p3b = power_plots[[4]] +
  theme(legend.key.size = unit(0.05, 'in'), legend.position=c(1, 1.25)) +
  guides(color=guide_legend(nrow=2, byrow=TRUE, title=NULL))
p3c = meta_eur_freq_fc('max')
p3low = known_novel_efo()
figure3 = function(output_format = 'png', create_subplots=F) {
  # table_size = if_else(create_subplots, 4, 2)
  
  if (create_subplots) {
    fig3_plots = list(p3a, p3b, p3c, p3low)
    for (i in seq_along(fig3_plots)) {
      output_type(output_format, paste0('figure3_panel', i, '.', output_format),
                  height=4, width=6)
      print(fig3_plots[i])
      dev.off()
    }
  } else {
    output_type(output_format, paste0('figure3.', output_format), height=6, width=7.5)
    # print(ggarrange(p1a, p1b, p1c, p1d, nrow = 2, ncol = 2))
    print(((p3a / p3b) | p3c) / p3low + plot_annotation(tag_levels = 'a'))
    dev.off()
  }
}
figure3()
figure3(create_subplots = T)

output_type('png', paste0('figure3_panel', i, '.png'), height=4, width=6)
print(p3a2)
dev.off()

quadrant_plot = function() {
  p3c = meta_eur_density_plot()
  sp3a2 = meta_eur_density_plot(annotate_zones=T)
  sp3b = generate_table()
}

output_type('png', 's_figure_X_het.png', height=4, width=6)
print(meta_eur_heterogeneity())
dev.off()

beta_correlation_by_p_value = function() {
  plot_color = 'royalblue'
  
  all_pairwise_data = load_ukb_file('pairwise_beta_corr_full.txt.bgz', 'effect_size/',
                                    force_cols=cols(pop1.n_controls=col_integer(),
                                                    pop2.n_controls=col_integer())) %>%
    filter(pop1.pop != pop2.pop)
  pairwise_data = all_pairwise_data %>% filter(pop1.pop > pop2.pop)
  
  plot_data = pairwise_data %>%
    filter(pairwise_corr.n > 20) %>%
    group_by(p_value_threshold) %>%
    summarize(sig_pos=sum(beta_p < 0.05 & beta > 0, na.rm=T),
              sig=sum(beta_p < 0.05, na.rm=T),
              prop_sig_pos = sig_pos / sig,
              sem = 1.96 * sqrt(prop_sig_pos * (1 - prop_sig_pos) / sig)) %>% 
    ungroup
  
  low_point = 5
  p = plot_data %>%
    ggplot + aes(x = p_value_threshold, y = prop_sig_pos,
                 ymin = prop_sig_pos - sem, ymax = prop_sig_pos + sem) + 
    geom_pointrange(color=plot_color) +
    scale_y_continuous(labels=percent, name='Proportion of betas positively correlated',
                       breaks=low_point:10/10, limits=c(low_point / 10, 1)) +
    geom_hline(yintercept=0.5, linetype='dashed') +
    scale_x_log10(name='p-value threshold in both populations')
  return(p)
}



chopping_block = function() {
  polygenicity = read_tsv('data/230123_sbrs_polygenicity_summary.txt')
  
  size = 1000000
  data_proc %>%
    filter(Pvalue_het < 2 & high_quality & sig_in == 'meta_only') %>%
    mutate(window=size * floor(pos / size)) %>%
    unite(col='Pheno', trait_type, phenocode, pheno_sex, coding, modifier, na.rm=T, remove=F, sep='-') %>%
    count(chrom, window, Pheno, trait_type, phenocode, pheno_sex, coding, modifier, description) %>%
    count(Pheno, trait_type, phenocode, pheno_sex, coding, modifier, description) %>%
    filter(n >= 10) %>%
    arrange(n) %>%
    left_join(polygenicity) %>%
    filter(R_GelmanRubin < 1.2) %>% 
    arrange(Mean)
  
}