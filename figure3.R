source('~/ukbb_pan_ancestry/constants.R')
library(gganimate)
library(gifski)

theme_set(theme_classic())

threshold = 1e-10
log_threshold = -log10(threshold)
orig_threshold = 5e-8
log_orig_threshold = -log10(orig_threshold)
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
p1 = plot_n_sig_per_pheno(F)
p2 = plot_n_sig_per_pheno()

print(p1 | p2)

meta_eur_freq_fc = function(pop = 'CSA', threshold=1e-10, clump_only=F) {
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
    geom_segment(x = 0, xend = log_threshold, y = log_threshold, yend = log_threshold, color = "black", linetype='dashed') +
    geom_segment(x = log_threshold, xend = log_threshold, y = 0, yend = log_threshold, color = "black", linetype='dashed') +
    labs(x = expression(paste(-log[10], "P (EUR)")), y = expression(paste(-log[10], "P (Meta-analysis)"))) +
    guides(color=guide_colorbar(title.position='top', title.hjust=1, direction = 'horizontal')) +
    theme(legend.position = c(1, 0), legend.justification = c(1, 0),
          legend.background=element_rect(fill = alpha("white", 0))) -> plt
  plt
  return(plt)
}

power_curve = function(return_plots=F, point_size = 1, EUR_fraction_to_sample = 0.01, meta_fraction_to_sample=0.5, clump_only=F) {
  staging = meta_eur_comparison %>%
    filter(modifier == 'irnt') %>%
    mutate(beta_EUR = abs(beta_EUR),
           beta_meta = abs(beta_meta),
           MAF_EUR = 0.5 - abs(0.5 - freq_EUR),
           MAF_max = 0.5 - abs(0.5 - freq_max))
  
  if (clump_only) {
    staging %<>% filter(!is.na(in_clump))
    EUR_fraction_to_sample = 1
    meta_fraction_to_sample = 1
  }
  
  staging %>% count(trait_type, phenocode, pheno_sex, coding, modifier) %>% count
  
  set.seed(663)
  plot_data = staging %>%
    filter(sig_in %in% c('EUR_only', 'both')) %>%
    sample_frac(EUR_fraction_to_sample) %>%
    union_all(staging %>% filter(sig_in == 'meta_only') %>% sample_frac(meta_fraction_to_sample))
  
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
p = power_curve(return_plots=T, point_size = 0.8)
p2 = power_curve(return_plots=T, point_size = 0.8, clump_only=T)

figure3 = function(output_format = 'png', create_subplots=F, clump_only=F) {
  power_plots = power_curve(return_plots=T, point_size = 0.8, meta_fraction_to_sample = 1, clump_only=clump_only)
  p3a = power_plots[[3]] +
    scale_color_identity(breaks = c("black", color_eur),
                         labels = c("Discovered only in meta-analysis", "Discovered only in EUR"),
                         guide='legend') +
    theme(legend.key.size = unit(0.05, 'in'), legend.position=c(1, 1.25)) +
    guides(color=guide_legend(nrow=2, byrow=TRUE, title=NULL))
  p3b = power_plots[[4]] +
    theme(legend.key.size = unit(0.05, 'in'), legend.position=c(1, 1.25)) +
    guides(color=guide_legend(nrow=2, byrow=TRUE, title=NULL))
  p3c = meta_eur_freq_fc('max', clump_only=clump_only)
  
  # table_size = if_else(create_subplots, 4, 2)
  
  if (create_subplots) {
    fig3_plots = list(p3a, p3b, p3c)
    for (i in seq_along(fig3_plots)) {
      output_type(output_format, paste0('figure3_panel', i, if_else(clump_only, '.clump',), '.', output_format),
                  height=4, width=6)
      print(fig3_plots[i])
      dev.off()
    }
  } else {
    output_type(output_format, paste0('figure3.', if_else(clump_only, 'clump.', ''), output_format), height=3, width=7.5)
    # print(ggarrange(p1a, p1b, p1c, p1d, nrow = 2, ncol = 2))
    print(((p3a / p3b) | p3c) + plot_annotation(tag_levels = 'a'))
    dev.off()
  }
}
figure3()
figure3(clump_only=T)
figure3(create_subplots = T)
