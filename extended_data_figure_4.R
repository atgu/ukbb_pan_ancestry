source('figure3.R')
library(ggExtra)

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

plot_n_sig_per_pheno_scatter = function(log_scale=F) {
  p = meta_eur_comparison %>%
    filter(Pvalue_meta > -log10(5e-8) & high_quality) %>%
    group_by_at(key_fields) %>%
    summarize(clumped=sum(!is.na(in_clump)), n=n()) %>%
    ggplot + aes(x = n, y = clumped, color = trait_type) +
    geom_point() + 
    trait_color_scale + 
    theme(legend.position=c(0.01, 0.99),
          legend.justification=c(0, 1))
    
  xlabel = 'Number of significant associations'
  ylabel = 'Number of independent significant associations'
  if (log_scale) {
    p + scale_x_log10(labels=comma, name=xlabel) +
      scale_y_log10(labels=comma, name=ylabel) %>%
      return
  } else {
    p + scale_x_continuous(labels=comma, name=xlabel) +
      scale_y_continuous(labels=comma, name=ylabel) %>%
    return
  }
}

n_variants_by_group_plot = function() {
  n_variants = tribble(
    ~ancestry, ~n_variants,
    'EUR', 9956319,
    'CSA', 10939516,
    'AFR', 17826645,
    'EAS', 18082070,
    'MID', 18256015,
    'AMR', 18296910,
  ) %>%
    mutate(ancestry=fct_relevel(ancestry, pops_by_sample_size))

  plus_pops = c('EUR', paste0('+', pops_by_sample_size)[pops_by_sample_size != 'EUR'])
  names(plus_pops) = pops_by_sample_size
  p = n_variants %>%
    ggplot + aes(x = ancestry, y = n_variants, fill = ancestry) + 
    geom_bar(stat='identity') +
    scale_y_continuous(labels=comma, name='Number of variants with AF > 1%\nin at least one genetic ancestry group') +
    scale_fill_manual(values=ukb_pop_colors, guide=F) +
    scale_x_discrete(labels=plus_pops, name=NULL)
  return(p)
}

meta_eur_comparison %>%
  filter(freq_EUR < 0.01 & sig_in != 'EUR_only') %>%
  count(chrom, pos, ref, alt) %>%
  count

extended_data_figure4 = function(output_format = 'png') {
  output_type(output_format, paste0('extended_data_figure4.', output_format), height=4.5, width=8.5)
  p1 = plot_n_sig_per_pheno_scatter(T) %>%
    ggMarginal(type='histogram', bins=50)
  p2 = n_variants_by_group_plot()
  print(ggarrange(p1, p2, ncol = 2, labels='auto'))
  dev.off()
}
extended_data_figure4()
