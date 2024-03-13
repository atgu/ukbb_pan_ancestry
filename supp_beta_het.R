source('figure3.R')

meta_eur_heterogeneity = function(threshold=1e-10) {
  plot_data = meta_eur_comparison %>%
    filter(Pvalue_meta > -log10(threshold) | Pvalue_EUR > -log10(threshold)) 
  
  high_col = muted('limegreen', l=50, c = 70)
  ggplot(data=plot_data %>% head) +
    scattermore::geom_scattermore(
      aes(Pvalue_EUR, Pvalue_meta, color = Pvalue_het),
      pointsize = 4, pixels = c(1024, 1024), alpha = 1,
      data = plot_data # %>% head(1000)
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    scale_x_continuous(trans = locusviz::trans_loglog_p(), labels=comma) +
    scale_y_continuous(trans = locusviz::trans_loglog_p(), labels=comma) +
    scale_color_gradient(high = high_col, low = 'gray90',
                          name=expression(paste(log[10], '(het P)'))) +
    geom_segment(x = 0, xend = log_threshold, y = log_threshold, yend = log_threshold, color = "black", linetype='dashed') +
    geom_segment(x = log_threshold, xend = log_threshold, y = 0, yend = log_threshold, color = "black", linetype='dashed') +
    labs(x = expression(paste(-log[10], "P (EUR)")), y = expression(paste(-log[10], "P (Meta-analysis)"))) +
    guides(color=guide_colorbar(title.position='top', title.hjust=1, direction = 'horizontal')) +
    theme(legend.position = c(1, 0), legend.justification = c(1, 0),
          legend.background=element_rect(fill = alpha("white", 0))) -> plt
  return(plt)
}

all_pairwise_data = load_ukb_file('pairwise_beta_corr_full.txt.bgz', 'effect_size/',
                                  force_cols=cols(pop1.n_controls=col_integer(),
                                                  pop2.n_controls=col_integer())) %>%
  filter(pop1.pop != pop2.pop)
pairwise_data = all_pairwise_data %>% filter(pop1.pop > pop2.pop)

beta_correlation_by_p_value = function() {
  plot_color = 'royalblue'
    
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

s_figure_32 = function(output_format = 'png') {
  p1 = meta_eur_heterogeneity()
  p2 = beta_correlation_by_p_value()
  output_type(output_format, paste0('s_figure_32.', output_format), height=3.75, width=7.5)
  print((p1 | p2) + plot_annotation(tag_levels = 'a'))
  dev.off()
}
s_figure_32()
