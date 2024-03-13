source('figure3.R')

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
    scale_x_continuous(trans = locusviz::trans_loglog_p(), label=comma) +
    scale_y_continuous(trans = locusviz::trans_loglog_p(), label=comma) +
    geom_segment(x = 0, xend = log_threshold, y = log_threshold, yend = log_threshold, color = "black", linetype='dashed') +
    geom_segment(x = log_threshold, xend = log_threshold, y = 0, yend = log_threshold, color = "black", linetype='dashed') +
    labs(x = expression(paste(-log[10], "P (EUR)")), y = expression(paste(-log[10], "P (Meta-analysis)"))) +
    scale_fill_gradient(low = low_color, high = high_color, labels=comma, name='Number of\ntrait-loci pairs') +
    theme(legend.position = c(0.01, 1), legend.justification = c(0, 1),
          legend.background=element_rect(fill = alpha("white", 0))) -> p
  return(p)
}
meta_eur_density_plot(annotate_zones=T)

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
        `Percent common in AFR (rare in EUR)`=100 * sum(common_AFR & !common_EUR)/n(),
        # mean_ldscore_EUR=mean(ldscore_EUR, na.rm=T),
        # mean_ldscore_AFR=mean(ldscore_AFR, na.rm=T),
        # mean_ldscore_resid_EUR=mean(ldscore_resid_EUR, na.rm=T),
        # mean_ldscore_resid_AFR=mean(ldscore_resid_AFR, na.rm=T),
        # mean_ldscore_diff_EUR_AFR=mean(ldscore_resid_EUR - ldscore_resid_AFR, na.rm=T)
      ) %>% gather(var, value, -sig_in) %>% 
      spread(sig_in, value) %>% arrange(var)
  ) %>% select(var, meta_only, both, EUR_only) -> summaries
  
  summaries = summaries %>% mutate(across(where(is.numeric), ~ gsub('\\.00$', '', sprintf("%.2f", .x))))
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

p1 = meta_eur_density_plot(annotate_zones=T)
p2 = generate_table()

extended_data_figure7 = function(output_format = 'png') {
  output_type(output_format, paste0('extended_data_figure7.', output_format), height=6, width=13)
  print(ggarrange(p1, p2, ncol = 2, labels='auto'))
  dev.off()
}
extended_data_figure7()


