source('constants.R')

ld_score_comparison = function(anc) {
  ld_score = read_tsv(gzfile(paste0('data/ukb_gnomad_ld.', anc, '.txt.bgz')))
  ld_score %>%
    ggplot + aes(x = gnomad_ld_score, y = ld_score) +
    scattermore::geom_scattermore(
      pointsize = 2, pixels = c(1024, 1024), alpha = 0.5,
      color = ukb_pop_colors[[anc]]
    ) + scale_x_continuous(label=comma, name='gnomAD LD score', breaks=1000*(0:4)) +
    scale_y_continuous(label=comma, name='UKB LD score') + 
    geom_abline(slope=1, intercept=0, linetype='dashed') %>%
    return
}

plots = map(pops_with_gnomad_data, ld_score_comparison)

s_figure_20 = function(output_format = 'png', create_subplots=F) {
  width = 9
  output_type(output_format, paste0('s_figure_20.', output_format), height=width/4, width=width)
  print(ggarrange(plotlist = plots, nrow=1, ncol=4, labels = 'auto'))
  dev.off()
}
s_figure_20()
