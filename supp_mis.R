source('constants.R')

pairwise_data = read_tsv(gzfile('data/pairwise_correlations_regressed.hq_phenos.txt.bgz'))

pairwise_hist = function(filter_r = 1) {
  pairwise_data %>%
    mutate(entry=if_else(abs(entry) <= filter_r, entry, sign(entry)*filter_r)) %>%
    filter(i < j) %>%
    ggplot + aes(x = entry) +
    geom_histogram(bins=100) + 
    scale_y_continuous(labels=comma, name='Pairs of phenotypes') +
    xlab('Correlation (r)') %>%
    return
}

s_figure_29 = function(output_format = 'png', create_subplots=F) {
  p1 = pairwise_hist()
  p2 = pairwise_hist(sqrt(0.1))
  output_type(output_format, paste0('s_figure_29.', output_format), height=2.5, width=7.5)
  print(ggarrange(p1, p2, nrow = 1, ncol = 2, labels="auto"))
  dev.off()
}
s_figure_29()
