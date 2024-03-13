source('constants.R')

polygenicity = read_tsv('data/230123_sbrs_polygenicity_summary.txt')

polygenicity %>% 
  count(R_GelmanRubin < 1.2)

p = polygenicity %>%
  filter(R_GelmanRubin < 1.2) %>% 
  mutate(trait_type = str_split_i(Pheno, '-', 1)) %>%
  ggplot + aes(x = Mean, fill = trait_type) +
  geom_histogram() + facet_wrap(~ trait_type) +
  trait_fill_scale + guides(fill=F) +
  xlab('Polygenicity') + ylab('Number of phenotypes')

s_figure_31 = function(output_format = 'png', create_subplots=F) {
  output_type(output_format, paste0('s_figure_31.', output_format), height=5, width=7.5)
  print(p)
  dev.off()
}
s_figure_31()
