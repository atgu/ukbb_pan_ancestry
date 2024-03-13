source('~/ukbb_pan_ancestry/constants.R')

h2_data = load_ukb_file('h2_estimates_all_flat_221123.tsv', parent_folder='h2/')

h2_by_group = function() {
  p = h2_data %>%
    # filter(qcflags.in_bounds_h2) %>%
    filter(qcflags.pass_all) %>%
    mutate(ancestry=fct_relevel(ancestry, pops_by_sample_size)) %>%
    group_by(ancestry, trait_type) %>%
    summarize(mean_observed_h2=mean(estimates.final.h2_observed, na.rm=T)) %>%
    mutate(trait_type=fct_reorder(trait_type, -mean_observed_h2)) %>%
    ggplot + aes(x = ancestry, fill = trait_type, y = mean_observed_h2) + 
    geom_bar(stat='identity', position='dodge') +
    trait_fill_scale + xlab(NULL) +
    scale_y_continuous(labels=percent, name=expression(paste('Mean observed ', h^2))) +
    theme(legend.position=c(0.01, 1), legend.justification=c(0, 1),
          legend.background=element_rect(fill = alpha("white", 0)))
  return(p)
}
h2_z_by_group = function() {
  high_level = 8
  h2_plot_data = h2_data %>%
    filter(qcflags.in_bounds_h2 & qcflags.normal_lambda & qcflags.normal_ratio) %>%
    mutate(z=if_else(estimates.final.h2_z > high_level, high_level, estimates.final.h2_z)) %>%
    filter(z > 0) %>%
    mutate(ancestry=fct_relevel(ancestry, pops_by_sample_size)) 
  h2_medians = h2_plot_data %>%
    group_by(ancestry, trait_type) %>%
    summarize(median_z=median(z))
  colors_darker = darken(ukb_pop_colors, amount=0.375)
  names(colors_darker) = names(ukb_pop_colors)
  p = h2_plot_data %>%
    ggplot + aes(x = z, fill = ancestry) + 
    geom_histogram() + geom_vline(xintercept=4, linetype='dashed') +
    geom_vline(aes(xintercept = median_z, color=ancestry), data=h2_medians, linetype='dotted') +
    coord_cartesian(ylim=c(NA, 250)) +
    scale_y_continuous(breaks=c(0, 125, 250)) + 
    # scale_x_continuous(breaks=c(0, 10, 20)) +
    scale_x_continuous(breaks=c(0, 4, 8)) +
    pop_fill_scale + scale_color_manual(values=colors_darker) +
    facet_grid(rows=vars(ancestry), cols=vars(trait_type)) +
    theme(strip.background = element_blank()) +
    guides(fill=F, color=F) + ylab('Number of phenotypes') + xlab('Heritability z')
  return(p)
}

p1 = h2_z_by_group()
p2 = h2_by_group()
extended_data_figure3 = function(output_format = 'png') {
  output_type(output_format, paste0('extended_data_figure3.', output_format), height=8.5, width=5.5)
  print(ggarrange(p1, p2, nrow = 2, labels='auto'))
  dev.off()
}
extended_data_figure3()
