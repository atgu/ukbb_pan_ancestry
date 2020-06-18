source('~/ukbb_pan_ancestry/constants.R')

ukb_gnomad_compare_options = list(xlab('gnomAD frequency'),
                                  ylab('UKB frequency'),
                                  scale_x_log10(label=pretty_axis_format),
                                  scale_y_log10(label=pretty_axis_format))
ukb_gnomad_ratio_name = 'UKB / gnomAD frequency ratio'
output_format = 'png'

load_all_ukb_gnomad_2d_hist_data = function() {
  map_df(pops_with_gnomad_data, function(x) load_ukb_file(paste0('gnomad_af_v_ukb_af_hist2d_', x, '.txt.bgz'), 'gnomad_comparison/') %>% mutate(pop=x)) %>%
    mutate(gnomad_freq = 10 ^ gnomad_freq, ukb_freq = 10 ^ ukb_freq) %>%
    return
}

data_2dhist = load_all_ukb_gnomad_2d_hist_data()

plot_ukb_gnomad_2dhist = function(pop_to_use, transformation = 'log10', save_plot = F, breaks=NULL) {
  if (save_plot) output_type(output_format, paste0('plots/', 'gnomad_vs_ukb_freq_', pop_to_use, '_', transformation, '.', output_format), height=4, width=5)
  if (is.null(breaks)) breaks = waiver()
  p = data_2dhist %>%
    filter(pop == pop_to_use) %>%
    ggplot + aes(x = gnomad_freq, y = ukb_freq, fill = n_variants) + 
    geom_tile() +
    scale_fill_continuous(low = low_color, high = high_color, name='Number of\nvariants', trans=transformation, breaks=breaks) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dotted') +
    ukb_gnomad_compare_options
  if (save_plot) {
    print(p)
    dev.off()
  }
  return(p)
}

map(pops_with_gnomad_data, function(x) plot_ukb_gnomad_2dhist(x, 'identity', TRUE))
map(pops_with_gnomad_data, function(x) plot_ukb_gnomad_2dhist(x, 'log10', TRUE))

# UKB estimate of gnomAD singleton frequency
plot_gnomad_singleton_freq_in_ukb = function(save_plot = F) {
  cut_point = 0.0011
  if (save_plot) output_type(output_format, paste0('plots/gnomad_singletons_in_ukb.', output_format), height=4, width=5.5)
  min_freq = data_2dhist %>% filter(pop == 'EUR') %$% min(gnomad_freq)
  p = data_2dhist %>% 
    filter(gnomad_freq == min_freq & pop == 'EUR') %>% 
    group_by(ukb_freq = if_else(ukb_freq > cut_point, cut_point, ukb_freq)) %>%
    summarize(n_variants = sum(n_variants)) %>%
    ggplot + aes(x = ukb_freq, y = n_variants) + 
    geom_bar(stat='identity', width=0.07, fill = color_nfe) + 
    geom_vline(xintercept = 1e-4, linetype = 'dashed') +
    scale_x_log10(name='UKB frequency estimate of gnomAD genomes singletons (0.01%)', label=pretty_axis_format,
                  breaks=c(2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2)
                  ) +
    scale_y_continuous(name='Number of variants', label=comma)
  if (save_plot) {
    print(p)
    dev.off()
  }
  return(p)
}
plot_gnomad_singleton_freq_in_ukb(T)
  

load_all_ukb_gnomad_scatter_data = function() {
  map_df(pops_with_gnomad_data, function(x) load_ukb_file(paste0('gnomad_af_v_ukb_af_scatter_', x, '.txt.bgz'), 'gnomad_comparison/',
                                         force_cols = cols(chrom=col_character())) %>% mutate(pop=x)) %>%
    mutate(gnomad_freq = 10 ^ gnomad_freq, ukb_freq = 10 ^ ukb_freq) %>%
    return
}
data_scatter = load_all_ukb_gnomad_scatter_data()

plot_off_diagonals = function(pop_to_use, color_by='pass_gnomad_genomes', n_points = 5000, true_first = TRUE, 
                              color_guide_title = '', save_plot=F) {
  if (save_plot) output_type(output_format, paste0('plots/', 'scatter_gnomad_vs_ukb_freq_', pop_to_use, '_by_', color_by, '.', output_format), height=4, width=5)
  
  data_scatter %>%
    filter(pop == pop_to_use) %>%
    group_by_at(vars(color_by)) %>%
    summarize(n_variants=n())
  
  dataset = data_scatter %>%
    filter(pop == pop_to_use) %>%
    filter_at(vars(color_by), any_vars(.)) %>%
    sample_n(n_points)
  
  not_dataset = data_scatter %>% 
    filter(pop == pop_to_use) %>%
    filter_at(vars(color_by), any_vars(!.)) %>% 
    sample_n(n_points)
  
  if (true_first) {
    plot_data = not_dataset %>%
      union_all(dataset)
  } else {
    plot_data = dataset %>%
      union_all(not_dataset)
  }
  p = plot_data %>%
    ggplot + aes(x = gnomad_freq, y = ukb_freq) + aes_string(color = color_by) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
    geom_point(alpha=0.5) + ukb_gnomad_compare_options
  
  if (color_guide_title != '') p = p + scale_color_few(name=color_guide_title, palette = 'Medium')
  
  if (save_plot) {
    print(p)
    dev.off()
  }
  return(p)
}
map(pops_with_gnomad_data, function(x) {
  plot_off_diagonals(x, n_points = 2000, true_first = F, color_guide_title = 'gnomAD\nPASS', save_plot=T)
  plot_off_diagonals(x, 'genotyped_variant', 1000, color_guide_title = 'Genotyped', save_plot=T)
  plot_off_diagonals(x, 'in_ld_with_hwe_problematic_variant', 2500, color_guide_title = 'In LD with\nbad HWE variant', save_plot=T)
})

sfigure_x1 = function(pop_to_use = 'AFR') {
  left = plot_ukb_gnomad_2dhist(pop_to_use, 'identity', breaks=c(1, 2.5e5, 5e+5))
  right = plot_off_diagonals(pop_to_use, n_points = 2000, true_first = F, color_guide_title = 'gnomAD\nPASS')
  
  output_type(output_format, paste0('plots/s_figure_x_ukb_gnomad_', pop_to_use, '.', output_format), height=3, width=6.5)
  print(ggpubr::ggarrange(left, right, nrow = 1, ncol = 2, legend='bottom', align='h', labels=c('a', 'b')))
  dev.off()
}

plot_ukb_gnomad_ratio_density = function(pop_to_use, color_by='pass_gnomad_genomes', color_guide_title='', save_plot=F) {
  if (save_plot) output_type(output_format, paste0('plots/', 'ukb_gnomad_ratio_density_', pop_to_use, '_by_', color_by, '.', output_format), height=4, width=5)
  
  p = data_scatter %>%
    filter(pop == pop_to_use)  %>%
    mutate(fold_change=ukb_freq / gnomad_freq) %>%
    ggplot + aes(x = fold_change) + 
    aes_string(group = color_by, fill = color_by) +
    geom_density(alpha = 0.5) + 
    scale_x_log10(name=ukb_gnomad_ratio_name)
  
  if (color_guide_title != '') {
    p = p + scale_fill_few(name=color_guide_title, palette = 'Dark')
  }
  
  if (save_plot) {
    print(p)
    dev.off()
  }
  return(p)
}
map(pops_with_gnomad_data, function(x) plot_ukb_gnomad_ratio_density(x, color_guide_title = 'gnomAD\nPASS', save_plot=T))

data_regional = load_ukb_file(paste0('n_variants_by_region_gnomad_status_full.txt.bgz'), 
                              'gnomad_comparison/', force_cols=cols(chrom=col_character()))


plot_snp_density = function(data_regional, pop_to_use = 'EUR', save_plot = F) {
  library(ggbio)
  library(GenomicRanges)
  if (save_plot) png(paste0('regional_missing_from_gnomad_', pop_to_use, '.png'), height=6, width=6, res=300, units='in')

  data(hg19Ideogram, package = "biovizBase")
  # Collapse bins further and get end coordinate
  binned = data_regional %>%
    filter(pop == pop_to_use) %>%
    group_by(gnomad_pass, chrom, start = as.integer(kb / 10) * 10 * 1000,
             end = pmin((as.integer(kb / 10) + 1) * 10000 - 1, seqlengths(hg19Ideogram)[paste0("chr", chrom)])) %>%
    summarize(n_variants = sum(n_variants)) %>% ungroup
  # Creating GRanges objects
  density <- binned %>% 
    filter(is.na(gnomad_pass)) %>%
    with(GRanges(paste0("chr", chrom), IRanges(start, end), strand="+", n_variants))
  # Consistent chromosome lengths/which chromosomes to use
  seqlengths(density) <- seqlengths(hg19Ideogram)[names(seqlengths(density))]
  density <- keepSeqlevels(trim(density), paste0("chr", c(1:22, "X")))
  #autoplot(density, layout = "karyogram")
  
  p = autoplot(density, layout = "karyogram", aes(color=n_variants, fill=n_variants, size=n_variants)) + # cheating a bit since it increases size of box based on frequency
    scale_colour_gradient(low = low_color, high = high_color, space = "Lab", na.value = low_color, guide = "colourbar") +
    scale_fill_gradient(low = low_color, high = high_color, space = "Lab", na.value = low_color, guide = "colourbar") +
    guides(size=F)
  if (save_plot) {
    print(p)
    dev.off()
  }
  return(p)
}
plot_snp_density(data_regional, 'EUR', T)

load_all_ukb_gnomad_cdf_data = function() {
  map_df(pops_with_gnomad_data, function(x) load_ukb_file(paste0('fold_change_ukb_over_gnomad_af_', x, '.txt.bgz'), 'gnomad_comparison/') %>% 
           mutate(pop=x, fold_change = 10 ^ values, percentile = 1 - ranks / max(ranks))) %>%
    return
}
data_cdf = load_all_ukb_gnomad_cdf_data()

plot_ukb_gnomad_ratio_cdf = function(zoom = 0, save_plot=F) {
  if (save_plot) output_type(output_format, paste0('plots/ukb_gnomad_ratio_cdf', if_else(zoom > 0, paste0('_', zoom), ''), '.', output_format), height=4, width=5)
  
  p = data_cdf %>%
    ggplot + aes(x = fold_change, y = percentile, group = pop, color = pop) + 
    geom_line(lwd = 1) + pop_color_scale +
    scale_x_log10(name=ukb_gnomad_ratio_name) +
    geom_vline(xintercept = 2, linetype = 'dashed')
  
  if (zoom != 0) {
    p = p + coord_cartesian(xlim=c(1 / zoom, zoom))
  }
  
  if (save_plot) {
    print(p)
    dev.off()
  }
  return(p)
}
plot_ukb_gnomad_ratio_cdf(zoom = 3, save_plot = T)
plot_ukb_gnomad_ratio_cdf(save_plot = T)

load_all_ukb_gnomad_titv_data = function() {
  map_df(pops_with_gnomad_data, function(x) load_ukb_file(paste0('titv_by_fold_change_ukb_over_gnomad_af_', x, '.txt.bgz'),
                                         'gnomad_comparison/') %>% mutate(pop=x))  %>%
    mutate(fold_change_bin=10 ^ fold_change_bin, freq_bin=10^freq_bin / 10, titv = tis / tvs) %>%
    return
}
data_titv = load_all_ukb_gnomad_titv_data()

plot_ukb_gnomad_titv = function(pop_to_use, min_freq = 0, save_plot=F) {
  fname = paste0('plots/ukb_gnomad_ratio_titv_', pop_to_use)
  if (min_freq > 0) fname = paste0(fname, '_min_freq_', min_freq)
  fname = paste0(fname, '.', output_format)
  
  if (pop_to_use == 'all') {
    if (save_plot) output_type(output_format, fname, height=4, width=6.5)
    p = data_titv %>% 
      filter(freq_bin > min_freq) %>%
      group_by(pop, fold_change_bin) %>%
      summarize(tis=sum(tis), tvs=sum(tvs), titv=tis/tvs, n_variants=tis+tvs) %>% ungroup %>%
      ggplot + aes(x = fold_change_bin, y = titv, size = n_variants, color=pop, fill=pop) + 
      annotate('rect', xmin=0.5, xmax=2, ymin=0, ymax=Inf, alpha=0.4, fill='lightblue') +
      geom_point() + #geom_line(lwd=0.5) +
      pop_color_scale + pop_fill_scale +
      scale_x_log10(name=ukb_gnomad_ratio_name) +
      geom_vline(xintercept = 1, linetype='dashed') +
      scale_size_continuous(name='Number of\nvariants') +
      ylab('Ti/Tv ratio')
  } else {
    if (save_plot) output_type(output_format, fname, height=4, width=5)
    p = data_titv %>%
      filter(pop == pop_to_use & freq_bin > min_freq) %>%
      group_by(fold_change_bin) %>%
      summarize(tis=sum(tis), tvs=sum(tvs), titv=tis/tvs) %>%
      ggplot + aes(x = fold_change_bin, y = titv, size = tis + tvs) + 
      geom_point(color = get_color(pop_to_use)) + 
      scale_x_log10(name=ukb_gnomad_ratio_name) +
      geom_vline(xintercept = 1, linetype='dashed') +
      scale_size_continuous(name='Number of\nvariants') +
      ylab('Ti/Tv ratio')
  }
  
  if (save_plot) {
    print(p)
    dev.off()
  }
  return(p)
}
map(pops_with_gnomad_data, function(x) plot_ukb_gnomad_titv(x, save_plot=T))
map(pops_with_gnomad_data, function(x) plot_ukb_gnomad_titv(x, min_freq=0.01, save_plot=T))
plot_ukb_gnomad_titv('all', min_freq=0.01, save_plot=T)

sfigure_x2 = function() {
  plots = map(pops_with_gnomad_data, function(x) plot_ukb_gnomad_titv(x, min_freq=0.01) +
                annotate('text', label = x, x = Inf, y = Inf, hjust = 1, vjust = 1, color=ukb_pop_colors[[x]], size=3)+
                annotate('rect', xmin=0.5, xmax=2, ymin=0, ymax=Inf, alpha=0.4, fill='lightblue') )
  
  output_type(output_format, paste0('plots/s_figure_x_ukb_gnomad_titv.', output_format), height=6, width=6.5)
  print(ggpubr::ggarrange(plotlist = plots, nrow = 2, ncol = 2, legend='bottom', align='h', labels=c('a', 'b', 'c', 'd'),
                          common.legend = T))
  dev.off()
}

plot_ukb_gnomad_ratio_cdf_by_freq = function(zoom = 0, min_freq = 0, save_plot=F) {
  fname = paste0('plots/ukb_gnomad_ratio_cdf_min_freq_', min_freq, if_else(zoom > 0, paste0('_', zoom), ''), '.', output_format)
  if (save_plot) output_type(output_format, fname, height=4, width=5)
  
  plot_data = data_titv %>%
    filter(freq_bin >= min_freq) %>%
    group_by(pop, fold_change_bin = if_else(is.na(fold_change_bin), 2000, fold_change_bin)) %>%
    summarize(tis=sum(tis), tvs=sum(tvs), titv=tis/tvs, n_variants=tis+tvs) %>% ungroup

  p = plot_data %>% group_by(pop) %>%
    mutate(n_variants=cumsum(n_variants),
           prop = 1 - n_variants / max(n_variants)) %>%
    ggplot + aes(x = fold_change_bin, y = prop, color = pop) + pop_color_scale +
    geom_vline(xintercept = 2, linetype='dashed') +
    geom_line(lwd=1) + scale_x_log10(name=ukb_gnomad_ratio_name) + ylab('Percentile')
  if (zoom != 0) {
    p = p + coord_cartesian(xlim=c(1 / zoom, zoom))
  }
  
  if (save_plot) {
    print(p)
    dev.off()
  }
  return(p)
}
plot_ukb_gnomad_ratio_cdf_by_freq(save_plot = T)
plot_ukb_gnomad_ratio_cdf_by_freq(zoom = 3, save_plot = T)
plot_ukb_gnomad_ratio_cdf_by_freq(min_freq = 0.01, save_plot = T)
plot_ukb_gnomad_ratio_cdf_by_freq(min_freq = 0.01, zoom = 3, save_plot = T)

plot_ukb_gnomad_ratio_summary = function(zoom = 0, min_freq = 0, plot_prop_instead_of_n = F, show_only_pop = '', save_plot=F) {
  fname = paste0('plots/ukb_gnomad_ratio_', case_when(show_only_pop != '' ~ paste0('pop_', show_only_pop),
                                                      plot_prop_instead_of_n ~ 'summary_prop', 
                                                      TRUE ~ 'summary'), '_min_freq_', min_freq, 
                 if_else(zoom > 0, paste0('_', zoom), ''), '.', output_format)
  if (save_plot) output_type(output_format, fname, height=if_else(show_only_pop != '', 5, 4), 
                     width=if_else(show_only_pop != '', 7, 5))
  
  plot_data = data_titv %>%
    filter(freq_bin >= min_freq) %>%
    mutate(fold_change_bin = if_else(is.na(fold_change_bin), 2000, fold_change_bin),
           gnomad_status = case_when(fold_change_bin == 2000 ~ 'Not in gnomAD',
                                     fold_change_bin == 1000 ~ 'In gnomAD, 0 in population',
                                     fold_change_bin >= 2 ~ 'In gnomAD, >2X frequency',
                                     TRUE ~ 'Well-calibrated'),
           n_variants = tis + tvs)
  
  if (show_only_pop == '') {
    p = plot_data %>%
      group_by(pop, gnomad_status) %>%
      summarize(n_variants = sum(n_variants)) %>%
      mutate(prop_variants = n_variants / sum(n_variants)) %>%
      ggplot + aes(x = pop, fill = gnomad_status) + 
      aes_string(y = if_else(plot_prop_instead_of_n, 'prop_variants', 'n_variants')) +
      geom_bar(stat='identity') + 
      scale_fill_few(name=NULL) + xlab('Population') +
      ylab(if_else(plot_prop_instead_of_n, 'Proportion of variants', 'Number of variants'))
  } else {
    p = plot_data %>%
      filter(pop == show_only_pop) %>%
      group_by(freq_bin, gnomad_status) %>%
      summarize(n_variants=sum(n_variants)) %>% ungroup %>%
      ggplot + aes(x = freq_bin, fill = gnomad_status, y = n_variants) + 
      geom_bar(stat='identity') + 
      scale_fill_few(name=NULL) + scale_x_log10(label=pretty_axis_format_binned_corrected, name='Frequency bin') +
      ylab('Number of variants') + theme(legend.position = 'bottom')
  }

  if (save_plot) {
    print(p)
    dev.off()
  }
  return(p)
}
plot_ukb_gnomad_ratio_summary(save_plot = T)
plot_ukb_gnomad_ratio_summary(min_freq = 0.01, save_plot = T)
plot_ukb_gnomad_ratio_summary(plot_prop_instead_of_n = T, save_plot = T)
map(pops_with_gnomad_data, function(x) plot_ukb_gnomad_ratio_summary(show_only_pop = x, save_plot = T))

# plot_not_in_gnomad = function(pop_to_use, plot_titv = F, save_plot=F) {
#   if (save_plot) output_type(output_format, paste0('plots/ukb_not_in_gnomad_', if_else(plot_titv, 'titv_', 'n_'), pop_to_use, '.', output_format), height=5, width=7)
#   
#   plot_data = data_titv %>%
#     filter(pop == pop_to_use) %>%
#     group_by(freq_bin, in_gnomad = !is.na(fold_change_bin)) %>%
#     summarize(tis=sum(tis), tvs=sum(tvs), titv=tis/tvs) %>% ungroup
#   print(paste('Pop:', pop_to_use))
#   print(plot_data %>% group_by(in_gnomad) %>% summarize(total=sum(tis) + sum(tvs)))
#   if (plot_titv) {
#     p = plot_data %>%
#       ggplot + aes(x = freq_bin, y = titv, size = tis + tvs, shape = in_gnomad) + 
#       geom_point(color = get_color(pop_to_use)) +
#       scale_size_continuous(name='Number of\nvariants') +
#       ylab('Ti/Tv ratio') + scale_shape_discrete(name='In gnomAD')
#   } else {
#     p = plot_data %>%
#       ggplot + aes(x = freq_bin, y = tis + tvs, fill = in_gnomad) + 
#       geom_bar(stat='identity', position='dodge') +
#       ylab('Number of variants') + scale_fill_few(name='In gnomAD', palette = 'Dark')
#   }
#   p = p + scale_x_log10(label=pretty_axis_format_binned, name='Frequency bin')
#   
#   if (save_plot) {
#     print(p)
#     dev.off()
#   }
#   return(p)
# }
# map(pops_with_gnomad_data, function(x) plot_not_in_gnomad(x, save_plot=T))
# map(pops_with_gnomad_data, function(x) plot_not_in_gnomad(x, plot_titv=T, save_plot=T))

load_all_ukb_gnomad_volcano_data = function() {
  map_df(pops_with_gnomad_data, function(x) load_ukb_file(paste0('or_vs_p_hist2d_', x, '.txt.bgz'),
                                                          'gnomad_comparison/') %>% mutate(pop=x)) %>%
    mutate(odds_ratio = 10 ^ odds_ratio) %>%
    return
}
data_volcano = load_all_ukb_gnomad_volcano_data()

compute_cumulative_volcano = function(data_volcano) {
  total_variants = data_volcano %$% sum(n_variants)
  data_volcano %>%
    complete(p_value, odds_ratio, fill = list(n_variants = 0)) %>%
    group_by(p_value) %>%
    arrange(desc(odds_ratio)) %>%
    mutate(cumulative_variants = cumsum(n_variants)) %>%
    group_by(odds_ratio) %>%
    arrange(desc(p_value)) %>%
    mutate(cumulative_variants = cumsum(cumulative_variants)) %>%
    ungroup %>% mutate(prop_variants_filtered = cumulative_variants / total_variants) %>%
    return
}
cumulative_volcano = map_df(pops_with_gnomad_data, 
                            function(x) compute_cumulative_volcano(data_volcano %>% filter(pop == x)) %>% mutate(pop = x))

plot_gnomad_enrichment_volcano = function(pop_to_use, cumulative = F, save_plot=F) {
  if (save_plot) output_type(output_format, paste0('plots/', 'or_vs_p_hist2d_', pop_to_use, if_else(cumulative, '_cumulative', ''), '.', output_format), height=4, width=5)
  
  if (!cumulative) {
    p = data_volcano %>%
      filter(pop == pop_to_use) %>%
      ggplot + aes(fill = n_variants)
  } else {
    p = cumulative_volcano %>% 
      filter(pop == pop_to_use) %>%
      ggplot + aes(fill = cumulative_variants, text = prop_variants_filtered)
  }
  p = p + scale_x_log10(name=ukb_gnomad_ratio_name) + ylab('-log10(p)') + 
    geom_tile() + aes(x = odds_ratio, y = p_value)
    scale_fill_continuous(low = low_color, high = high_color, name='Number of\nvariants', trans='log10')
  
  if (save_plot) {
    print(p)
    dev.off()
  }
  return(p)
}
plot_gnomad_enrichment_volcano('EUR', save_plot=T)
ggplotly(plot_gnomad_enrichment_volcano('EUR', T, save_plot=T))

plot_n_variants_by_gnomad_enrichment_fixed_p = function(pop_to_use, neglog_p_cutoff = 2, save_plot=F) {
  if (save_plot) output_type(output_format, paste0('plots/ukb_gnomad_enrichment_fixed_p.', output_format), height=4, width=5)
  
  p = cumulative_volcano %>%
    filter(p_value == neglog_p_cutoff) %>%
    ggplot + aes(x = odds_ratio, y = prop_variants_filtered, color = pop, group = pop) +
    geom_line(lwd=1.5) + pop_color_scale +
    scale_x_log10(name=ukb_gnomad_ratio_name) +
    geom_vline(xintercept = 2, linetype = 'dashed') +
    ylab('Proportion variants filtered')
  
  if (save_plot) {
    print(p)
    dev.off()
  }
  return(p)
}
plot_n_variants_by_gnomad_enrichment_fixed_p(save_plot=T)

cross_population_data = load_ukb_file(paste0('n_variants_by_failing_pops_full.txt.bgz'), 'gnomad_comparison/')

plot_cross_population = function(save_plot=F) {
  library(ComplexUpset)
  library(UpSetR)
  
  if (save_plot) output_type(output_format, paste0('plots/ukb_gnomad_failing_cross_population.', output_format), height=4, width=6.5, onefile = T)
  plot_data = cross_population_data %>%
    filter(!is.na(failing_pops)) %>%
    mutate(AFR=grepl('AFR', failing_pops),
           AMR=grepl('AMR', failing_pops),
           EAS=grepl('EAS', failing_pops),
           EUR=grepl('EUR', failing_pops))
  
  temp = cross_population_data %>%
    filter(!is.na(failing_pops)) %>%
    transmute(failing_pops=gsub(',', '&', failing_pops), n)
  
  data = temp$n
  names(data) = temp$failing_pops
  
  # TODO: appears to be setting colors incorrectly
  # p = ComplexUpset::upset(fromExpression(data), c('AFR', 'AMR', 'EAS', 'EUR'),
  #                         queries=list(
  #                           upset_query(set='AFR', fill=color_afr),
  #                           upset_query(set='AMR', fill=color_amr),
  #                           upset_query(set='EAS', fill=color_eas),
  #                           upset_query(set='EUR', fill=color_eur)
  #                         ))
  
  p = upset(fromExpression(data))
  if (save_plot) {
    print(p)
    dev.off()
  }
  return(p)
}

plot_cross_population(save_plot=T)

sfigure_x3 = function(pop = 'EUR') {
  cross_population_data %>% filter(!in_gnomad_genomes)
  cross_population_data %>% filter(in_gnomad_genomes & !is.na(failing_pops)) %>% summarize(total_variants=sum(n))
  
  left = plot_ukb_gnomad_ratio_summary(show_only_pop = pop)
  right = plot_ukb_gnomad_ratio_summary()
  
  output_type(output_format, paste0('plots/s_figure_x_missing_from_gnomad.', output_format), height=4.5, width=9.75)
  print(ggpubr::ggarrange(left, right, nrow = 1, ncol = 2, legend='bottom', align='h', labels=c('a', 'b'),
                          common.legend = T))
  dev.off()
}

cross_population_data_by_pop = load_ukb_file(paste0('n_variants_by_failing_pops_by_pop_full.txt.bgz'), 'gnomad_comparison/')
cross_population_data_by_pop %>%
  count(pop, high_quality=is.na(failing_pops) & in_gnomad_genomes, wt = n_variants) %>%
  ggplot + aes(x = pop, y = n, fill = high_quality) + 
  geom_bar(stat='identity') +
  scale_fill_few('Dark', name='High quality') + 
  xlab('Population') + ylab('Number of variants')


ukb_ans = c("AFR" = 9220 * 2, "EAS" = 2900 * 2, "AMR" = 1148 * 2, "EUR" = 458937 * 2)
gnomad_ans = c("AFR" = 4359 * 2, "EAS" = 780 * 2, "AMR" = 424 * 2, "EUR" = 4299 * 2)

simulate_ukb_gnomad_ratio = function(ukb_an = 917874, gnomad_an = 6142, pop = 'EUR_8-75857876', y_axis_ukb = F, 
                                     width = 0.02, cutoff = 1e-6, save_plot = F) {
  if (save_plot) output_type(output_format, paste0('plots/simulated_ukb_gnomad_', pop, '.', output_format), height=4, width=5)
  freqs = 10 ^ seq(-4, 0, width)

  p_values_by_freq = map_df(freqs, function(u) map_df(freqs, function(g) {
    data.frame(ukb_freq=u, gnomad_freq=g, 
                 p_value=chisq.test(matrix(c(ukb_an * u, ukb_an * (1-u), gnomad_an * g, gnomad_an * (1-g)), ncol = 2))$p.value)
    }))
  
  plot_data = p_values_by_freq %>%
    mutate(p = if_else(p_value < cutoff, -log10(cutoff), -log10(p_value)))
  
  if (y_axis_ukb) {
    p = p_values_by_freq %>%
      ggplot + aes(x = gnomad_freq, y = ukb_freq, fill = p) +
      geom_tile() + 
      scale_x_log10(name='gnomAD frequency') + scale_y_log10(name='UKB frequency')
  } else {
    p = plot_data %>%
      ggplot + aes(x = gnomad_freq, y = gnomad_freq / ukb_freq, fill = p) +
      geom_tile(height=width, width=width) + 
      scale_y_log10(name=ukb_gnomad_ratio_name, breaks=c(1, 2, 5, 10, 20, 50)) + 
      scale_x_log10(name=paste0('gnomAD frequency (', pop, ')'), label=pretty_axis_format) +
      geom_hline(yintercept = 2, linetype='dashed') +
      scale_fill_continuous(low = low_color, high = high_color, name='-log10(p)') +
      coord_cartesian(ylim=c(1, 20))
  }
  
  if (save_plot) {
    print(p)
    dev.off()
  }
  return(p)
}
simulate_ukb_gnomad_ratio()
map(pops_with_gnomad_data, function(x) simulate_ukb_gnomad_ratio(ukb_an = ukb_ans[[x]], gnomad_an = gnomad_ans[[x]], pop = x, save_plot=T))




