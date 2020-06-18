source('~/ukbb_pan_ancestry/constants.R')

output_format = 'png'

load_all_lambdas_data = function(by='ac') {
  map_df(pops, function(x) load_ukb_file(paste0('lambdas_by_', by, '_', x, '.txt.bgz'), 'lambda/',
                                                          force_cols = cols(coding=col_character(),
                                                                            phenocode=col_character(), 
                                                                            modifier=col_character())
                                                          ) %>% mutate(pop=x)) %>%
    filter_at(vars(paste0('lambda_gc_by_', if_else(by == 'ac', 'case_ac', 'af'))), all_vars(!is.na(.))) %>%
    return
}

lambda_data_ac = load_all_lambdas_data('ac')


plot_lambdas_by_case_ac = function(plot_n_variants = F, plot_prop_variants = F, plot_loess = F, 
                                   plot_pointrange = F, save_plot = F) {
  metric = case_when(plot_n_variants ~ 'n', 
                     plot_prop_variants ~ 'prop',
                     TRUE ~ 'lambda')
  if (save_plot) output_type(output_format, paste0('plots/', metric, '_by_case_ac',
                            if_else(plot_loess, '_loess', ''), '.', output_format), height=3, width=6.5)
  
  plot_data = lambda_data_ac %>%
    group_by(pop, ac) %>%
    summarize(n_variants = mean(n_variants_by_case_ac, na.rm=T),
              median_lambda = median(lambda_gc_by_case_ac, na.rm=T),
              mad_lambda = mad(lambda_gc_by_case_ac, na.rm=T),
              n_phenos=n(),
              sem = 1.96 * mad_lambda * 1.4826 / sqrt(n_phenos)) %>% ungroup
  plot_data = plot_data %>%
    left_join(plot_data %>% filter(ac == 0) %>% select(pop, total_variants=n_variants)) %>%
    mutate(proportion_variants = n_variants / total_variants)
  
  if (plot_n_variants | plot_prop_variants) {
    p = plot_data %>%
      ggplot + aes(x = ac, color = pop) + aes_string(y = if_else(plot_prop_variants, 'proportion_variants', 'n_variants')) +
      geom_point() +
      pop_color_scale +
      scale_x_log10(name='Allele Count >= x') + 
      ylab(if_else(plot_prop_variants, 'Proportion of variants', 'Number of variants'))
  } else {
    if (plot_loess) {
      p = lambda_data_ac %>%
        filter(lambda_gc_by_case_ac > 0.8 & lambda_gc_by_case_ac < 1.2) %>%
        ggplot + aes(x = ac + 0.5, y = lambda_gc_by_case_ac, group = pop, color = pop, fill = pop) + 
        geom_smooth() + pop_fill_scale + ylab('Lambda')
    } else {
      p = plot_data %>%
        ggplot + aes(x = ac, y = median_lambda, color = pop,
                     ymin = median_lambda - sem, ymax = median_lambda + sem)
      if (plot_pointrange) {
        pj = position_jitter(width=0.1, height=0)
        p = p + geom_pointrange(position = pj) + geom_line(position = pj)
      } else {
        p = p + geom_point() + geom_line()
      }
      p = p +
        coord_cartesian(ylim=c(0.97, 1.03)) +
        ylab('Median lambda')
    }
    p = p + pop_color_scale +
      geom_hline(yintercept = 1, linetype='dashed') +
      scale_x_log10(name='Allele Count in cases', breaks=c(1,2,5,10,20,50,100)) 
  }
  if (save_plot) {
    print(p)
    dev.off()
  }
  return(p)
}
plot_lambdas_by_case_ac(save_plot=T)
plot_lambdas_by_case_ac(plot_loess=T, save_plot=T)
plot_lambdas_by_case_ac(plot_n_variants=T, save_plot=T)
plot_lambdas_by_case_ac(plot_prop_variants=T, save_plot=T)

# Sanity checking number of variants per pop total
lambda_data_ac %>%
  filter(ac == 0) %>%
  group_by(pop) %>%
  summarize(max_variants_per_pop=max(n_variants_by_case_ac),
            min_variants_per_pop=min(n_variants_by_case_ac))

lambda_data_ac %>%
  filter(ac == 0 & n_sig_by_case_ac > 0) %>%
  ggplot + aes(x = n_sig_by_case_ac, group = pop, fill = pop) +
  geom_histogram(position='dodge') + pop_fill_scale +
  scale_x_log10(name='Number of significant hits per phenotype') + 
  scale_y_continuous(name='Number of phenotypes')

lambda_data_ac %>%
  # filter(lambda_gc_by_case_ac > 0.8 & lambda_gc_by_case_ac < 1.2) %>%
  filter(ac == 0 & n_sig_by_case_ac > 10) %>%
  ggplot + aes(y = n_sig_by_case_ac, x = n_variants_by_case_ac, color = pop, 
               description = description, coding_description = coding_description) +
  geom_point() + pop_color_scale +
  scale_y_log10(name='Number of significant hits per phenotype') + 
  scale_x_log10(name='Number of variants assessed') -> p

chart_link = api_create(p, filename = "n_significant_by_pheno")

p
ggplotly(p)

lambda_data_af = load_all_lambdas_data('af')

lambda_data_af %>%
  filter(af == 0) %>%
  ggplot + aes(x = n_variants_by_af, fill=pop) + geom_histogram()

plot_lambdas_by_af = function(by_pop = F, save_plot = F) {
  if (save_plot) output_type(output_format, paste0('plots/lambda_by_af', if_else(by_pop, '_by_pop', ''), '.', output_format), height=4, width=5)
  
  if (by_pop) {
    p = lambda_data_af %>%
      group_by(pop, af) %>%
      summarize(median_lambda = median(lambda_gc_by_af, na.rm=T),
                mad_lambda = mad(lambda_gc_by_af, na.rm=T)) %>%
      ggplot + aes(x = af, y = median_lambda, color = pop) + 
      geom_point() + geom_line() + coord_cartesian(ylim=c(0.95, 1.05))
  } else {
    p = lambda_data_af %>%
      group_by(af) %>%
      summarize(median_lambda = median(lambda_gc_by_af, na.rm=T),
                mad_lambda = mad(lambda_gc_by_af, na.rm=T)) %>%
      ggplot + aes(x = af, y = median_lambda, ymin = median_lambda - mad_lambda, ymax = median_lambda + mad_lambda) + 
      geom_pointrange() 
  }
  
  p = p + pop_color_scale +
    geom_hline(yintercept = 1, linetype='dashed') +
    scale_x_log10(label=pretty_axis_format, name='Allele Frequency')
  if (save_plot) {
    print(p)
    dev.off()
  }
  return(p)
}
plot_lambdas_by_af(save_plot = T)
plot_lambdas_by_af(by_pop = T, save_plot = T)


