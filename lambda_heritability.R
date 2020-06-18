source('~/ukbb_pan_ancestry/constants.R')

output_format = 'png'

final_lambda_data = load_ukb_file('lambdas_full.txt.bgz', 'lambda/',
                                  force_cols=cols(phenocode=col_character(), coding=col_character())) %>%
  mutate(lambda_gc = pheno_data.lambda_gc, pop = pheno_data.pop,
         n_cases=pheno_data.n_cases, n_controls=pheno_data.n_controls,
         lambda1000 = if_else(is.na(n_controls), 1 + (lambda_gc - 1) * 2000 / n_cases,
                              1 + (lambda_gc - 1) * (1 / n_cases + 1 / n_controls) * 500))

lambda_pre_qc = load_all_lambdas_data('ac') %>% filter(ac == 0) %>% 
  mutate(lambda_gc = lambda_gc_by_case_ac) %>% 
  union_all(load_all_lambdas_data('af') %>% filter(af == 0) %>% mutate(lambda_gc = lambda_gc_by_af)) %>% 
  mutate(lambda1000 = if_else(is.na(n_controls), 1 + (lambda_gc - 1) * 2000 / n_cases,
                              1 + (lambda_gc - 1) * (1 / n_cases + 1 / n_controls) * 500))


# Raw lambdas before QC
lambda_pre_qc %>%
  filter(lambda_gc < 10) %>%
  group_by(pop) %>%
  summarize(median_lambda = median(lambda_gc, na.rm=T))

# Post QC
final_lambda_data %>% 
  filter(lambda_gc < 10) %>%
  group_by(pop) %>%
  summarize(median_lambda=median(lambda_gc))

final_lambda_data %>%
  group_by(pop) %>%
  summarize(outliers=sum(lambda_gc > 2))

final_lambda_data %>% 
  filter(pop == 'EUR') %>%
  filter(lambda_gc < 2) %>%
  ggplot + aes(x = lambda_gc) + geom_histogram() +
  xlab('Lambda GC') + ylab('Number of phenotypes')

plot_final_lambda_gc_distr = function(save_plot = F) {
  if (save_plot) output_type(output_format, paste0('plots/lambda_gc_distribution.', output_format), height=4, width=5)
  
  final_lambda_data %>%
    filter(lambda_gc > 1.2) %>%
    count(pop)
  

  final_lambda_data %>%
    filter(lambda_gc > 1.2) %>%
    count(pop, trait_type)
  
  low_cutoff = 0.5
  high_cutoff = 2
  filter_lambdas = function(x) {
    return(case_when(x > high_cutoff ~ high_cutoff, 
              x < low_cutoff ~ low_cutoff,
              TRUE ~ x))
  }
  
  final_lambda_data %>% 
    count(lambda_cutoffs = lambda_gc >= 1.2 | lambda_gc <= 0.8,
          lambda1000_cutoffs = lambda1000 >= 1.2 | lambda1000 <= 0.8)
  
  final_lambda_data %>% 
    count(lambda_cutoffs = lambda_gc >= 2 | lambda_gc <= 0.5,
          lambda1000_cutoffs = lambda1000 >= 2 | lambda1000 <= 0.5)
  
  p = final_lambda_data %>% 
    mutate(lambda_gc = filter_lambdas(lambda_gc),
           lambda1000 = filter_lambdas(lambda1000)) %>%
    ggplot + aes(x = lambda1000, group = pop, fill = pop) + geom_histogram() +
    pop_fill_scale + scale_x_log10() +
    xlab('Lambda GC') + ylab('Number of phenotypes')
  p
  if (save_plot) {
    print(p)
    dev.off()
  }
  return(p)
}
plot_final_lambda_gc_distr(T)

filter_phenos = function(df_to_filter, full_df, input_pop='all') {
  if (input_pop == 'all') {
    x = full_df %>% count(trait_type, phenocode, pheno_sex, coding, modifier) %>% filter(n > 5)
  } else{
    x = full_df %>% filter(pop == input_pop)
  }
  x %<>% select(trait_type, phenocode, pheno_sex, coding, modifier)
  return(semi_join(df_to_filter, x))
}

# final_lambda_data %>%
#   filter(pop == 'EUR') %>%
#   filter_phenos(final_lambda_data, 'AFR') %>%
#   filter(lambda_gc < 2) %>%
#   ggplot + aes(x = lambda_gc) + 
#   geom_histogram() + 
#   xlab('Lambda GC') + ylab('Number of phenotypes')

plot_lambda_gc_by_case_count = function(pre_qc = F, plot_lambda_1000 = F, save_plot = F) {
  if (save_plot) output_type(output_format, paste0('plots/lambda_', if_else(plot_lambda_1000, '1000', 'gc'),
                                                   '_by_case_count', 
                                                   if_else(pre_qc, '_pre_qc', ''), '.', output_format), height=4, width=5)
  plot_data = if (pre_qc) lambda_pre_qc else final_lambda_data
  
  p = plot_data %>% 
    filter(pop == 'EUR') %>%
    filter(lambda_gc < 2) %>%
    ggplot + aes(x = n_cases, color = trait_type, description = description, coding_description = coding_description) + 
    aes_string(y = if_else(plot_lambda_1000, 'lambda1000', 'lambda_gc')) +
    geom_point() + scale_x_log10(labels=comma, name='Number of cases') + 
    scale_color_few('Dark', name='Trait type') +
    ylab(if_else(plot_lambda_1000, 'Lambda 1000', 'Lambda GC'))
  
  if (save_plot) {
    print(p)
    dev.off()
  }
  return(p)
}
plot_lambda_gc_by_case_count(pre_qc = T, save_plot = T)
plot_lambda_gc_by_case_count(plot_lambda_1000 = T, save_plot = T)
p = plot_lambda_gc_by_case_count(pre_qc = F, save_plot = T)
ggplotly(p)

lambda_with_qc = final_lambda_data %>%
  mutate(bad_lambda=(lambda_gc > 1.2) | (lambda_gc < 0.8),
         too_many_sig=pheno_data.n_sig_variants > 0.5 * pheno_data.n_variants)

# All "too many sig" are captured in "bad lambda" (no matter the cutoff, even at 0.1 and 10)
lambda_with_qc %>%
  count(pop, bad_lambda, too_many_sig)
  
plot_n_phenos_by_passing_qc = function(save_plot = F) {
  if (save_plot) output_type(output_format, paste0('plots/n_phenos_by_passing_qc.', output_format), height=4, width=5)

  p = lambda_with_qc %>%
    group_by(trait_type, phenocode, pheno_sex, coding, modifier) %>%
    summarize(n_good=sum(!bad_lambda),
              n_total=n()) %>% ungroup %>%
    count(n_total, n_good) %>%
    ggplot + aes(x=n_total, y = n_good, fill = n, label=n) + geom_tile() +
    scale_x_continuous(breaks=0:6, name='Number of populations where phenotype was run') +
    scale_y_continuous(breaks=0:6, name='Number of populations with calibrated lambda') +
    geom_text(color='white') +
    scale_fill_continuous(low = low_color, high = high_color, name='Number of\nphenotypes', trans='log10')
  if (save_plot) {
    print(p)
    dev.off()
  }
  return(p)  
}
plot_n_phenos_by_passing_qc(T)

lambda_with_qc %>%
  filter(!bad_lambda) %>%
  count(pheno_data.n_sig_variants == 0)

lambda_with_qc %>% 
  filter(!bad_lambda) %>%
  mutate(rank = min_rank(pheno_data.n_sig_variants)) %>%
  ggplot + aes(x = pheno_data.n_sig_variants - 1, y = rank) + 
  geom_line() + geom_point() + scale_x_log10(name='Number of significant hits (p < 5e-8)') +
  ylab('Number of phenotypes') + coord_cartesian(ylim=c(9000, 15300))

lambda_with_qc %>% 
  filter(!bad_lambda) %>%
  group_by(pop) %>%
  mutate(rank = min_rank(pheno_data.n_sig_variants)) %>% ungroup %>%
  ggplot + aes(x = pheno_data.n_sig_variants - 1, y = rank, color = pop) + 
  pop_color_scale + ylab('Number of phenotypes') +
  geom_line(lwd=1) + scale_x_log10(name='Number of significant hits (p < 5e-8)')

lambda_with_qc %>% 
  filter(!bad_lambda) %>%
  group_by(pop) %>%
  mutate(rn = n() - min_rank(pheno_data.n_sig_variants)) %>% ungroup %>%
  # filter(pheno_data.n_sig_variants > 0) %>%
  ggplot + aes(x = pheno_data.n_sig_variants, y = rn, color = pop) + 
  pop_color_scale +
  geom_line(lwd=1) + scale_y_log10(name='Number of phenotypes') +
  scale_x_log10(name='Number of significant hits (p < 5e-8) >= x')

lambda_with_qc %>% 
  filter(!bad_lambda) %>%
  group_by(pop) %>%
  mutate(rn = (n() - min_rank(pheno_data.n_sig_variants)) / n()) %>% ungroup %>%
  filter(pheno_data.n_sig_variants > 0) %>%
  ggplot + aes(x = pheno_data.n_sig_variants, y = rn, color = pop) + 
  pop_color_scale +
  geom_line(lwd=1) + scale_y_continuous(name='Proportion of phenotypes', label=function(x) percent(x, accuracy=10)) +
  scale_x_log10(name='Number of significant hits (p < 5e-8) >= x')

final_lambda_data %>% 
  ggplot + aes(x = n_cases, y = pheno_data.n_variants, color = pop) +
  geom_point() + scale_x_log10(name='Number of cases') +
  pop_color_scale + ylab('Number of variants passing QC')

final_lambda_data %>% 
  filter(pop == 'EUR') %>%
  ggplot + aes(x = n_cases, y = pheno_data.n_variants, color = trait_type, description = description) +
  geom_point() + scale_x_log10(name='Number of cases') +
  ylab('Number of variants passing QC') -> p
ggplotly(p)

function() {
  final_lambda_data %>%
    filter_phenos(final_lambda_data) %>%
    group_by(pop) %>%
    summarize(mean_heritability = mean(pheno_data.heritability))
  
  final_lambda_data %>%
    #   filter_phenos(final_lambda_data) %>%
    select(trait_type, phenocode, pheno_sex, coding, modifier, pop, pheno_data.heritability) %>%
    pivot_wider(names_from=pop, values_from=pheno_data.heritability) %>%
    ggplot + aes(x = AFR, y = EUR, color = trait_type) + 
    scale_color_few('Dark', name='Trait type') +
    geom_point() + geom_abline(slope = 1, intercept = 0, linetype = 'dashed')
  
  final_lambda_data %>%
    select(trait_type, phenocode, pheno_sex, coding, modifier, pop, pheno_data.heritability) %>%
    pivot_wider(names_from=pop, values_from=pheno_data.heritability) %>%
    ggplot + aes(x = AFR, y = CSA, color = trait_type) + 
    scale_color_few('Dark', name='Trait type') +
    geom_point() + geom_abline(slope = 1, intercept = 0, linetype = 'dashed')
  
  final_lambda_data %>%
    select(trait_type, phenocode, pheno_sex, coding, modifier, pop, pheno_data.heritability) %>%
    pivot_wider(names_from=pop, values_from=pheno_data.heritability) %>% 
    select(AFR:MID) %>%
    cor(use='complete.obs')
  
  final_lambda_data %>%
    group_by(trait_type, phenocode, pheno_sex, coding, modifier) %>%
    summarize(mean_heritability=mean(pheno_data.heritability)) %>% 
    mutate(quant_trait = trait_type %in% c('biomarkers', 'continuous')) %>%
    ggplot + aes(x = mean_heritability) +
    aes(fill = trait_type) + scale_fill_few(name='Trait type') +
    # aes(fill = quant_trait) +
    geom_density(alpha=0.2) 
  
  final_lambda_data %>%
    filter(pop == 'EUR' & trait_type %in% c('biomarkers', 'continuous')) %>%
    ggplot + aes(x = pheno_data.heritability) +
    aes(fill = trait_type) + scale_fill_few(name='Trait type') +
    # aes(fill = quant_trait) +
    geom_density(alpha=0.2) 
}
