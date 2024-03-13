source('figure2.R')

ukb_31063_h2 <- read_tsv('data/ukb31063_h2_topline.02Oct2019.tsv') %>%
  mutate(ancestry = 'EUR') %>%
  transmute(phenotype, ancestry, intercept_r2 = intercept, 
            intercept_se_r2 = intercept_se,
            ratio_r2 = ratio, ratio_se_r2 = ratio_se, h2_liab_r2 = h2_liability,
            h2_liab_se_r2 = h2_liability_se, h2_z_r2 = h2_z)

# S Figure 21

round2_comp = function(sldsc=F) {
  h2_col = if_else(sldsc, 'estimates.sldsc_25bin.h2_liability', 'estimates.ldsc.h2_liability')
  ylab = if_else(sldsc, 'Pan-UKB EUR S-LDSC (liability)', 'Pan-UKB EUR LDSC (liability)')
  se_col = paste0(h2_col, '_se')
  p = h2_data %>% 
    filter(ancestry == 'EUR') %>%
    inner_join(ukb_31063_h2, by=c('code' = 'phenotype')) %>%
    select(h2_liab_r2, h2_liab_se_r2, h2y:={{h2_col}}, h2yse:={{se_col}}) %>%
    ggplot +
    aes(x = h2_liab_r2, y = h2y,
        alpha = sqrt(1/(h2yse^2) + 1/h2_liab_se_r2^2)) +
    geom_point() + geom_abline(slope=1, intercept=0, linetype='dashed') +
    guides(alpha=F) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=F) +
    xlab(expression(paste('Round 2 ', 'h'^'2',' (S-LDSC liability)'))) + ylab(ylab)
  return(p)
}

sig_h2_by_ancestry = function() {
  p = h2_data %>%
    filter(!is.na(estimates.sldsc_25bin.h2_z)) %>%
    mutate(ancestry=fct_relevel(ancestry, pops_by_sample_size),
           Significance = estimates.sldsc_25bin.h2_z >= 4) %>%
    ggplot + aes(x = ancestry, fill = ancestry, alpha = Significance) +
    pop_fill_scale + scale_y_continuous(labels=comma, name='Number of phenotypes') +
    scale_alpha_discrete(labels=c('TRUE' = expression(paste('h'^'2', 'z'>='4')),
                                  'FALSE' = expression(paste('h'^'2', 'z'<'4')))) +
    geom_histogram(stat='count') + guides(fill=F) +
    xlab('Genetic ancestry group') +
    theme(legend.position=c(1, 1), legend.justification = c(1, 1))
  return(p)
}

h2_data %>%
  filter(!is.na(estimates.sldsc_25bin.h2_z)) %>%
  group_by(ancestry == 'EUR') %>%
  summarize(sig=sum(estimates.sldsc_25bin.h2_z >= 4), n=n())

s_figure_21 = function(output_format = 'png', create_subplots=F) {
  pa = round2_comp()
  pb = round2_comp(T)
  pc = sig_h2_by_ancestry()
  output_type(output_format, paste0('s_figure_21.', output_format), height=3, width=7.5)
  print(ggarrange(pa, pb, pc, nrow = 1, ncol = 3, labels="auto"))
  dev.off()
}
s_figure_21()

# S Figure 22

iter_tab <- read_tsv(paste0('data/final_results_iterations.tsv'))
iter_tab_30 <- read_tsv(paste0('data/final_results_iterations_30rvectors.tsv'))
iter_tab_50 <- read_tsv(paste0('data/final_results_iterations_50rvectors.tsv'))

process_iter = function(iter) {
  iter %>%
    mutate(pheno = str_extract(phenotype_id,'.+(?=-iter[0-9]{1,2})'),
           pheno_trunc = str_trunc(pheno,40,'right'),
           iter_num = str_extract(phenotype_id, '(?<=iter)[0-9]{1,2}$'),
           type = str_split_fixed(pheno,'-',2)[,1]) %>%
    mutate(iter_num = as.numeric(iter_num)) %>%
    group_by(ancestry, pheno) %>%
    summarize(N = n(),
              mean_h2 = mean(h2),
              mean_se = mean(h2_se),
              empirical_sd = sd(h2),
              empirical_se = empirical_sd/sqrt(N),
              se_ratio = empirical_sd / mean_se,
              type = unique(type)) %>% ungroup %>% return
}

combined_iter_tab <- bind_rows(process_iter(iter_tab) %>% mutate(n_iter = 10),
                               process_iter(iter_tab_30) %>% mutate(n_iter = 30),
                               process_iter(iter_tab_50) %>% mutate(n_iter = 50))

run_rhemc = function() {
  p = ggplot(combined_iter_tab) +
    geom_col(aes(x = pheno, y=se_ratio, fill = factor(n_iter),
                 group = factor(n_iter)), width = 0.6,
             position = position_dodge(0.7)) +
    scale_fill_brewer(type="seq", palette = "Blues") +
    geom_hline(yintercept = 1, linetype='dashed') +
    facet_grid(ancestry~type, scales = 'free') +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          strip.background = element_blank(),
          panel.grid = element_blank(),
          legend.position = 'bottom') +
    guides(color=F) +
    xlab('Phenotype') + ylab('Empirical SD\n(normalized by estimator SE)') +
    labs(fill = '# random vectors')
  return(p)
}
s_figure_22 = function(output_format = 'png', create_subplots=F) {
  output_type(output_format, paste0('s_figure_22.', output_format), height=5, width=7.5)
  print(run_rhemc())
  dev.off()
}
s_figure_22()

# S Figure 23

# RHEmc vs round 2
rhemc_vs_round2 = function() {
  h2_data %>% 
    filter(ancestry == 'EUR') %>%
    inner_join(ukb_31063_h2, by=c('code' = 'phenotype')) %>%
    filter(h2_liab_r2 > 0 | estimates.rhemc_8bin.h2_liability > 0) %>%
    ggplot +
    aes(x = h2_liab_r2, y = estimates.rhemc_8bin.h2_liability,
        xmin = h2_liab_r2 - h2_liab_se_r2,
        xmax = h2_liab_r2 + h2_liab_se_r2,
        ymin = estimates.rhemc_8bin.h2_liability - estimates.rhemc_8bin.h2_liability_se,
        ymax = estimates.rhemc_8bin.h2_liability + estimates.rhemc_8bin.h2_liability_se,
        color = trait_type) +
    geom_pointrange() + geom_errorbarh() +
    xlab('S-LDSC heritability (Round 2)') +
    ylab('RHE-mc heritability (Pan-UKB EUR)') +
    coord_cartesian(xlim=c(0,0.6), ylim=c(0,0.6)) +
    trait_color_scale +
    geom_abline(slope = 1, intercept = 0, linetype='dotted') +
    theme(legend.position=c(0.01, 1), legend.justification = c(0, 1),
          legend.background=element_rect(fill = alpha("white", 0)),
          legend.key.size = unit(0.15, 'in'),) %>%
    return
}

s_figure_23 = function(output_format = 'png', create_subplots=F) {
  output_type(output_format, paste0('s_figure_23.', output_format), height=5, width=5)
  print(rhemc_vs_round2())
  dev.off()
}
s_figure_23()

# S Figure 24
heritability_method_comparison = function() {
  methods = c('ldsc', 'sldsc_25bin', 'rhemc_25bin', 'rhemc_8bin', 'rhemc_25bin_50rv') # , 'saige')
  h2_fields = paste0('estimates.', methods, '.h2_observed')
  h2_se_fields = paste0('estimates.',methods,'.h2_observed_se')
  z_fields = paste0('estimates.',methods,'.h2_z')
  measure_se_to_measure = h2_fields
  names(measure_se_to_measure) = h2_se_fields
  method_names = c('LDSC', 'S-LDSC (25 bin)', 'RHE-mc (25 bin)', 'RHE-mc (8 bin)', 'RHE-mc (25 bin, 50 rv)') # , 'saige')
  h2_to_method = method_names
  names(h2_to_method) = h2_fields
  
  phenos <- h2_data %>%
    filter(!is.na(estimates.rhemc_8bin.h2_observed)) %>%
    filter(ancestry == 'EUR') %>%
    pull(phenotype_id)
  h2_mod_for_methods_comparison <- h2_data %>%
    filter(phenotype_id %in% phenos) %>%
    # mutate(estimates.saige.h2_observed = estimates.saige.h2,
    #        estimates.saige.h2_observed_se = NA,
    #        estimates.saige.h2_z = NA) %>%
    filter(ancestry != 'MID') %>%
    left_join(phenotypes)
  
  h2_pivoted <- h2_mod_for_methods_comparison %>% 
    pivot_longer(cols = h2_fields, names_to = 'measure', values_to = 'h2')
  
  z_pivoted <- h2_mod_for_methods_comparison %>% 
    pivot_longer(cols = z_fields, names_to = 'measure', values_to = 'z')
  
  h2_se_pivoted <- h2_mod_for_methods_comparison %>% 
    pivot_longer(cols = h2_se_fields, names_to = 'measure_se', values_to = 'h2_se') %>%
    mutate(measure = measure_se_to_measure[measure_se]) %>%
    left_join(y = h2_pivoted %>% transmute(phenotype_id, ancestry, measure, h2))
  
  trait <- unique(h2_mod_for_methods_comparison[h2_mod_for_methods_comparison$trait_type == 'biomarkers',]$phenotype_id)[1:5]
  
  p = h2_se_pivoted %>% 
    filter(phenotype_id %in% trait) %>%
    mutate(measure = factor(measure, levels = h2_fields, ordered=T)) %>%
    mutate(description2=if_else(grepl(' [a-z]+', description, perl=T), str_replace(description, ' ', '\n'), description)) %>%
    ggplot + aes(x = measure, y=h2, ymin = h2 + h2_se, ymax = h2 - h2_se, color = ancestry) +
    geom_hline(yintercept=0, linetype='dashed') + geom_pointrange(fatten=2.5) +
    facet_grid(description2 ~ ancestry) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.background = element_blank()) +
    coord_cartesian(ylim = c(-0.25, 1)) +
    scale_x_discrete(labels = h2_to_method) + pop_color_scale +
    guides(color='none') + xlab(NULL) + ylab(expression(paste('h'^'2',' point estimates (\u00B1 1 se)')))
  return(p)
}
s_figure_24 = function(output_format = 'png', create_subplots=F) {
  output_type(output_format, paste0('s_figure_24.', output_format), height=6.75, width=5.5)
  p = heritability_method_comparison()
  print(p)
  dev.off()
}
s_figure_24()

# S Figure 25

h2_z_discrete = function() {
  p1 = h2_data %>%
    count(ancestry, h2z = floor(pmax(pmin(estimates.final.h2_z, 8), 0))) %>%
    arrange(desc(h2z)) %>%
    group_by(ancestry) %>%
    mutate(total=cumsum(n)) %>%
    ggplot + aes(x = h2z, y = total, color = ancestry) +
    geom_line() + geom_point() + pop_color_scale +
    scale_y_continuous(labels=comma, name='Number of phenotypes passing cutoff') +
    xlab('Heritability z score cutoff') +
    theme(legend.position = c(1, 1), legend.justification = c(1, 1))
  
  # Rahul's solution
  lapply(0:8, function(z) {
    h2_data %>%
      filter(estimates.final.h2_z >= z) %>%
      group_by(phenotype_id) %>%
      summarize(n_ancestries = n()) %>%
      mutate(z_cutoff = z) %>%
      return()
  }) %>% bind_rows() %>%
    group_by(z_cutoff, n_ancestries) %>%
    summarize(n_traits = n()) %>% ungroup() %>% 
    right_join(y=expand_grid(z_cutoff=0:8, n_ancestries=1:6), by=c('z_cutoff','n_ancestries')) %>%
    mutate(n_anc = factor(n_ancestries, levels=1:6, ordered=T),
           total = ifelse(is.na(n_traits), 0, n_traits), h2z=z_cutoff) %>%
    ggplot + aes(x = h2z, y = total, color = n_anc, group = n_anc) +
    geom_line() + geom_point() +
    scale_y_continuous(labels=comma, name='Number of phenotypes passing cutoff') +
    xlab('Heritability z score cutoff') +
    theme(legend.position = c(1, 1), legend.justification = c(1, 1))
  
  p2 = h2_data %>%
    mutate(h2z = floor(pmax(pmin(estimates.final.h2_z, 8), 0))) %>%
    count(h2z, trait_type, phenocode, pheno_sex, coding, modifier) %>%
    arrange(desc(h2z)) %>%
    group_by(trait_type, phenocode, pheno_sex, coding, modifier) %>%
    mutate(n_anc=cumsum(n)) %>% ungroup %>% count(h2z, n_anc) %>%
    arrange(desc(h2z)) %>%
    group_by(n_anc=as.factor(n_anc)) %>% mutate(total=cumsum(n)) %>%
    ggplot + aes(x = h2z, y = total, color = n_anc, group = n_anc) +
    geom_line() + geom_point() +
    scale_y_continuous(labels=comma, name='Number of phenotypes passing cutoff') +
    xlab('Heritability z score cutoff') +
    scale_color_brewer(name='Number of\nancestries', type="seq", palette = "Spectral") +
    theme(legend.position = c(1, 1), legend.justification = c(1, 1))
  return(list(p1, p2))
}

s_figure_25 = function(output_format = 'png', create_subplots=F) {
  output_type(output_format, paste0('s_figure_25.', output_format), height=3.5, width=7.5)
  ps = h2_z_discrete()
  print(ggarrange(plotlist=ps, nrow = 1, ncol = 2, labels="auto"))
  dev.off()
}
s_figure_25()

# S Figure 27

h2_table_full = function() {
  h2_table_dat = h2_data %>%
    group_by(ancestry) %>%
    summarize(`\nGWAS run\n`=sum(qcflags.GWAS_run),
              `\nDefined heritability estimate\n`=sum(qcflags.defined_h2),
              `\nHeritability z >= 4\n`=sum(qcflags.defined_h2 & qcflags.significant_z),
              `\nRemove AMR due to limited power\n`=sum(if_else(ancestry == 'AMR', 0, qcflags.defined_h2 & qcflags.significant_z)),
              `Heritability within bounds\nfor all ancestries\n`=sum(if_else(ancestry == 'AMR', 0, qcflags.defined_h2 & qcflags.significant_z & qcflags.in_bounds_h2)),
              `\nLambda GC > 0.9 for all ancestries\n`=sum(if_else(ancestry == 'AMR', 0, qcflags.defined_h2 & qcflags.significant_z & qcflags.in_bounds_h2 & qcflags.normal_lambda)),
              `S-LDSC ratio < 0.3 or ratio z < 4\nin all of EUR, CSA, or AFR\n`=sum(if_else(ancestry == 'AMR', 0, qcflags.defined_h2 & qcflags.significant_z & qcflags.in_bounds_h2 & qcflags.normal_lambda & qcflags.normal_ratio)),
              `Passes all filters in EUR and\nat least 1 other ancestry group\n`=sum(if_else(ancestry == 'AMR', 0, qcflags.defined_h2 & qcflags.significant_z & qcflags.in_bounds_h2 & qcflags.normal_lambda & qcflags.normal_ratio & qcflags.EUR_plus_1))
    ) %>% pivot_longer(-ancestry, names_to="Phenotype QC filters") %>%
    pivot_wider(names_from='ancestry') %>%
    rowwise %>%
    mutate(Total = sum(c_across(AFR:MID))) %>%
    select("Phenotype QC filters", "Total", all_of(pops_by_sample_size))
  return(h2_table_dat)
}  
generate_heritability_table_full = function(size = 10) {
  tab = h2_table_full()
  tab2 = tribble(
    ~`Phenotype QC filters`, ~Total, ~EUR, ~CSA, ~AFR, ~EAS, ~MID, ~AMR,
    "Phenotype QC filters", "Total", "EUR", "CSA", "AFR", "EAS", "MID", "AMR"
  )
  tab = union_all(tab2, tab %>% mutate_if(is.double, comma))
  header_cell_color = 'white'
  header_bg_color = "darkgray"
  thm = ttheme(
    padding = unit(c(2, 2), "mm"),
    tbody.style = tbody_style(color = "black", fill = 'white', size=size,
                              hjust = as.vector(matrix(c(0, rep(1, ncol(tab)-1)), ncol = ncol(tab), nrow = nrow(tab), byrow = TRUE)),
                              x = as.vector(matrix(c(0.02, rep(0.9, ncol(tab)-1)), ncol = ncol(tab), nrow = nrow(tab), byrow = TRUE)),
                              parse=TRUE)
  )
  ggtab <- ggtexttable(tab, rows=NULL, cols=NULL, theme = thm)
  ggtab
  ggtab <- table_cell_font(ggtab, row = 1, column = 1, size=size, face = "bold")
  ggtab <- table_cell_bg(ggtab, row = 1, column = 2, fill=header_bg_color, color=header_bg_color)
  ggtab <- table_cell_font(ggtab, row = 1, column = 2, color=header_cell_color, size=size)
  ggtab
  for (i in 2:7) {
    ggtab <- table_cell_bg(ggtab, row = 1, column = i+1,
                           fill=ukb_pop_colors[pops_by_sample_size[i-1]],
                           color=ukb_pop_colors[pops_by_sample_size[i-1]])
    ggtab <- table_cell_font(ggtab, row = 1, column = i+1, color=header_cell_color, size=size, face = "bold")
  }
  ggtab
  return(ggtab)
}

plot_h2_metric_by_ancestry = function(metric = 'estimates.sldsc_25bin.ratio', legend=F) {
  metric_names = c('estimates.sldsc_25bin.ratio' = 'S-LDSC ratio',
                   'lambda_gc' = 'Lambda GC')
  ylim1 = if_else(metric == 'estimates.sldsc_25bin.ratio', 0, 0.75)
  ylim2 = if_else(metric == 'estimates.sldsc_25bin.ratio', 2, 2.5)
  phenotypes %>%
    pivot_longer(all_of(paste0('lambda_gc_', pops_by_sample_size))) %>%
    mutate(ancestry=str_sub(name, start=-3), lambda_gc=value) %>%
    select(all_of(key_fields), 'ancestry', 'lambda_gc') %>%
    inner_join(h2_data) %>%
    filter(qcflags.GWAS_run, qcflags.defined_h2, qcflags.significant_z, 
           qcflags.in_bounds_h2, qcflags.normal_lambda, ancestry != 'AMR') %>%
    ggplot + aes(x = ancestry) + aes_string(y = metric) + 
    geom_violin(aes(fill = ancestry), color=NA, scale = 'width', alpha=0.2) +
    coord_cartesian(ylim=c(ylim1, ylim2)) + pop_fill_scale +
    geom_jitter(aes(color = qcflags.normal_ratio), alpha=0.5, width=0.2, size=0.75) +
    theme(legend.position = 'bottom') + guides(fill=F) +
    ylab(metric_names[[metric]]) -> p
  if (legend) {
    p = p + scale_color_discrete(name='S-LDSC ratio < 0.3 or ratio z < 4\nin all of EUR, CSA, or AFR')
  } else {
    p = p + guides(color=F)
  }
  return(p)
}

s_figure_27 = function(output_format = 'png', create_subplots=F) {
  output_type(output_format, paste0('s_figure_27.', output_format), height=6.5, width=11.5)
  pa = generate_heritability_table_full()
  pb = plot_h2_metric_by_ancestry('lambda_gc', legend=F)
  pc = plot_h2_metric_by_ancestry(legend=T)
  print((pa | (pb / pc)) + plot_annotation(tag_levels = 'a'))
  # print(ggarrange(pa, pb, pc, nrow = 1, ncol = 3, labels="auto", common.legend=T))
  dev.off()
}
s_figure_27()

sig_heritability_by_trait_type = function() {
  eur = h2_data %>%
    filter(ancestry == 'EUR') %>%
    mutate(EUR_z=estimates.final.h2_z >= 4) %>%
    select(all_of(key_fields), 'EUR_z')
  h2_data %>%
    left_join(eur) %>%
    group_by_at(all_of(c(key_fields, 'EUR_z'))) %>%
    summarize(n_ancestries=sum(estimates.final.h2_z >= 4), .groups='drop') %>%
    count(trait_type, EUR_z, n_ancestries) %>%
    filter(n_ancestries > 0 & !is.na(EUR_z)) %>%
    group_by(trait_type, EUR_z) %>%
    arrange(desc(n_ancestries)) %>%
    mutate(n=cumsum(n)) %>% ungroup %>%
    ggplot + aes(x = n_ancestries, y = n, alpha = EUR_z, fill = trait_type) +
    trait_fill_scale + geom_bar(stat='identity', color='black') +
    scale_alpha_discrete(name=expression(paste('EUR S-LDSC Z' >= 4))) +
    facet_wrap(~trait_type, scales='free_y', labeller=as_labeller(trait_type_names)) +
    xlab(expression(paste('Number of ancestries with Z' >= 4))) +
    guides(fill=F) + theme(legend.position = 'bottom') %>%
    return
}

s_figure_28 = function(output_format = 'png', create_subplots=F) {
  output_type(output_format, paste0('s_figure_28.', output_format), height=4.5, width=5.5)
  p = sig_heritability_by_trait_type()
  print(p)
  dev.off()
}
s_figure_28()

upset_by_combo = function() {
  library(UpSetR)
  p = h2_data %>%
    filter(qcflags.pass_all) %>%
    select_at(all_of(c(key_fields, 'ancestry'))) %>%
    mutate(observed = 1) %>%
    pivot_wider(names_from = ancestry, values_from = observed, values_fill = list(observed = 0)) %>%
    as.data.frame %>%
    upset(sets = pops_by_sample_size[1:5], keep.order = TRUE,
          sets.bar.color = ukb_pop_colors[pops_by_sample_size[1:5]])
  
  h2_data %>%
    filter(qcflags.pass_all) %>%
    group_by_at(all_of(key_fields)) %>%
    summarize(n=n()) %>%
    filter(n == 5) %>%
    left_join(phenotypes) %>%
    select(description)
  
  h2_data %>%
    filter(qcflags.pass_all) %>%
    group_by_at(all_of(key_fields)) %>%
    summarize(n=n()) %>% nrow
  
  return(p)
}

s_figure_29 = function(output_format = 'png', create_subplots=F) {
  output_type(output_format, paste0('s_figure_29.', output_format), height=3.5, width=5.5)
  p = upset_by_combo()
  print(p)
  dev.off()
}
s_figure_29()
