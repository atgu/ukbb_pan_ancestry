source('constants.R')

pheno_summary = load_ukb_file('pheno_summary_full.txt.bgz', 'phenos/')

# 1c: full distribution?

by_pop = pheno_summary %>%
  count(pheno_data.pop) %>%
  transmute(pop=fct_reorder(pheno_data.pop, -n), n=n)

n_phenos_by_pop = function(table_size = 2.5) {
  inset_table = pheno_summary %>%
    group_by(pheno_data.pop) %>%
    summarize(n_samples=max(pheno_data.n_cases + pheno_data.n_controls, na.rm=T)) %>% ungroup %>%
    arrange(desc(n_samples)) %>%
    transmute(`Genetic ancestry`=ukb_pop_names[pheno_data.pop], `Num. individuals`=n_samples)

  p = by_pop %>%
    ggplot + aes(x = pop, y = n, fill = pop) + 
    geom_bar(stat='identity') +
    ylab('Number of phenotypes') +
    xlab(NULL) +
    scale_y_continuous(labels=comma) + 
    scale_fill_manual(values=ukb_pop_colors, guide=F) +
    # scale_x_discrete(labels=ukb_pop_names) +
    # theme(axis.text.x = element_text(angle = 35, hjust=1)) + 
    geom_text(aes(label=comma(n)), vjust=0, nudge_y = 80, size=2.5) +
    annotate(geom = "table", x = 6.5, y = max(by_pop$n), label = list(inset_table), 
             vjust = 1, hjust = 1, table.theme = ttheme_gtplain, table.hjust = 1,
             size=table_size)
  return(p)
}

n_phenos_by_pop_combo = function(legend_size=0.15) {
  # multi_pop_names = c('Any', paste('Any', 2:5), 'All 6')
  multi_pop_names = c(paste0(1:5, '+'), '6')
  names(multi_pop_names) = as.character(1:6)
  
  by_pop_cumulative = pheno_summary %>%
    group_by_at(key_fields) %>%
    dplyr::count(name='n_pops') %>% ungroup %>%
    dplyr::count(n_pops) %>%
    arrange(desc(n_pops)) %>%
    mutate(n=cumsum(n)) %>% ungroup %>%
    arrange(n_pops) %>%
    mutate(pop=as.factor(n_pops))
  
  trait_type_order = pheno_summary %>%
    select_at(key_fields) %>% distinct %>%
    dplyr::count(trait_type) %>%
    arrange(desc(n)) %>%
    pull(trait_type) %>% rev
  
  by_pop_cumulative_by_trait = pheno_summary %>%
    group_by_at(key_fields) %>%
    dplyr::count(trait_type, name='n_pops') %>% ungroup %>%
    dplyr::count(trait_type, n_pops) %>%
    group_by(trait_type) %>%
    arrange(desc(n_pops)) %>%
    mutate(n=cumsum(n)) %>% ungroup %>%
    arrange(n_pops) %>%
    mutate(trait_type=fct_relevel(trait_type, trait_type_order))
  
  p = by_pop_cumulative_by_trait %>%
    ggplot + aes(x = as.character(n_pops), y = n, fill = trait_type) + 
    geom_bar(stat='identity') +
    ylab('Number of phenotypes') +
    scale_fill_manual(values=trait_type_colors, name='Trait type', labels=trait_type_names) +
    scale_x_discrete(labels=multi_pop_names) +
    scale_y_continuous(labels=comma) + 
    xlab('Number of genetic ancestries') +
    theme(# axis.text.x = element_text(angle = 35, hjust=1),
          legend.key.size = unit(legend_size, 'in')) + 
    geom_text(aes(label=comma(n)), vjust=0, nudge_y = 80, size=2.5,
              data=by_pop_cumulative %>% mutate(trait_type=NA_character_)) +
    theme(legend.position=c(1, 1), legend.justification = c(1, 1))
  return(p)
}

sig_hits_by_pheno = load_ukb_file('sig_hits_pops_by_pheno_full.txt.bgz', 'sig_hits/', use_local=F)

get_sig_hits_data = function(filter_to_6pop=T) {
  run_phenos = pheno_summary %>%
    select_at(c(key_fields, 'pheno_data.pop')) %>% rename(pop=pheno_data.pop)
  phenos_6pop = run_phenos %>% group_by_at(key_fields) %>% 
    count() %>% filter(n == 6) %>%
    ungroup %>% select(-n)
  
  plot_data = sig_hits_by_pheno %>%
    inner_join(run_phenos)
  
  if (filter_to_6pop) {
    plot_data = plot_data %>%
    inner_join(phenos_6pop)
  }
  plot_data = plot_data %>%
    mutate(sig_pops_by_pheno_clumped=if_else(is.na(clumped), 0, clumped),
           sig_pops_by_pheno_total=if_else(is.na(total), 0, total)) %>%
    group_by(pop) %>%
    mutate(sig_hits_clumped_percentile=ntile(sig_pops_by_pheno_clumped, 100),
           sig_hits_total_percentile=ntile(sig_pops_by_pheno_total, 100),
    ) %>% ungroup
  return(plot_data)
}

get_sig_hits_data() %>%
  group_by(pop) %>%
  summarize(mean_number_clumped_hits=mean(sig_pops_by_pheno_clumped),
            per10_number_clumped_hits=quantile(sig_pops_by_pheno_clumped, 0.9),
            per25_number_clumped_hits=quantile(sig_pops_by_pheno_clumped, 0.75)) %>%
  group_by(pop == 'EUR') %>%
  summarize_if(is.numeric, sum)

sig_hits_cdf = function(only_clumped=T, include_amr=F) {
  plot_data = get_sig_hits_data()
  
  means = plot_data %>%
    group_by(pop) %>%
    summarize(mean_clumped_hits=mean(sig_pops_by_pheno_clumped)) %>%
    left_join(plot_data) %>%
    mutate(diff=abs(sig_pops_by_pheno_clumped-mean_clumped_hits)) %>%
    group_by(pop) %>%
    slice(which.min(diff))
  
  p = plot_data %>%
    filter(include_amr | (pop != 'AMR')) %>%
    ggplot + aes_string(x = if_else(only_clumped, 'sig_pops_by_pheno_clumped', 'sig_pops_by_pheno_total'),
                        y = if_else(only_clumped, 'sig_hits_clumped_percentile', 'sig_hits_total_percentile')) +
    aes(group = pop, color = pop) +
    geom_line(lwd=1) + ylab('Phenotype percentile') +
    # geom_point(aes(x = mean_clumped_hits), data=means) +
    scale_x_log10(name=paste0('Number of ', if_else(only_clumped, 'independent ', ''), 'significant associations'), label=comma) + 
    scale_color_manual(values=ukb_pop_colors, name=NULL) + 
    theme(legend.position = c(1, 0.01), legend.justification = c(1, 0),
          plot.margin=margin(r=22))
  
  # plot_data %>%
  #   ggplot + aes(x = sig_pops_by_pheno_total, y = sig_hits_total_percentile, group = pop, color = pop) +
  #   geom_line(lwd=1) + ylab('Phenotype percentile') +
  #   scale_x_log10(name='Number of significant hits for phenotype') +
  #   pop_color_scale
  return(p)
}

n_hits_per_pop = function(include_amr=F) {
  plot_data = get_sig_hits_data()
  
  p = plot_data %>%
    filter(include_amr | (pop != 'AMR')) %>%
    group_by(pop) %>%
    summarize(mean_clumped_hits=mean(sig_pops_by_pheno_clumped)) %>%
    left_join(by_pop) %>%
    mutate(pop=fct_reorder(pop, -n)) %>%
    # TODO: harmonize qc data
    # filter(pheno_data.heritability.qcflags.pass_all & pheno_data.pop != 'AMR') %>%
    ggplot + aes(x = pop, y = mean_clumped_hits, fill = pop) +
    geom_bar(stat='identity', color=NA) + scale_y_continuous(label=comma) +
    # scale_y_break(c(5, 288)) + 
    # scale_y_log10() +
    scale_y_continuous(trans = locusviz::trans_loglog_p(2), breaks=c(0, 1, 2, 5, 10, 20, 50, 100, 200)) +
    pop_fill_scale +
    xlab('Genetic ancestry group') + ylab('Mean number of independent\nsignificant associations per phenotype') +
    guides(fill='none')
  return(p)
}

n_eff_vs_sig_hits = function() {
  p = sig_hits_by_pheno %>%
    inner_join(run_phenos) %>%
    inner_join(phenos_6pop) %>%
    mutate(sig_pops_by_pheno_clumped=if_else(is.na(sig_pops_by_pheno_clumped), 0, sig_pops_by_pheno_clumped),
           sig_pops_by_pheno_total=if_else(is.na(sig_pops_by_pheno_total), 0, sig_pops_by_pheno_total),
           n_eff=if_else(!is.na(n_controls), 4 / (1 / n_cases + 1 / n_controls), n_cases)) %>%
    group_by(pop, n_eff_binned=10^(round(log10(n_eff)))) %>%
    summarize(sig_pops_by_pheno_clumped=mean(sig_pops_by_pheno_clumped)) %>%
    ggplot + aes(x = n_eff_binned, y = sig_pops_by_pheno_clumped, color = pop) +
    geom_line(lwd=1) + scale_x_log10(name='Effective sample size', labels=comma) +
    scale_y_log10(name='Mean number of significant\nclumped hits',
                  labels=function(x) if_else(x < 1, comma(x, 0.1), comma(x, 1))) +
    guides(color=F) + pop_color_scale +
    theme(legend.position = c(1, 0.01), legend.justification = c(1, 0))
  return(p)
}


figure1 = function(output_format = 'png', create_subplots=F) {
  table_size = if_else(create_subplots, 4, 2.75)
  p1a = n_phenos_by_pop(table_size)
  p1b = n_phenos_by_pop_combo()
  # p1c = n_hits_per_pop(include_amr=T) + guides(fill='none')
  p1c = sig_hits_cdf(include_amr=T) + guides(color='none')
  if (create_subplots) {
    fig1_plots = list(p1a, p1b, p1c)
    for (i in seq_along(fig1_plots)) {
      output_type(output_format, paste0('figure1_panel', i, '.', output_format),
                  height=3.5, width=5)
      print(fig1_plots[i])
      dev.off()
    }
  } else {
    output_type(output_format, paste0('figure1.', output_format), height=3.25, width=10)
    print(ggarrange(p1a, p1b, p1c, nrow = 1, ncol = 3, labels="auto"))
    # print((p1a + p1b) / (p1c + p1d))
    dev.off()
  }
}
figure1()
figure1(create_subplots = T)



