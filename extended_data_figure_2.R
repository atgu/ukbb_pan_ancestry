source('figure2.R')

rhemc_v_sldsc = function(legend_size = 0.25) {
  rhe_sldsc_data = h2_data %>%
    filter(estimates.sldsc_25bin.h2_liability > 0 | estimates.rhemc_8bin.h2_liability > 0) %>%
    filter(ancestry == 'EUR' & !is.na(estimates.rhemc_8bin.h2_z))
  
  
  # dat = matrix(c(plot_data[[x_point]], plot_data[[x_se]], plot_data[[y_point]], plot_data[[y_se]]), ncol=4)
  # york_res = york(dat)
  # print(paste("York slope =", york_res$b[[1]], "; intercept =", york_res$a[[1]], "; p =", york_res$p.value))
  res = lm(estimates.rhemc_8bin.h2_liability ~ estimates.sldsc_25bin.h2_liability, rhe_sldsc_data)
  print(summary(res))
  rhe_sldsc_data %>%
    ggplot +
    aes(x = estimates.sldsc_25bin.h2_liability, y = estimates.rhemc_8bin.h2_liability,
        xmin = estimates.sldsc_25bin.h2_liability - estimates.sldsc_25bin.h2_liability_se,
        xmax = estimates.sldsc_25bin.h2_liability + estimates.sldsc_25bin.h2_liability_se,
        ymin = estimates.rhemc_8bin.h2_liability - estimates.rhemc_8bin.h2_liability_se,
        ymax = estimates.rhemc_8bin.h2_liability + estimates.rhemc_8bin.h2_liability_se,
        color = trait_type) +
    geom_pointrange() + geom_errorbarh() +
    xlab('S-LDSC heritability (EUR)') +
    ylab('RHE-mc heritability (EUR)') +
    trait_color_scale +
    theme(legend.position = c(0.8, 0.8),
          legend.key.size = unit(legend_size, 'in'),
          # legend.background=element_rect(fill = alpha("white", 0))
    ) +
    geom_abline(slope = res$coefficients[[2]], intercept = res$coefficients[[1]], linetype='dashed') +
    geom_abline(slope = 1, intercept = 0, linetype='dotted') -> p
  return(p)
}

heritability_histograms = function() {
  upper_bound = 10
  lower_bound = -5
  h2_data %>%
    filter(ancestry != 'EUR') %>%
    select(`S-LDSC`=estimates.ldsc.h2_z, `RHE-mc`=estimates.final.h2_z) %>%
    pivot_longer(cols = everything()) %>%
    mutate(value=case_when(value < lower_bound ~ lower_bound,
                           value > upper_bound ~ upper_bound,
                           TRUE ~ value)) %>%
    ggplot + aes(x = value, group = name, fill = name) +
    geom_histogram(bins=100, alpha=0.5, position='identity') +
    geom_vline(xintercept=4, linetype='dashed') +
    theme(legend.position = c(0.8, 0.9)) +
    scale_fill_manual(name=NULL, values=c('RHE-mc' = muted('red', l=50),
                                          'S-LDSC' = muted('blue', l=50))) +
    xlab('Heritability Z-score') + ylab('Number of ancestry-trait pairs') -> p
  # p
  return(p)
}

p1 = rhemc_v_sldsc(legend_size=0.1)
p2 = heritability_histograms()
p3 = heritability_correlations(filter_to_pass=F, filter_EUR_z=T)
extended_data_figure2 = function(output_format = 'png') {
  output_type(output_format, paste0('extended_data_figure2.', output_format), height=3, width=9.5)
  print(ggarrange(p1, p2, p3, ncol = 3, labels='auto'))
  dev.off()
}
extended_data_figure2()
