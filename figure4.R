source('~/ukbb_pan_ancestry/constants.R')
library(patchwork)
library(gganimate)
library(gifski)

meta_eur_comparison_all = load_ukb_file('meta_eur_comparison.tsv.bgz', subfolder = 'known_novel/', force_cols = cols(coding=col_character()))
meta_eur_comparison_all %>% filter(N_pops > 1) %>% mutate(
  # nlog10p(meta) - nlog10p(EUR) = log10(P[EUR]/P[Meta])
  ratio = Pvalue_meta - Pvalue_EUR,
  freq_max = pmax(freq_AFR, freq_AMR, freq_CSA, freq_EAS, freq_MID),
  sig_in=case_when(Pvalue_meta > log_threshold & Pvalue_EUR < log_threshold ~ 'meta_only',
                   Pvalue_meta < log_threshold & Pvalue_EUR > log_threshold ~ 'EUR_only',
                   Pvalue_meta > log_threshold & Pvalue_EUR > log_threshold ~ 'both',
                   TRUE ~ 'neither'),
  orig_sig_in=case_when(Pvalue_meta > log_orig_threshold & Pvalue_EUR < log_orig_threshold ~ 'meta_only',
                        Pvalue_meta < log_orig_threshold & Pvalue_EUR > log_orig_threshold ~ 'EUR_only',
                        TRUE ~ 'both')
) %>% filter(Pvalue_meta > log_orig_threshold | Pvalue_EUR > log_orig_threshold) -> meta_eur_comparison


power_curve = function(return_plots=F, point_size = 1, EUR_fraction_to_sample = 0.01, meta_fraction_to_sample=0.5, clump_only=F,
                       x_shape=F) {
  staging = meta_eur_comparison %>%
    filter(modifier == 'irnt') %>%
    mutate(beta_EUR = abs(beta_EUR),
           beta_meta = abs(beta_meta),
           MAF_EUR = 0.5 - abs(0.5 - freq_EUR),
           MAF_max = 0.5 - abs(0.5 - freq_max))
  
  if (clump_only) {
    staging %<>% filter(!is.na(in_clump))
    EUR_fraction_to_sample = 1
    meta_fraction_to_sample = 1
  }
  
  staging %>% count(trait_type, phenocode, pheno_sex, coding, modifier) %>% count
  
  set.seed(663)
  plot_data = staging %>%
    filter(sig_in %in% c('EUR_only', 'both')) %>%
    sample_frac(EUR_fraction_to_sample) %>%
    union_all(staging %>% filter(sig_in == 'meta_only') %>% sample_frac(meta_fraction_to_sample))
  
  settings = list(ylim(0, 1.15), xlab('MAF'), ylab('Beta'),
                  theme(legend.position = c(1, 1), legend.justification = c(1, 1)))
  p1 = plot_data %>%
    filter(sig_in != 'meta_only') %>%
    ggplot + aes(x = MAF_EUR, y = beta_EUR) +
    geom_point(color=color_eur, size = point_size) +
    settings + xlab('EUR MAF') + ylab('EUR beta')

  p2 = plot_data %>%
    filter(sig_in != 'meta_only') %>%
    ggplot + aes(x = MAF_EUR, y = beta_EUR) +
    geom_point(color=color_eur, size = point_size) +
    geom_point(aes(x = MAF_EUR, y = beta_EUR, color=max_maf_pop), size = point_size,
               data = plot_data %>% filter(sig_in == 'meta_only' & Pvalue_meta > 12)) +
    settings + xlab('EUR MAF') + ylab('EUR beta')
  
  p2_black = plot_data %>%
    filter(sig_in != 'meta_only') %>%
    ggplot + aes(x = MAF_EUR, y = beta_EUR, color=color_eur) +
    geom_point(size = point_size) + 
    geom_point(aes(x = MAF_EUR, y = beta_EUR, color='black'), size = point_size,
               data = plot_data %>% filter(sig_in == 'meta_only' & Pvalue_meta > 12)) +
    settings + xlab('EUR MAF') + ylab('EUR beta') + scale_color_identity()
  
  data2 = plot_data %>% filter(sig_in == 'meta_only' & Pvalue_meta > 12) %>% 
    mutate(x = chrom == 'X')
  p3 = plot_data %>%
    filter(sig_in != 'meta_only') %>%
    ggplot + aes(x = MAF_EUR, y = beta_EUR) +
    geom_point(color=color_eur, size = point_size)
  
  if (x_shape) {
    p3 = p3 + geom_point(aes(x = MAF_max, y = beta_meta, color=max_maf_pop), size = point_size, data = data2)
  } else {
    p3 = p3 + geom_point(aes(x = MAF_max, y = beta_meta, color=max_maf_pop,
                             shape = x), size = point_size, data = data2) +
      guides(shape=F)
  }
  p3 = p3 + settings + xlab('Ancestry-specific MAF') + ylab('Meta-analysis beta') + pop_color_scale
  p4 = plot_data %>%
    filter(sig_in != 'meta_only') %>%
    ggplot + aes(x = MAF_EUR, y = beta_EUR) +
    geom_point(color=color_eur, size = point_size) +
    geom_segment(aes(x = MAF_EUR, y = beta_EUR, xend=MAF_max, yend = beta_meta, color=max_maf_pop),
                 arrow=arrow(length=unit(0.05, "inches")),
                 data = plot_data %>% filter(sig_in == 'meta_only' & Pvalue_meta > 12)) +
    settings + xlab('Ancestry-specific MAF') + ylab('Meta-analysis beta') + pop_color_scale
  
  if (return_plots) {
    return(list(p1, p2, p2_black, p3, p4))
  }
  width = 4.5
  height = 3.5
  png('EUR_discovery.png', width = width, height = height, res = 300, units='in')
  print(p1)
  dev.off()
  png('meta_discovery_in_EUR.png', width = width, height = height, res = 300, units='in')
  print(p2)
  dev.off()
  png('meta_discovery_in_EUR_anon.png', width = width, height = height, res = 300, units='in')
  print(p2_black)
  dev.off()
  png('meta_discovery_in_meta.png', width = width, height = height, res = 300, units='in')
  print(p3)
  dev.off()
  
  before_state = 1
  after_state = 2
  test = plot_data %>% filter(sig_in != 'meta_only') %>% mutate(state=before_state, max_maf_pop = 'EUR') %>%
    union_all(plot_data %>% filter(sig_in != 'meta_only') %>% mutate(state=after_state, max_maf_pop = 'EUR')) %>%
    union_all(plot_data %>% filter(sig_in == 'meta_only' & Pvalue_meta > 12) %>% mutate(state=before_state)) %>%
    union_all(plot_data %>% filter(sig_in == 'meta_only' & Pvalue_meta > 12) %>% mutate(state=after_state) %>%
                mutate(MAF_EUR = MAF_max, beta_EUR = beta_meta))

  # possible to change from EUR color to other color by creating more states and interpolating manually
  xlabel = c('1' = 'EUR MAF', '2' = 'Ancestry-specific MAF')
  ylabel = c('1' = 'EUR beta', '2' = 'Meta-analysis beta')
  p = test %>%
    ggplot + aes(x = MAF_EUR, y = beta_EUR, color=max_maf_pop) +
    geom_point() +
    transition_states(state, wrap = F) +
    settings + 
    ease_aes('quadratic-out') +
    labs(x = '{xlabel[[closest_state]]}', y = '{ylabel[[closest_state]]}')
  
  ren = animate(p, duration = 3, fps = 10, width = width, height = height, res = 300, units='in',
                renderer = gifski_renderer(loop = FALSE), rewind=FALSE)
  anim_save("meta_discovery.gif", ren)
  ren_loop = animate(p, duration = 3, fps = 10, width = width, height = height, res = 300, units='in',
                renderer = gifski_renderer(loop = TRUE), rewind=FALSE)
  anim_save("meta_discovery_loop.gif", ren_loop)
}
p = power_curve(return_plots=T, point_size = 0.8)
p2 = power_curve(return_plots=T, point_size = 0.8, clump_only=T)

figure4 = function(output_format = 'png', create_subplots=F, clump_only=F) {
  power_plots = power_curve(return_plots=T, point_size = 0.8, meta_fraction_to_sample = 1, clump_only=clump_only)
  p4a = power_plots[[3]] +
    scale_color_identity(breaks = c("black", color_eur),
                         labels = c("Discovered only in meta-analysis", "Discovered only in EUR"),
                         guide='legend') +
    theme(legend.key.size = unit(0.05, 'in'), legend.position=c(1, 1.25)) +
    guides(color=guide_legend(nrow=2, byrow=TRUE, title=NULL))
  p4b = power_plots[[4]] +
    theme(legend.key.size = unit(0.05, 'in'), legend.position=c(1, 1.25)) +
    guides(color=guide_legend(nrow=2, byrow=TRUE, title=NULL))
  
  # table_size = if_else(create_subplots, 4, 2)
  
  if (create_subplots) {
    fig4_plots = list(p4a, p4b)
    for (i in seq_along(fig4_plots)) {
      output_type(output_format, paste0('figure4_panel', i, if_else(clump_only, '.clump', ''), '.', output_format),
                  height=4, width=6)
      print(fig4_plots[i])
      dev.off()
    }
  } else {
    output_type(output_format, paste0('figure4.', if_else(clump_only, 'clump.', ''), output_format), height=3, width=3.75)
    # print(ggarrange(p1a, p1b, p1c, p1d, nrow = 2, ncol = 2))
    print((p4a / p4b) + plot_annotation(tag_levels = 'a'))
    dev.off()
  }
}
figure4()
figure4(clump_only=T)
figure4(create_subplots = T)
