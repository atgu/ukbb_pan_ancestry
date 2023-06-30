source('~/ukbb_pan_ancestry/constants.R')
library(googlesheets4)

theme_set(theme_classic())

phenotypes = read_tsv(gzfile('data/phenotype_manifest.tsv.bgz'))
# gs://ukb-diverse-pops-public-free/h2/h2_estimates_all_flat_221123.tsv
h2_data = load_ukb_file('h2_estimates_all_flat_221123.tsv', parent_folder='h2/')

filter_to_pass = function(x, anc) {
  return(filter_at(x, vars(starts_with('phenotype_qc_') & ends_with(anc)), any_vars(. == 'PASS')))
}

compute_neff = function(x, anc) {
  case = paste0('n_cases_', anc)
  control = paste0('n_controls_', anc)
  x[paste0('n_eff_', anc)] = ifelse(is.na(x[[control]]), x[[case]], 4/(1/x[[case]] + 1/x[[control]]))
  return(x)
}

heritability_table = function() {
  tab = tribble(
    ~`Phenotype QC filters`, ~Total, ~EUR, ~CSA, ~AFR, ~EAS, ~MID,
    "Significant and bounded
    heritability", 2679, 937, 1146, 248, 56, 292,
    'Calibrated~\u03BB["GC"]', 2661, 930, 1143, 245, 54, 289,
    "Controlled S-LDSC ratio", 2537, 881, 1098, 234, 51, 273,
    "Passes all filters in EUR &
    ≥1 other genetic ancestry", 1312, 527, 441, 110, 38, 196
  )
  tab2 = tribble(
    ~`Phenotype QC filters`, ~Total, ~EUR, ~CSA, ~AFR, ~EAS, ~MID,
    "Phenotype QC filters", "Total", "EUR", "CSA", "AFR", "EAS", "MID"
  )
  tab = union_all(tab2, tab %>% mutate_if(is.double, as.character))
  return(tab)
}

generate_heritability_table = function(size = 10) {
  tab = heritability_table()
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
  for (i in 2:6) {
    ggtab <- table_cell_bg(ggtab, row = 1, column = i+1,
                           fill=ukb_pop_colors[pops_by_sample_size[i-1]],
                           color=ukb_pop_colors[pops_by_sample_size[i-1]])
    ggtab <- table_cell_font(ggtab, row = 1, column = i+1, color=header_cell_color, size=size, face = "bold")
  }
  ggtab
  return(ggtab)
}

heritability_correlations = function(anc_x='EUR', anc_y='CSA', type='observed',
                                     indep_only=TRUE, remove_questionnaire=TRUE, return_plot=T,
                                     omit_guide=T, omit_type=T) {
  field_x = paste0(if_else(anc_x == 'EUR', 'sldsc_25bin_h2_', 'rhemc_25bin_50rv_h2_'), type)
  field_y = paste0(if_else(anc_y == 'EUR', 'sldsc_25bin_h2_', 'rhemc_25bin_50rv_h2_'), type)
  label_x = paste0('Heritability in ', anc_x, '\n(', if_else(anc_x == 'EUR', 'S-LDSC', 'RHE-mc'), if_else(omit_type, "", paste(";", type)), ')')
  label_y = paste0('Heritability in ', anc_y, '\n(', if_else(anc_y == 'EUR', 'S-LDSC', 'RHE-mc'), if_else(omit_type, "", paste(";", type)), ')')
  x_point = paste0(field_x, '_', anc_x)
  y_point = paste0(field_y, '_', anc_y)
  x_se = paste0(field_x, '_se_', anc_x)
  y_se = paste0(field_y, '_se_', anc_y)
  plot_data = get_phenos(anc_x=anc_x, anc_y=anc_y, indep_only=indep_only, remove_questionnaire=remove_questionnaire)
  print(paste(nrow(plot_data), 'phenotypes remain'))
  if (nrow(plot_data) == 0) {
    if (return_plot) {
      return(NULL) 
    } else {
      return(data.frame(estimate=-1, p=-1, n_phenos=0))
    }
  }
  dat = matrix(c(plot_data[[x_point]], plot_data[[x_se]], plot_data[[y_point]], plot_data[[y_se]]), ncol=4)
  york_res = york(dat)
  print(paste("York slope =", york_res$b[[1]], "; intercept =", york_res$a[[1]], "; p =", york_res$p.value))
  res = cor.test(plot_data[[x_point]], plot_data[[y_point]], method='spearman')
  p = plot_data %>%
    filter(rhemc_25bin_50rv_h2_liability_CSA + rhemc_25bin_50rv_h2_liability_se_CSA < 1 &
             sldsc_25bin_h2_liability_EUR + sldsc_25bin_h2_liability_se_EUR < 1) %>%
    ggplot + aes(x = .data[[x_point]], y = .data[[y_point]],
                 xmin = .data[[x_point]] - .data[[x_se]],
                 xmax = .data[[x_point]] + .data[[x_se]],
                 ymin = .data[[y_point]] - .data[[y_se]],
                 ymax = .data[[y_point]] + .data[[y_se]],
                 text = description, color=trait_type) +
    geom_pointrange() + geom_errorbarh() +
    xlab(label_x) + ylab(label_y) + trait_color_scale +
    geom_abline(slope=1, intercept=0, linetype='dotted') +
    geom_abline(slope=york_res$b[[1]], intercept=york_res$a[[1]], linetype='dashed') +
    ylim(c(0, NA))
  if (return_plot) {
    if (omit_guide) {
      p = p + guides(color=F)
    }
    return(p)
  } else {
    return(data.frame(estimate=york_res$b[[1]], p=york_res$p.value, n_phenos=nrow(plot_data)))
    # return(data.frame(estimate=res$estimate[[1]], p=res$p.value, n_phenos=nrow(plot_data)))
  }
}
heritability_correlations(type='observed', remove_questionnaire = F)

get_phenos = function(anc_x='EUR', anc_y='CSA', indep_only=TRUE, remove_questionnaire=FALSE) {
  phenotypes %>% 
    filter_to_pass(anc_x) %>%
    filter_to_pass(anc_y) %>%
    # filter(num_pops_pass_qc > 2) %>%
    filter(!remove_questionnaire | !grepl('question', description_more, ignore.case = TRUE)) %>%
    filter(!indep_only | in_max_independent_set) %>%
    return
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
    geom_histogram(bins=200, alpha=0.5, position='identity') +
    geom_vline(xintercept=4, linetype='dashed') +
    theme(legend.position = c(0.9, 0.9)) +
    scale_fill_discrete(name=NULL) +
      xlab('Heritability Z-score') + ylab('Number of\nancestry-trait pairs') -> p
  return(p)
}
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
    xlab('S-LDSC heritability') +
    ylab('RHE-mc heritability') +
    trait_color_scale +
    theme(legend.position = c(0.8, 0.8),
          legend.key.size = unit(legend_size, 'in')
          ) +
    geom_abline(slope = res$coefficients[[2]], intercept = res$coefficients[[1]], linetype='dashed') +
    geom_abline(slope = 1, intercept = 0, linetype='dotted') -> p
  return(p)
}

all_pairs = map_dfr(pops, function(x) {
  map_dfr(pops, function(y) heritability_correlations(x, y, remove_questionnaire = F, return_plot=F) %>%
    mutate(pop1=x, pop2=y)
  )
})
all_pairs %>%
  filter(n_phenos > 2 & pop1 < pop2)

phenotypes %>%
  filter_to_pass('EUR') %>%
  ggplot + aes(x = n_eff_EUR, y = 1 / (sldsc_25bin_h2_liability_se_EUR ^ 2)) + geom_point()

figure2 = function(output_format = 'png') {
  p2a = rhemc_v_sldsc(legend_size=0.1)
  p2b = heritability_histograms()
  p2c = generate_heritability_table(8.5)
  p2d = heritability_correlations(remove_questionnaire = F)
  output_type(output_format, paste0('figure2.', output_format), height=5, width=7.5)
  print(ggarrange(p2a, p2b, p2c, p2d, nrow=2, ncol = 2, labels="auto"))
  dev.off()
}
figure2()
