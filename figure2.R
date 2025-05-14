source('~/ukbb_pan_ancestry/constants.R')
library(googlesheets4)
library(IsoplotR)

theme_set(theme_classic())

phenotypes = load_ukb_file('phenotype_manifest.tsv.bgz', subfolder='phenos/')
# gs://ukb-diverse-pops-public-free/h2/h2_estimates_all_flat_221123.tsv
h2_data = load_ukb_file('h2_estimates_all_flat_221123.tsv', parent_folder='h2/') %>%
  mutate(phenotype_id = paste_noempty(trait_type, phenocode, pheno_sex, coding, modifier, sep='-')) %>%
  mutate(code2 = str_remove(str_extract(phenotype_id, '(?=.+)-[A-Za-z0-9_]+'),'^-')) %>% 
  mutate(post_str = str_extract(phenotype_id, '(?<=[A-Za-z0-9]{1,30}-[A-Za-z0-9]{1,30}-(both_sexes|female|male)-)[A-Za-z0-9_]+')) %>%
  mutate(code = ifelse(!is.na(post_str), paste0(code2,'_',post_str), code2))

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
  tab2 = tribble(
    ~`Phenotype QC filters`, ~Total, ~EUR, ~CSA, ~AFR, ~EAS, ~MID,
    "Phenotype QC filters", "Total", "EUR", "CSA", "AFR", "EAS", "MID"
  )
  h2_table_dat = h2_data %>%
    filter(ancestry != 'AMR') %>%
    group_by(ancestry) %>%
    summarize(`Significant and bounded\nheritability`=sum(qcflags.defined_h2 & qcflags.significant_z & qcflags.in_bounds_h2),
              `lambda`=sum(qcflags.defined_h2 & qcflags.significant_z & qcflags.in_bounds_h2 & qcflags.normal_lambda),
              `Controlled S-LDSC ratio`=sum(qcflags.defined_h2 & qcflags.significant_z & qcflags.in_bounds_h2 & qcflags.normal_lambda & qcflags.normal_ratio),
              `Passes all filters in EUR and\n≥1 other ancestry group`=sum(qcflags.pass_all)
    ) %>% pivot_longer(-ancestry, names_to="Phenotype QC filters") %>%
    pivot_wider(names_from='ancestry') %>%
    rowwise %>%
    mutate(Total = sum(c_across(AFR:MID))) %>%
    select("Phenotype QC filters", "Total", all_of(pops_by_sample_size[1:5])) %>%
    mutate(`Phenotype QC filters` = if_else(`Phenotype QC filters` == 'lambda', 'Calibrated~\u03BB["GC"]', `Phenotype QC filters`))
  tab = union_all(tab2, h2_table_dat %>% mutate_if(is.numeric, comma))
  return(tab)
}

phenotypes %>%
  filter(grepl('AFR', pops_pass_qc) & grepl('CSA', pops_pass_qc) & grepl('EUR', pops_pass_qc)) %>%
  View

# polygenicity %>% filter(grepl('30600', Pheno) | grepl('30620', Pheno) | grepl('30870', Pheno) | grepl('30890', Pheno) | grepl('20002-both_sexes-1226', Pheno) | grepl('30080', Pheno) | grepl('30100', Pheno) | grepl('30270', Pheno) | grepl('30300', Pheno) | grepl('3063', Pheno))

h2_data %>%
  group_by(ancestry) %>%
  summarize(significant_phenos=sum(estimates.final.h2_z > 4, na.rm=T))

h2_data %>%
  filter(estimates.final.h2_z > 4) %>%
  group_by_at(key_fields) %>%
  summarize(n=n()) %>%
  ungroup %>% count(n) %>%
  mutate(prop=nn/sum(nn))

h2_data %>%
  filter(qcflags.pass_all) %T>% 
  {print(nrow(.))} %>%
  group_by_at(key_fields) %>%
  summarize(n=n()) %>%
  ungroup %>% count(n) %>%
  mutate(prop=nn/sum(nn))

phenotypes %>% filter(num_pops_pass_qc > 1 & in_max_independent_set) %T>%
  {print(nrow(.))} %$%
  sum(num_pops_pass_qc)

h2_data %>%
  filter(qcflags.pass_all) %>%
  group_by(trait_type) %>%
  summarize(mean_observed=mean(estimates.final.h2_observed),
            mean_EUR_observed=sum(estimates.final.h2_observed*(ancestry == 'EUR'))/sum(ancestry == 'EUR'))

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

get_phenos = function(anc_x='EUR', anc_y='CSA', indep_only=TRUE, remove_questionnaire=FALSE,
                      filter_to_pass=TRUE, filter_EUR_z=FALSE) {
  if (filter_to_pass) {
    dat = phenotypes %>% 
      filter_to_pass(anc_x) %>%
      filter_to_pass(anc_y)
  } else {
    field_x = paste0(if_else(anc_x == 'EUR', 'sldsc_25bin_h2_', 'rhemc_25bin_50rv_h2_'), 'observed_', anc_x)
    field_y = paste0(if_else(anc_y == 'EUR', 'sldsc_25bin_h2_', 'rhemc_25bin_50rv_h2_'), 'observed_', anc_y)
    dat = phenotypes %>%
      filter(!is.na(get(field_x)) & !is.na(get(field_y)))
  }
  dat %>%
    filter(!filter_EUR_z | (sldsc_25bin_h2_z_EUR > 4)) %>%
    filter(!remove_questionnaire | !grepl('question', description_more, ignore.case = TRUE)) %>%
    filter(!indep_only | in_max_independent_set) %>%
    return
}

heritability_correlations = function(anc_x='EUR', anc_y='CSA', type='observed',
                                     indep_only=TRUE, remove_questionnaire=F,
                                     filter_to_pass=TRUE, filter_EUR_z=FALSE,
                                     return_plot=T, omit_guide=T, omit_type=T, ymax=NA) {
  field_x = paste0(if_else(anc_x == 'EUR', 'sldsc_25bin_h2_', 'rhemc_25bin_50rv_h2_'), type)
  field_y = paste0(if_else(anc_y == 'EUR', 'sldsc_25bin_h2_', 'rhemc_25bin_50rv_h2_'), type)
  label_x = paste0('Heritability in ', anc_x, '\n(', if_else(anc_x == 'EUR', 'S-LDSC', 'RHE-mc'), if_else(omit_type, "", paste(";", type)), ')')
  label_y = paste0('Heritability in ', anc_y, '\n(', if_else(anc_y == 'EUR', 'S-LDSC', 'RHE-mc'), if_else(omit_type, "", paste(";", type)), ')')
  x_point = paste0(field_x, '_', anc_x)
  y_point = paste0(field_y, '_', anc_y)
  x_se = paste0(field_x, '_se_', anc_x)
  y_se = paste0(field_y, '_se_', anc_y)
  plot_data = get_phenos(anc_x=anc_x, anc_y=anc_y, indep_only=indep_only,
                         remove_questionnaire=remove_questionnaire,
                         filter_to_pass=filter_to_pass, filter_EUR_z=filter_EUR_z)
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
  # res = cor.test(plot_data[[x_point]], plot_data[[y_point]])
  # print(res)
  # res = cor.test(plot_data[[x_point]], plot_data[[y_point]], method='spearman')
  # print(res)
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
    coord_cartesian(ylim=c(0, ymax))
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
heritability_correlations()

output_type('png', 'eur_afr_heritability.png', height=3.5, width=4.5)
print(heritability_correlations(filter_to_pass=T, anc_y='AFR'))
dev.off()
# all_pairs = map_dfr(pops, function(x) {
#   map_dfr(pops, function(y) heritability_correlations(x, y, remove_questionnaire = F, return_plot=F) %>%
#     mutate(pop1=x, pop2=y)
#   )
# })
# all_pairs %>%
#   filter(n_phenos > 2 & pop1 < pop2)

# phenotypes %>%
#   filter_to_pass('EUR') %>%
#   compute_neff('EUR') %>%
#   ggplot + aes(x = n_eff_EUR, y = 1 / (sldsc_25bin_h2_liability_se_EUR ^ 2)) + geom_point()

passing_traits = function() {
  h2_data %>%
    filter(qcflags.pass_all) %>%
    select_at(key_fields) %>%
    distinct %>%
    count(trait_type) %>%
    ggplot + aes(x = trait_type, y = n, fill=trait_type) + 
    geom_bar(stat='identity') +
    trait_fill_scale + ylab('Number of traits') +
    scale_x_discrete(labels=trait_type_names, name=NULL) +
    guides(fill=F) + theme(axis.text.x = element_text(angle = 35, hjust=1))
}
figure2 = function(output_format = 'png') {
  p2a = generate_heritability_table(8.5)
  p2b = heritability_correlations(remove_questionnaire = F)
  p2c = passing_traits()
  output_type(output_format, paste0('figure2.', output_format), height=2.5, width=8.5)
  print(ggarrange(p2a, p2b, p2c, nrow=1, ncol = 3, labels="auto", widths = c(1.5, 1, 1)))
  dev.off()
}
figure2()
