source('~/gnomad_lof/R/constants.R')
library(tidyverse)
library(scales)

data = read_delim('data/phenos_icd.tsv', delim = '\t') %>%
  filter(!(pop %in% c('OCE', 'oth'))) %>%
  mutate(pop=fct_recode(tolower(pop), sas = 'csa', mde = 'mid'))

data %>%
  filter(n_cases_all > 100) %>%
  count(pop)

data %>%
  filter(truncated) %>%
  group_by(pop) %>%
  mutate(n_pheno=min_rank(-n_cases_all)) %>%
  filter(n_cases_all > 0) %>%
  ggplot + aes(y = n_pheno, x = n_cases_all, color = pop) + 
  geom_line(lwd=1) + theme_bw() +
  ylab('Number of phenotypes') +
  # geom_vline(xintercept = 500, linetype='dashed') +
  scale_x_log10(label=comma, name='>= Number of cases') +
  annotation_logticks(sides='b') + 
  scale_color_manual(values = pop_colors, labels=pop_names, name='Population')
  # scale_y_log10(label=comma_format(accuracy=1))

# data %>%
#   group_by(icd_code, truncated, meaning) %>%
#   summarize_if(is.numeric, sum) %>% ungroup %>%
#   filter(n_cases_all > 100) %>%
#   count(truncated)
  

data %>%
  filter(pop == 'EUR' & n_cases_all > 100) %>%
  count(truncated)
