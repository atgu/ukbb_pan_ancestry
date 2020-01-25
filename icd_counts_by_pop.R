source('~/gnomad_lof/R/constants.R')
library(tidyverse)
library(tidylog)
library(plotly)
library(scales)
library(tictoc)

data = read_delim('data/phenos_full.tsv', delim = '\t',
                  col_types=cols(pheno=col_character(), coding=col_character())) %>%
  filter(!(pop %in% c('OCE', 'oth')) & !is.na(pheno)) %>%
  mutate(pop=fct_recode(tolower(pop), sas = 'csa', mde = 'mid'))

data %>%
  filter(n_cases_by_pop > 200) %>%
  count(pop)

data %>%
  group_by(pop) %>%
  mutate(n_pheno=min_rank(-n_cases_by_pop)) %>%
  filter(n_cases_by_pop > 0) %>%
  ggplot + aes(y = n_pheno, x = n_cases_by_pop, color = pop) + 
  geom_line(lwd=1) + theme_bw() +
  ylab('Number of phenotypes') +
  # geom_vline(xintercept = 500, linetype='dashed') +
  scale_x_log10(label=comma, name='>= Number of cases') +
  annotation_logticks(sides='b') + 
  scale_color_manual(values = pop_colors, labels=pop_names, name='Population') -> p

ggplotly(p)
# scale_y_log10(label=comma_format(accuracy=1))

summary_data = data %>%
  filter(n_cases_by_pop > 1000 & coding != 'raw') %>%
  group_by(pheno, coding, meaning) %>%
  summarize(max_mean = max(stats.mean),
            max_pop = pop[which.max(stats.mean)],
            min_mean = min(stats.mean),
            min_pop = pop[which.min(stats.mean)],
            max_stdev = max(stats.stdev),
            max_pop_stdev = pop[which.max(stats.stdev)],
            min_stdev = min(stats.stdev),
            min_pop_stdev = pop[which.min(stats.stdev)],
            n=n()) %>% ungroup %>%
  filter(n > 1)

summary_data %>%
  arrange(desc((max_mean - min_mean) / (max_stdev))) %>%
  head(20)

summary_data %>%
  arrange(desc((max_stdev - min_stdev) / (min_stdev))) %>%
  head(20)

# data %>%
#   group_by(icd_code, truncated, meaning) %>%
#   summarize_if(is.numeric, sum) %>% ungroup %>%
#   filter(n_cases_all > 100) %>%
#   count(truncated)

pairwise_data = read_delim(gzfile('data/pairwise_correlations.txt.bgz'), delim = '\t',
                           col_types = cols(i_data.pheno=col_character(), j_data.pheno=col_character(),
                                            i_data.coding=col_character(), j_data.coding=col_character()))

full_matrix = pairwise_data %>%
  select(i, j, entry) %>%
  pivot_wider(names_from = j, values_from = entry)

pairwise_data %>%
  filter(i_data.data_type == 'prescriptions' & j_data.data_type == 'prescriptions' & !is.nan(entry)) %>%
  select(i, j, entry) %>%
  pivot_wider(names_from = j, values_from = entry) %>%
  select(-i) %>%
  as.matrix -> test_matrix

size = 200
d = as.dist(1 - test_matrix[1:size, 1:size] ^ 2)
hc1 <- hclust(d, method = "complete")
plot(hc1, cex = 0.5, hang = -1)

size = 100
tic(); corrplot(test_matrix[1:size, 1:size], method='circle', order='hclust', tl.col='black', addgrid.col = NA); toc()


