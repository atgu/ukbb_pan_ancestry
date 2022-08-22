library(tidyverse)
library(RColorBrewer)
library(maptools)
library(cowplot)

setwd('/Users/alicia/daly_lab/ukbb_diverse_pops/pca')

pop_assign <- read_delim(gzfile('globalref_ukbb_pca_pops_rf_50.txt.gz'), delim='\t') %>%
  select(s, pop) %>%
  mutate(pop=case_when(pop=='oth' ~ 'Other',
                       TRUE ~ pop))
pop_assign$s <- as.character(pop_assign$s)
ethnicity <- read_tsv('../data/ukb31063.ethnicity_birth.txt')
ethnicity_ancestry <- pop_assign %>%
  left_join(ethnicity, by=c('s'='userId')) %>%
  mutate(Ethnicity=case_when(ethnicity=='-3' ~ 'Prefer not to answer',
                             ethnicity=='-1' ~ 'Do not know',
                             ethnicity=='1' ~ 'White',
                             ethnicity=='1001' ~ 'British',
                             ethnicity=='1002' ~ 'Irish',
                             ethnicity=='1003' ~ 'Any other white background',
                             ethnicity=='2' ~ 'Mixed',
                             ethnicity=='2001' ~ 'White and Black Caribbean',
                             ethnicity=='2002' ~ 'White and Black African',
                             ethnicity=='2003' ~ 'White and Asian',
                             ethnicity=='2004' ~ 'Any other mixed background',
                             ethnicity=='3' ~ 'Asian or Asian British',
                             ethnicity=='3001' ~ 'Indian',
                             ethnicity=='3002' ~ 'Pakistani',
                             ethnicity=='3003' ~ 'Bangladeshi',
                             ethnicity=='4' ~ 'Black or Black British',
                             ethnicity=='4001' ~ 'Caribbean',
                             ethnicity=='4002' ~ 'African',
                             ethnicity=='4003' ~ 'Any other Black background',
                             ethnicity=='5' ~ 'Chinese',
                             ethnicity=='6' ~ 'Other ethnic group',
                             TRUE ~ 'NA'
    
  ))
table(ethnicity_ancestry$ethnicity, ethnicity_ancestry$pop)
table(ethnicity_ancestry$continent, ethnicity_ancestry$pop)
conts <- unique(pop_assign$pop)

immigrant_country_ancestry <- function(cont) {
  a <- as.data.frame(table(ethnicity_ancestry$pop, ethnicity_ancestry$country)) %>% 
    subset(Var1==cont) %>%
    arrange(desc(Freq)) %>%
    filter(Freq>0)
  print(a)
  return(a)
}

anc_country <- map_dfr(conts, immigrant_country_ancestry)
write.table(anc_country, 'country_immigration_ancestry.txt', quote=F, row.names=F, sep='\t')


# Plot global PCs ---------------------------------------------------------

ref <- read.table(gzfile('globalref_ukbb_scores.txt.bgz'), header=T)
ref_info <- read.csv('/Users/alicia/Dropbox (Partners HealthCare)/martin_lab/projects/hgdp_tgp/tgp_hgdp.csv', header=T)
ref <- ref %>% 
  left_join(ref_info, by=c('s'='Sample.ID')) %>%
  mutate(Population=Genetic.region) %>%
  select(c(s, starts_with('PC'), Population))
  
ukbb <- read.table(gzfile('ukbb_globalref_scores.txt.bgz'), header=T) %>%
  mutate(Population='UKBB')
ukbb$s <- as.character(ukbb$s)
ref_ukbb <- bind_rows(ukbb, ref)

ukb_pop_colors = c('AFR' = '#941494', 'AMR' = '#ED1E24', 'CSA' = '#FF9912', 'EAS' = '#108C44', 'EUR' = '#6AA5CD', 'MID' = '#EEA9B8', 'OCE' = '#a6761d', 'UKBB' = 'black', 'Other' = '#ABB9B9')
ref_ukbb$Population <- factor(ref_ukbb$Population, levels=as.character(names(ukb_pop_colors)))

plot_global_pca <- function(pcs, pop_color, pop_shape=NA, first_pc='PC1', second_pc='PC2', legend_name='Population') {
  if(is.na(pop_shape)) {
    pca_pop <- ggplot(pcs, aes_string(x=first_pc, y=second_pc, color=legend_name)) + geom_point(alpha=0.8)
  } else {
    pca_pop <- ggplot(pcs, aes_string(x=first_pc, y=second_pc, color=legend_name, shape=legend_name)) + geom_point(alpha=0.8)
  }
  pca_pop <- pca_pop +
    scale_color_manual(values=pop_color, name=legend_name) +
    scale_shape_manual(values=pop_shape, name=legend_name) +
    theme_classic() +
    theme(text = element_text(size=14),
          axis.text = element_text(color='black'),
          legend.text = element_text(size=10))
  
  
  x_lim = ggplot_build(pca_pop)$layout$panel_scales_x[[1]]$range$range
  y_lim = ggplot_build(pca_pop)$layout$panel_scales_y[[1]]$range$range
  
  return(list(pca_pop, x_lim, y_lim))
}

global_pcs_1_2 <- plot_global_pca(ref_ukbb, ukb_pop_colors)
global_pcs_3_4 <- plot_global_pca(ref_ukbb, ukb_pop_colors, first_pc = 'PC3', second_pc = 'PC4')  
global_pcs_5_6 <- plot_global_pca(ref_ukbb, ukb_pop_colors, first_pc = 'PC5', second_pc = 'PC6') 

plot_pca_density <- function(pcs, x_lim, y_lim, first_pc='PC1', second_pc = 'PC2') {
  dens_pca <- ggplot(pcs, aes_string(x=first_pc, y=second_pc)) +
    geom_hex(bins=50) +
    scale_fill_gradientn(trans = "log", breaks=c(1,20,400,8000,163000), name='Count',
                         colours = rev(brewer.pal(5,'Spectral'))) +
    theme_classic() +
    lims(x=x_lim, y=y_lim) +
    theme(text = element_text(size=14),
          axis.text = element_text(color='black'),
          legend.text = element_text(size=10))
  return(dens_pca)
}

ukb_dens_pcs_1_2 <- plot_pca_density(ukbb, global_pcs_1_2[[2]], global_pcs_1_2[[3]])
ukb_dens_pcs_3_4 <- plot_pca_density(ukbb, global_pcs_3_4[[2]], global_pcs_3_4[[3]], first_pc='PC3', second_pc='PC4')
ukb_dens_pcs_5_6 <- plot_pca_density(ukbb, global_pcs_5_6[[2]], global_pcs_5_6[[3]], first_pc='PC5', second_pc='PC6')

ukbb_assign <- ukbb %>%
  left_join(pop_assign) %>%
  mutate(Population=pop) %>%
  select(-pop) %>%
  filter(!is.na(Population))

global_pcs_1_2_assign <- plot_global_pca(ukbb_assign, ukb_pop_colors)
global_pcs_3_4_assign <- plot_global_pca(ukbb_assign, ukb_pop_colors, first_pc = 'PC3', second_pc = 'PC4')
global_pcs_5_6_assign <- plot_global_pca(ukbb_assign, ukb_pop_colors, first_pc = 'PC5', second_pc = 'PC6')

global_pcs <- plot_grid(global_pcs_1_2[[1]], global_pcs_3_4[[1]], global_pcs_5_6[[1]], 
                        ukb_dens_pcs_1_2, ukb_dens_pcs_3_4, ukb_dens_pcs_5_6,
                        global_pcs_1_2_assign[[1]], global_pcs_3_4_assign[[1]], global_pcs_5_6_assign[[1]], 
                        labels=LETTERS[1:9], nrow=3)

ggsave('ukbb_ref_agg_dens_rf_pca.png', global_pcs, width=12.5, height=10)


# Plot PCA x self-reported ethnicity info ---------------------------------

blues <- brewer.pal(5, 'Blues')[2:5] #1-1003
reds <- brewer.pal(6, 'Reds')[2:6] #2-2004
oranges <- brewer.pal(6, 'Oranges')[2:6] #3-3004
purples <- brewer.pal(5, 'Purples')[2:5] #4-4003
greys <- brewer.pal(4, 'Greys')[2:4] #-3, -1, 6
ethnicity_colors = c(greys[1:2], blues, reds, oranges, purples, '#108C44', greys[3])
names(ethnicity_colors) <- c('Prefer not to answer', 'Do not know', 'White', 'British', 'Irish', 'Any other white background',
                             'Mixed', 'White and Black Caribbean', 'White and Black African', 'White and Asian', 'Any other mixed background',
                             'Asian or Asian British',
                             'Indian', 'Pakistani', 'Bangladeshi', 'Any other Asian background',
                             'Black or Black British', 'Caribbean', 'African', 'Any other Black background',
                             'Chinese',
                             'Other ethnic group')
ethnicity_shapes <- c(rep(4, 2), rep(15, 4), rep(16, 5), rep(17, 5), rep(18, 4), 3, 4)
names(ethnicity_shapes) <- names(ethnicity_colors)
ethnicity_ancestry$s <- as.character(ethnicity_ancestry$s)
ukbb_ethnicity <- ukbb %>%
  left_join(ethnicity_ancestry %>% filter(Ethnicity!='NA'))
ukbb_ethnicity$Ethnicity <- factor(ukbb_ethnicity$Ethnicity, levels = names(ethnicity_colors))

global_pcs_1_2_eth <- plot_global_pca(ukbb_ethnicity, ethnicity_colors, ethnicity_shapes, legend_name='Ethnicity')
global_pcs_3_4_eth <- plot_global_pca(ukbb_ethnicity, ethnicity_colors, ethnicity_shapes, first_pc = 'PC3', second_pc = 'PC4', legend_name='Ethnicity')
global_pcs_5_6_eth <- plot_global_pca(ukbb_ethnicity, ethnicity_colors, ethnicity_shapes, first_pc = 'PC5', second_pc = 'PC6', legend_name='Ethnicity')

legend_b <- get_legend(global_pcs_1_2_eth[[1]] + theme(legend.position='bottom'))
anc_eth <- plot_grid(global_pcs_1_2_eth[[1]] + guides(color=F, shape=F), 
                     global_pcs_3_4_eth[[1]] + guides(color=F, shape=F), 
                     global_pcs_5_6_eth[[1]] + guides(color=F, shape=F),
          labels=LETTERS[1:3], nrow=1)
anc_eth_legend <- plot_grid(anc_eth, legend_b, ncol=1, rel_heights=c(1, .4))
ggsave('ukbb_pca_eth.png', anc_eth_legend, width=12, height=7)
# 
# 
# ukbb$s <- as.character(ukbb$s)
# rf <- read.table('ukbb_pca_pops_rf.txt.gz', header=T) %>% select(c('s', 'pop'))
# rf$s <- as.character(rf$s)
# in_gwas <- read.table('../ukb31063.gwas_samples.both_sexes.txt', header=T)
# pigment <- read.table('../skin_color_tanning.txt.gz', header=T, sep='\t')
# pigment$s <- as.character(pigment$s)
# ukbb2 <- ukbb %>%
#   mutate(in_gwas=ifelse(s %in% in_gwas$s, TRUE, FALSE)) %>%
#   left_join(rf, by='s') %>%
#   left_join(pigment, by='s')
# 
# tgp_pops <- read.table('../integrated_call_samples_v3.20130502.ALL.panel', header=T)
# tgp <- merge(tgp, tgp_pops, by.x='s', by.y='sample', all=T)
# 
# ukbb_tgp <- ukbb %>%
#   bind_rows(tgp) %>%
#   mutate(all_pop=ifelse(is.na(super_pop), 'UKBB', as.character(super_pop)))
# 
# brewer_vec <- brewer.pal(7, 'Set1')
# brewer_vec <- c(brewer_vec, 'black', 'black', brewer_vec[4])
# names(brewer_vec) <- c('EUR', 'EAS', 'AMR', 'SAS', 'AFR', 'MID', 'OCE', 'UKBB', 'oth', 'CSA')


# Plot location maps ------------------------------------------------------

data(wrld_simpl)
world <- fortify(wrld_simpl)

# Read the latitude/longitdue/plotting data for reference populations
pop_pos <- read.csv('/Users/alicia/Dropbox (Partners HealthCare)/daly_lab/UKBB-Diverse-Pops/data/pop_plot_info.csv', header=T) %>%
  filter(Population !='ASW')

plot_cont_map <- function(cont_name, lon_lim, lat_lim, rand_col=FALSE) {
  pop_pos_plot <- subset(pop_pos, Continent == cont_name)
  pop_pos_plot$Population <- factor(pop_pos_plot$Population, levels=as.character(pop_pos_plot$Population))
  if(rand_col) {
    color_vec <- colorRampPalette(brewer.pal(4, 'Spectral'))(length(pop_pos_plot$Population))
  }
  else {
    color_vec <- as.character(pop_pos_plot$Color)
  }
  shape_vec <- rep_len(c(21:25), length.out = length(color_vec))
  names(color_vec) <- pop_pos_plot$Population
  names(shape_vec) <- pop_pos_plot$Population
  
  # plot the map of Africa with data points labeled
  p_map <- ggplot() +
    geom_polygon(data = world, aes(long, lat, group=group), fill='lightyellow', color='lightgrey') +
    geom_point(data = pop_pos_plot, aes(Longitude, Latitude, color=Population, fill=Population, shape=Population), size=3) +
    coord_fixed(xlim = lon_lim, ylim = lat_lim) +
    labs(x='Longitude', y='Latitude') +
    theme_classic() +
    scale_fill_manual(name = "Population",
                      values = color_vec) +
    scale_color_manual(name = "Population",
                       values = color_vec) +
    scale_shape_manual(name = "Population",
                       values = shape_vec) +
    theme(panel.background = element_rect(fill = "lightblue"),
          plot.background = element_rect(fill = "transparent", color = NA),
          #legend.position='bottom',
          text = element_text(size=14),
          axis.text = element_text(color='black'),
          legend.text = element_text(size=10))
  return(list(p_map, color_vec, shape_vec))
}

afr <- plot_cont_map('AFR', c(-20,50), c(-35,35))
amr <- plot_cont_map('AMR', c(-140,-35), c(-50,65), rand_col=TRUE)
csa <- plot_cont_map('CSA', c(60,95), c(5,45), rand_col=TRUE)
eas <- plot_cont_map('EAS', c(78,148), c(0,70), rand_col=TRUE)
eur <- plot_cont_map('EUR', c(-25,40), c(34,71), rand_col=TRUE)
mid <- plot_cont_map('MID', c(0,60), c(10,50), rand_col=TRUE)

ggsave('afr_ref_map.pdf', afr[[1]], width=10, height=10)
ggsave('csa_ref_map.pdf', csa[[1]])
ggsave('eas_ref_map.pdf', eas[[1]])
ggsave('mid_ref_map.pdf', mid[[1]])
ggsave('amr_ref_map.pdf', amr[[1]])
ggsave('eur_ref_map.pdf', eur[[1]])
p2 <- p1 + guides(fill=F, color=F, shape=F)


# Load population PCA and covariate info ----------------------------------

# NOTE: change order here to correspond to order in color_vec
load_ref_pcs <- function(ref_pcs, ref_fam, pop_color) {
  ref_pcs <- read_delim(gzfile(ref_pcs), delim='\t')
  fam <- read.table(ref_fam, col.names=c('pop', 's', 'dad', 'mom', 'sex', 'pheno'))  %>%
    select(pop, s, sex)
  ref_data <- merge(ref_pcs, fam, by='s')
  ref_data$pop <- factor(ref_data$pop, levels = names(pop_color))
  return(ref_data)
}

ref_afr <- load_ref_pcs('AFR_HGDP_1kG_AGVP_maf005_geno05_unrel_ukbb_scores.txt.bgz', 'AFR_HGDP_1kG_AGVP_maf005_geno05_unrel.fam', afr[[2]])
ref_amr <- load_ref_pcs('AMR_HGDP_1kG_maf005_geno05_unrel_ukbb_scores.txt.bgz', 'AMR_HGDP_1kG_maf005_geno05_unrel.fam', amr[[2]])
ref_csa <- load_ref_pcs('CSA_HGDP_1kG_maf005_geno05_unrel_ukbb_scores.txt.bgz', 'CSA_HGDP_1kG_maf005_geno05_unrel.fam', csa[[2]])
ref_eas <- load_ref_pcs('EAS_HGDP_1kG_maf005_geno05_unrel_ukbb_scores.txt.bgz', 'EAS_HGDP_1kG_maf005_geno05_unrel.fam', eas[[2]])
ref_eur <- load_ref_pcs('EUR_HGDP_1kG_maf005_geno05_unrel_ukbb_scores.txt.bgz', 'EUR_HGDP_1kG_maf005_geno05_unrel.fam', eur[[2]])
ref_mid <- load_ref_pcs('MID_HGDP_1kG_maf005_geno05_unrel_ukbb_scores.txt.bgz', 'MID_HGDP_1kG_maf005_geno05_unrel.fam', mid[[2]])

# Load UKB population data ------------------------------------------------

load_ukb <- function(cont_name, filename) {
  ukb_pop <- read_delim(gzfile(filename), delim='\t') %>%
    left_join(pop_assign) %>%
    filter(pop==cont_name) %>%
    left_join(ethnicity, by=c('s'='userId')) #####
}
ukb_afr <- load_ukb('AFR', 'ukbb_AFR_HGDP_1kG_AGVP_maf005_geno05_unrel_scores.txt.bgz')
ukb_csa <- load_ukb('CSA', 'ukbb_CSA_HGDP_1kG_maf005_geno05_unrel_scores.txt.bgz')
ukb_eas <- load_ukb('EAS', 'ukbb_EAS_HGDP_1kG_maf005_geno05_unrel_scores.txt.bgz')
ukb_eur <- load_ukb('EUR', 'ukbb_EUR_HGDP_1kG_maf005_geno05_unrel_scores.txt.bgz')
ukb_mid <- load_ukb('MID', 'ukbb_MID_HGDP_1kG_maf005_geno05_unrel_scores.txt.bgz')
ukb_amr <- load_ukb('AMR', 'ukbb_AMR_HGDP_1kG_maf005_geno05_unrel_scores.txt.bgz')

#ukb_afr %>% dplyr::count(country) %>% arrange(desc(n)) %>% head(11)


# Plot PCA ----------------------------------------------------------------

p_afr <- ggplot(afr, aes(x=PC1, y=PC2, color=pop)) +
  geom_point(data=ukb_afr, color='grey') +
  geom_point() +
  scale_color_manual(values=color_vec, name='Population') +
  theme_classic() +
  theme(text = element_text(size=16))

ggsave('afr_cont_projection.pdf', p_afr, width=8, height=6)

plot_pca_ref_ukb <- function(ref_pop, ukb_pop, pop_color, pop_shape, first_pc='PC1', second_pc='PC2') {
  pca_pop <- ggplot(ref_pop, aes_string(x=first_pc, y=second_pc, color='pop', fill='pop', shape='pop')) +
    geom_point(data=ukb_pop, color='grey', fill='grey', shape=21) +
    geom_point() +
    scale_color_manual(values=pop_color, name='Population') +
    scale_fill_manual(values=pop_color, name='Population') +
    scale_shape_manual(values=pop_shape, name='Population') +
    guides(color=F, fill=F, shape=F) +
    theme_classic() +
    theme(text = element_text(size=12))
  
  x_lim <- ggplot_build(pca_pop)$layout$panel_scales_x[[1]]$range$range
  y_lim <- ggplot_build(pca_pop)$layout$panel_scales_y[[1]]$range$range
  pca_density <- ggplot(ukb_pop, aes_string(x=first_pc, y=second_pc)) +
    geom_hex(bins=50) +
    scale_fill_gradientn(trans='sqrt', name='Count',
                         colours = rev(brewer.pal(5,'Spectral'))) +
    lims(x=x_lim, y=y_lim) +
    theme_classic() +
    theme(text = element_text(size=12))
  
  return(list(pca_pop, pca_density))
}

save_pca_plot <- function(pop, pop_name, ref_pop, ukb_pop, base_height, base_width) {
  p_pop_1_2 <- plot_pca_ref_ukb(ref_pop, ukb_pop, pop[[2]], pop[[3]], 'PC1', 'PC2')
  p_pop_3_4 <- plot_pca_ref_ukb(ref_pop, ukb_pop, pop[[2]], pop[[3]], 'PC3', 'PC4')
  p_pop_5_6 <- plot_pca_ref_ukb(ref_pop, ukb_pop, pop[[2]], pop[[3]], 'PC5', 'PC6')
  my_plot_1_2=plot_grid(p_pop_1_2[[1]], p_pop_1_2[[2]], rel_widths=c(1, 1.15))
  my_plot_3_4=plot_grid(p_pop_3_4[[1]], p_pop_3_4[[2]], rel_widths=c(1, 1.15))
  my_plot_5_6=plot_grid(p_pop_5_6[[1]], p_pop_5_6[[2]], rel_widths=c(1, 1.15))
  my_plot = plot_grid(pop[[1]], my_plot_1_2, my_plot_3_4, my_plot_5_6, ncol=1, labels=c('A', 'B', 'C', 'D'), rel_heights=c(1.5, 1, 1, 1))
  save_plot(paste0(pop_name, '_cont_projection_1-6.png'), my_plot, base_height = base_height, base_width = base_width)
  # save_plot(filename=paste0(pop_name, '_cont_projection_1_2.png'), plot=my_plot, base_height = 5, base_width=10)
  # save_plot(paste0(pop_name, '_cont_projection_3_4.png'), plot_grid(p_pop_3_4[[1]], p_pop_3_4[[2]], rel_widths=c(1, 1.15)), base_height = 5, base_width=10)
  # save_plot(paste0(pop_name, '_cont_projection_5_6.png'), plot_grid(p_pop_5_6[[1]], p_pop_5_6[[2]], rel_widths=c(1, 1.15)), base_height = 5, base_width=10)
}

save_pca_plot(afr, 'afr', ref_afr, ukb_afr, base_height=12, base_width=8)
save_pca_plot(amr, 'amr', ref_amr, ukb_amr, base_height=12, base_width=8)
save_pca_plot(csa, 'csa', ref_csa, ukb_csa, base_height=12, base_width=8)
save_pca_plot(eas, 'eas', ref_eas, ukb_eas, base_height=12, base_width=8)
save_pca_plot(mid, 'mid', ref_mid, ukb_mid, base_height=12, base_width=8)
save_pca_plot(eur, 'eur', ref_eur, ukb_eur, base_height=12, base_width=8)

# p_afr_1_2 <- plot_pca_ref_ukb(ref_afr, ukb_afr, afr[[2]], afr[[3]], 'PC1', 'PC2')
# p_afr_3_4 <- plot_pca_ref_ukb(ref_afr, ukb_afr, afr[[2]], afr[[3]], 'PC3', 'PC4')
# p_afr_5_6 <- plot_pca_ref_ukb(ref_afr, ukb_afr, afr[[2]], afr[[3]], 'PC5', 'PC6')
# my_plot <- plot_grid(p_afr_1_2[[1]], p_afr_1_2[[2]], rel_widths=c(1, 1.15))
# save_plot('afr_cont_projection_1_2.png', my_plot, base_height = 5, base_width=10)
# save_plot('afr_cont_projection_3_4.png', plot_grid(p_afr_3_4[[1]], p_afr_3_4[[2]], rel_widths=c(1, 1.15)), base_height = 5, base_width=10)
# save_plot('afr_cont_projection_5_6.png', plot_grid(p_afr_5_6[[1]], p_afr_5_6[[2]], rel_widths=c(1, 1.15)), base_height = 5, base_width=10)


# Within pop PCA (no ref) -------------------------------------------------

setwd('/Users/alicia/daly_lab/ukbb_diverse_pops/pca/ukb_within_continent')

read_pca <- function(pop_name, rel_unrel) {
  if(rel_unrel == 'rel') {
    pca <- read.table(gzfile(paste0(pop_name, '_rel_scores.txt.bgz')), header=T) %>%
      mutate(pop=pop_name, rel=rel_unrel)
  } else {
    pca <- read.table(gzfile(paste0(pop_name, '_scores.txt.bgz')), header=T) %>%
      mutate(pop=pop_name, rel=rel_unrel)
  }
  return(pca)
}

pops <- c('AFR', 'AMR', 'CSA', 'MID', 'EAS', 'EUR')

age_sex <- read.table(gzfile('uk_round2_allSamples_phenos_phesant.6148_5.tsv.gz'), header=T, sep='\t') %>%
  select(userId, age, sex)

bind_rels <- function(pop) {
  pop_rel <- read_pca(pop, 'rel')
  pop_unrel <- read_pca(pop, 'unrel')
  pop_bind <- pop_rel %>% bind_rows(pop_unrel)
}

afr <- bind_rels('AFR')
amr <- bind_rels('AMR')
csa <- bind_rels('CSA')
mid <- bind_rels('MID')
eas <- bind_rels('EAS')
eur <- bind_rels('EUR')

bind_pops <- afr %>%
  bind_rows(amr) %>%
  bind_rows(csa) %>%
  bind_rows(mid) %>%
  bind_rows(eas) %>%
  bind_rows(eur) %>%
  left_join(age_sex, by=c('s'='userId')) %>%
  mutate(age2 = age^2, age_sex = age*sex, age2_sex = age^2 * sex)

write.table(bind_pops, 'within_pop_pc_covs.txt', quote=F, row.names=F, sep='\t')

plot_pca_density <- function(dataset, first_pc, second_pc) {
  pc_biplot <- ggplot(dataset, aes_string(x=first_pc, y=second_pc)) +
    geom_hex(bins=50) +
    scale_fill_gradientn(trans = "log", breaks=c(1,20,400,8000,163000), name='Count',
                         colours = rev(brewer.pal(5,'Spectral'))) +
    theme_classic() +
    theme(text = element_text(size=16))
  return(pc_biplot)
}

pop_ellipse <- function(df, num_ellipses) {
  # get mean and SD of each PC among each pop
  pc_nams <- paste("PC",1:10,sep="")
  mean_pcs <- colMeans(df[,pc_nams])
  sd_pcs <- apply(df[,pc_nams],2,sd)
  # compute centroid distance for each individual
  centroid_dist <- rep(0,nrow(df))
  for(i in 1:num_ellipses) {
    centroid_dist <- centroid_dist + (df[,pc_nams[i]]-mean_pcs[i])^2/(sd_pcs[i]^2)
  }
  pop_dist <- df %>%
    mutate(centroid_dist=centroid_dist)
  return(pop_dist)
}

pop_centroid <- function(ind_dist, cutpoint0, cutpoint1) {
  pop_cut <- subset(ind_dist, centroid_dist < cutpoint0)
  p_centroid <- ggplot(pop_cut, aes(x=centroid_dist)) + 
    geom_histogram(bins=50) + 
    labs(title=paste0('Sample size: ', nrow(subset(ind_dist, centroid_dist < cutpoint0)), ' -> ', nrow(subset(ind_dist, centroid_dist < cutpoint1)))) +
    xlab('Centroid distance') + 
    ylab('Count') +
    geom_vline(xintercept=cutpoint1) +
    theme_bw() +
    theme(text = element_text(size=16))
  return(list(p=p_centroid, pop_cut=pop_cut))
}

save_filt_plots <- function(pop_name, pop_dist, cutpoint0, cutpoint1) {
  p_centroid0 = pop_centroid(pop_dist, cutpoint0, cutpoint1)
  ggsave(paste0(pop_name, '_within_pop_centroid_nofilt.pdf'), p_centroid0$p, height=7, width=7)
  p2 <- plot_pca_density(subset(pop_dist, centroid_dist < cutpoint0), 'PC1', 'PC2')
  ggsave(paste0(pop_name, '_within_pop_nofilt_pc1_2.png'), p2)
  p3 <- plot_pca_density(subset(pop_dist, centroid_dist < cutpoint0), 'PC3', 'PC4')
  ggsave(paste0(pop_name, '_within_pop_nofilt_pc3_4.png'), p3)
  p4 <- plot_pca_density(subset(pop_dist, centroid_dist < cutpoint0), 'PC5', 'PC6')
  ggsave(paste0(pop_name, '_within_pop_nofilt_pc5_6.png'), p4)
  p_centroid1 = pop_centroid(pop_dist, cutpoint1, cutpoint1)
  ggsave(paste0(pop_name, '_within_pop_centroid_filt.pdf'), p_centroid1$p, height=7, width=7)
  p6 <- plot_pca_density(subset(pop_dist, centroid_dist < cutpoint1), 'PC1', 'PC2')
  ggsave(paste0(pop_name, '_within_pop_filt_pc1_2.png'), p6)
  p7 <- plot_pca_density(subset(pop_dist, centroid_dist < cutpoint1), 'PC3', 'PC4')
  ggsave(paste0(pop_name, '_within_pop_filt_pc3_4.png'), p7)
  p8 <- plot_pca_density(subset(pop_dist, centroid_dist < cutpoint1), 'PC5', 'PC6')
  ggsave(paste0(pop_name, '_within_pop_filt_pc5_6.png'), p8)
  my_plot=plot_grid(p_centroid0$p, p2, p3, p4, p_centroid1$p, p6, p7, p8, nrow=2)
  save_plot(paste0(pop_name, '_within_pop.png'), my_plot, base_height=10, base_width = 18)
  return(p_centroid1$pop_cut)
}

csa_cut <- save_filt_plots('csa', csa_dist <- pop_ellipse(csa, 3), 1000, 3) #3, 3
afr_cut <- save_filt_plots('afr', afr_dist <- pop_ellipse(afr, 3), 1000, 2)
eas_cut <- save_filt_plots('eas', eas_dist <- pop_ellipse(eas, 3), 1000, 7.5)
amr_cut <- save_filt_plots('amr', amr_dist <- pop_ellipse(amr, 3), 1000, 4.8)
mid_cut <- save_filt_plots('mid', mid_dist <- pop_ellipse(mid, 5), 1000, 15)
eur_cut <- save_filt_plots('eur', eur_dist <- pop_ellipse(eur, 5), 1000, 10)


pop_cuts <- csa_cut %>%
  bind_rows(afr_cut, eas_cut, amr_cut, mid_cut, eur_cut) %>%
  select(s, pop)

write.table(pop_cuts, 'ukb_diverse_pops_pruned.tsv', row.names=F, sep='\t', quote=F)

