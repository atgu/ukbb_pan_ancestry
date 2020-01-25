library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cowplot)

setwd('~/daly_lab/skin_pigmentation/pca/')

# Plot UKBB with 1000 Genomes ---------------------------------------------

tgp <- read.table('tgp_pca_scores.txt.bgz', header=T)
ukbb <- read.table('ukbb_pca_scores.txt.bgz', header=T)
ukbb$s <- as.character(ukbb$s)
rf <- read.table('ukbb_pca_pops_rf.txt.gz', header=T) %>% select(c('s', 'pop'))
rf$s <- as.character(rf$s)
in_gwas <- read.table('../ukb31063.gwas_samples.both_sexes.txt', header=T)
pigment <- read.table('../skin_color_tanning.txt.gz', header=T, sep='\t')
pigment$s <- as.character(pigment$s)
ukbb2 <- ukbb %>%
  mutate(in_gwas=ifelse(s %in% in_gwas$s, TRUE, FALSE)) %>%
  left_join(rf, by='s') %>%
  left_join(pigment, by='s')

tgp_pops <- read.table('../integrated_call_samples_v3.20130502.ALL.panel', header=T)
tgp <- merge(tgp, tgp_pops, by.x='s', by.y='sample', all=T)

ukbb_tgp <- ukbb %>%
  bind_rows(tgp) %>%
  mutate(all_pop=ifelse(is.na(super_pop), 'UKBB', as.character(super_pop)))

brewer_vec <- brewer.pal(7, 'Set1')
brewer_vec <- c(brewer_vec, 'black', 'black', brewer_vec[4])
names(brewer_vec) <- c('EUR', 'EAS', 'AMR', 'SAS', 'AFR', 'MID', 'OCE', 'UKBB', 'oth', 'CSA')
plot_pca <- function(dataset, first_pc, second_pc, color_pop) {
  pc_biplot <- ggplot(dataset, aes_string(x=first_pc, y=second_pc, color=color_pop)) +
    geom_point(alpha=0.2) +
    guides(color=guide_legend(override.aes = list(alpha=1))) +
    theme_classic() +
    scale_color_manual(values=brewer_vec, name='Population')
  return(pc_biplot)
}

pc1_2 <- plot_pca(ukbb_tgp, 'PC1', 'PC2', 'all_pop')
pc3_4 <- plot_pca(ukbb_tgp, 'PC3', 'PC4', 'all_pop')
pc5_6 <- plot_pca(ukbb_tgp, 'PC5', 'PC6', 'all_pop')
pc7_8 <- plot_pca(ukbb_tgp, 'PC7', 'PC8', 'all_pop')

pc_agg <- plot_grid(pc1_2, pc3_4, pc5_6, labels=c('A', 'B', 'C'), ncol=3)
save_plot('pc1-5_2.png', pc_agg, base_height=3.5, base_width=12)
ggsave('pc7_8.png', pc7_8)

plot_pca_density <- function(dataset, first_pc, second_pc) {
  pc_biplot <- ggplot(dataset, aes_string(x=first_pc, y=second_pc)) +
    geom_hex(bins=50) +
    scale_fill_gradientn(trans = "log", breaks=c(1,20,400,8000,163000), name='Count',
                         colours = rev(brewer.pal(5,'Spectral'))) +
    theme_classic()
  return(pc_biplot)
}

pc1_2_dens <- plot_pca_density(subset(ukbb_tgp, all_pop=='UKBB'), 'PC1', 'PC2')
pc3_4_dens <- plot_pca_density(subset(ukbb_tgp, all_pop=='UKBB'), 'PC3', 'PC4')
pc5_6_dens <- plot_pca_density(subset(ukbb_tgp, all_pop=='UKBB'), 'PC5', 'PC6')
pc_agg_dens <- plot_grid(pc1_2_dens, pc3_4_dens, pc5_6_dens, labels=c('A', 'B', 'C'), ncol=3)
save_plot('pc1-5_dens.png', pc_agg_dens, base_height=3.5, base_width=12)


# Plot UKBB random forest assigned pops (TGP) -----------------------------

rf <- read.table('ukbb_pca_pops_rf.txt.gz', header=T)

pc1_2 <- plot_pca(rf, 'PC1', 'PC2', 'pop')
pc3_4 <- plot_pca(rf, 'PC3', 'PC4', 'pop')
pc5_6 <- plot_pca(rf, 'PC5', 'PC6', 'pop')
pc7_8 <- plot_pca(rf, 'PC7', 'PC8', 'pop')

pc_agg <- plot_grid(pc1_2, pc3_4, pc5_6, labels=c('A', 'B', 'C'), ncol=3)
save_plot('pc1-5_rf.png', pc_agg, base_height=3.5, base_width=12)


# Plot UKBB with 1000 Genomes and HGDP ------------------------------------

setwd('~/daly_lab/ukbb_diverse_pops/pca/')
stringsAsFactors = FALSE

ref <- read.table('globalref_ukbb_scores.txt.bgz', header=T)
ukbb <- read.table('ukbb_globalref_scores.txt.bgz', header=T)
ref$s <- as.character(ref$s)
ukbb$s <- as.character(ukbb$s)
metadata <- read.csv('../data/tgp_hgdp.csv', header=T)
ref_data <- bind_rows(ukbb, ref) %>%
  left_join(metadata, by=c('s'='Sample.ID')) %>%
  mutate(pop=case_when(is.na(Genetic.region)~'UKBB',
                             TRUE~as.character(Genetic.region)))
ref_data <- ref_data %>%
  mutate(all_pop=ifelse(is.na(Genetic.region), 'UKBB', as.character(Genetic.region)))

pc1_2 <- plot_pca(ref_data, 'PC1', 'PC2', 'pop')
pc3_4 <- plot_pca(ref_data, 'PC3', 'PC4', 'pop')
pc5_6 <- plot_pca(ref_data, 'PC5', 'PC6', 'pop')
pc7_8 <- plot_pca(ref_data, 'PC7', 'PC8', 'pop')

pc_agg <- plot_grid(pc1_2, pc3_4, pc5_6, labels=c('A', 'B', 'C'), ncol=3)
save_plot('pc1-5_rf.png', pc_agg, base_height=3.5, base_width=12)

pc1_2_dens <- plot_pca_density(subset(ref_data, all_pop=='UKBB'), 'PC1', 'PC2')
pc3_4_dens <- plot_pca_density(subset(ref_data, all_pop=='UKBB'), 'PC3', 'PC4')
pc5_6_dens <- plot_pca_density(subset(ref_data, all_pop=='UKBB'), 'PC5', 'PC6')
pc_agg_dens <- plot_grid(pc1_2_dens, pc3_4_dens, pc5_6_dens, labels=c('A', 'B', 'C'), ncol=3)
save_plot('pc1-5_dens.png', pc_agg_dens, base_height=3.5, base_width=12)


# Compare reference panel RF ----------------------------------------------

rf_tgp <- read.table('/Users/alicia/daly_lab/skin_pigmentation/pca/ukbb_pca_pops_rf_50.txt.gz', header=T) %>% mutate(ref='tgp')
rf_hgdp <- read.table('/Users/alicia/daly_lab/ukbb_diverse_pops/pca/globalref_ukbb_pca_pops_rf_50.txt.gz', header=T) %>% mutate(ref='hgdp')
rf_tgp$pop <- factor(rf_tgp$pop, levels=c('EUR', 'SAS', 'AFR', 'EAS', 'AMR', 'oth'))
rf_hgdp$pop <- factor(rf_hgdp$pop, levels=c('EUR', 'CSA', 'AFR', 'EAS', 'AMR', 'oth', 'MID', 'OCE'))
rf <- merge(rf_tgp, rf_hgdp, by='s')
table(rf$pop.x, rf$pop.y)

ethnicity <- read.table('../data/ethnicity.tsv.bgz', header=T)
rf_ethnicity <- merge(rf_hgdp, ethnicity, by.x='s', by.y='userId')



# Compare # PCs used in RF ------------------------------------------------

rf_hgdp50 <- read.table('/Users/alicia/daly_lab/ukbb_diverse_pops/pca/globalref_ukbb_pca_pops_rf_50_20PCs.txt.gz', header=T) %>% mutate(ref='20PCs')
rf <- merge(rf_hgdp, rf_hgdp50, by='s')
table(rf$pop.x, rf$pop.y)

# Plot UKBB Africans with TGP + AGVP Africans -----------------------------

ukbb_afr <- read.table('ukbb_afr_pca_scores.txt.bgz', header=T)
ukbb_afr$pop <- 'UKBB'
ukbb_afr <- ukbb_afr %>% select(pop, s, PC1:PC20)

agvp_tgp <- read.table('afr_pca_scores.txt.bgz', header=T)
agvp_tgp_fam <- read.table('TGP_AGVP.postQC.autosomes.geno10.unrel.fam')
colnames(agvp_tgp_fam) <- c('pop', 's', 'pat_id', 'mat_id', 'sex', 'is_case')
agvp_tgp <- merge(agvp_tgp, agvp_tgp_fam, by='s') %>% select(pop, s, PC1:PC20)

ukbb_afr_comb <- rbind(ukbb_afr, agvp_tgp)
pop_info <- read.csv('/Users/alicia/daly_lab/GINGER/agvp/african_pop_summary3.csv', header=T)
color_vec <- c(as.character(pop_info$Color), '#bdbdbd', '#636363', 'black')
names(color_vec) <- c(as.character(pop_info$Population), 'ASW', 'ACB', 'UKBB')
shape_vec <- c(pop_info$Shape, 21, 21, 21)
names(shape_vec) <- c(as.character(pop_info$Population), 'ASW', 'ACB', 'UKBB')

ukbb_afr_comb$pop <- factor(ukbb_afr_comb$pop, levels=c(as.character(pop_info$Population), 'ACB', 'ASW', 'UKBB'))

afr_biplot <- function(first_pc, second_pc, legend_right=FALSE) {
  pca_biplot <- ggplot(ukbb_afr_comb, aes_string(x=first_pc, y=second_pc, color='pop', shape='pop', fill='pop')) + 
  geom_point(alpha=0.8) +
  guides(color=guide_legend(override.aes = list(alpha=1))) +
  scale_color_manual(name = "Population", values = color_vec) +
  scale_fill_manual(name = "Population", values = color_vec) +
  scale_shape_manual(name = "Population", values = shape_vec) +
  theme_classic()
  if(legend_right) {
    pca_biplot <- pca_biplot
  } else {
    pca_biplot <- pca_biplot + guides(fill=F, color=F, shape=F)
  }
  return(pca_biplot)
}

p_legend <- afr_biplot('PC1', 'PC2', legend_right=T)
pop_legend <- g_legend(p_legend)
pc1_2_afr <- p_legend + guides(fill=F, color=F, shape=F)
pc3_4_afr <- afr_biplot('PC3', 'PC4') 

afr_pca_density <- function(dataset, first_pc, second_pc, legend_right=T) {
  pc_biplot <- ggplot(dataset, aes_string(x=first_pc, y=second_pc)) +
    geom_hex(bins=50) +
    scale_fill_gradientn(trans = "log", breaks=c(1,5,20,75,200), name='Count',
                         colours = rev(brewer.pal(5,'Spectral')), limits=c(1,378)) +
    #guides(color=guide_legend(override.aes = list(alpha=1))) +
    theme_classic()
  if(legend_right) {
    pc_biplot <- pc_biplot
  } else {
    pc_biplot <- pc_biplot + guides(fill=F)
  }
  #scale_color_manual(values=brewer_vec, name='Population')
  return(pc_biplot)
}

p_legend_dens <- afr_pca_density(subset(ukbb_afr_comb, pop=='UKBB'), 'PC1', 'PC2', legend_right=T)
pop_legend_dens <- g_legend(p_legend_dens)
pc1_2_afr_dens <- afr_pca_density(subset(ukbb_afr_comb, pop=='UKBB'), 'PC1', 'PC2', legend_right=F)
pc3_4_afr_dens <- afr_pca_density(subset(ukbb_afr_comb, pop=='UKBB'), 'PC3', 'PC4', legend_right=F)

plot_grid(pc1_2_afr, pc3_4_afr, pc1_2_afr_dens, pc3_4_afr_dens, align='h')
p_top <- plot_grid(pc1_2_afr, pc3_4_afr)
p_bottom <- plot_grid(pc1_2_afr_dens, pc3_4_afr_dens)

pca_density_afr <- ggdraw() +
  draw_plot(p_top, x=0, y=0.5, width=0.7, height=0.5) +
  draw_plot(pop_legend, x=0.7, y=0.5, width=0.3, height=0.5) +
  draw_plot(p_bottom, x=0, y=0, width=0.7, height=0.5) +
  draw_plot(pop_legend_dens, x=0.7, y=0, width=0.3, height=0.5)

save_plot('ukbb_afr_pca_dens.pdf', pca_density_afr, base_height=7, base_width=10)
ggplotly(p1)


# AFR biplots -------------------------------------------------------------
p_1_2 <- afr_biplot('PC1', 'PC2', legend_right=T)
p_3_4 <- afr_biplot('PC3', 'PC4', legend_right=T)
p_5_6 <- afr_biplot('PC5', 'PC6', legend_right=T)
p_7_8 <- afr_biplot('PC7', 'PC8', legend_right=T)
p_9_10 <- afr_biplot('PC9', 'PC10', legend_right=T)
p_11_12 <- afr_biplot('PC11', 'PC12', legend_right=T)
p_13_14 <- afr_biplot('PC13', 'PC14', legend_right=T)
p_15_16 <- afr_biplot('PC15', 'PC16', legend_right=T)
p_17_18 <- afr_biplot('PC17', 'PC18', legend_right=T)
p_19_20 <- afr_biplot('PC19', 'PC20', legend_right=T)

pdf('ukbb_afr_pca_biplots.pdf', onefile=T)
p_1_2
p_3_4
p_5_6
p_7_8
p_9_10
p_11_12
p_13_14
p_15_16
p_17_18
p_19_20
dev.off()

# In mega-gwas ------------------------------------------------------------

in_gwas <- read.table('ukb31063.gwas_samples.both_sexes.txt', header=T)

rf2 <- rf %>%
  mutate(in_gwas=ifelse(s %in% in_gwas$s, TRUE, FALSE))

ggplot(ukbb_samp, aes(x=PC1, y=PC2)) +
  geom_hex() +
  theme_classic()


# Pigmentation ------------------------------------------------------------
setwd('~/daly_lab/skin_pigmentation/')

skin <- table(ukbb2$pop, ukbb2$X1717)
colnames(skin) <- c('Prefer no', 'Don\'t know', 'Very Fair', 'Fair', 'Light olive', 'Dark olive', 'Brown', 'Black')
skin <- skin[,3:8]
skin_norm <- t(apply(skin, 1, function(x)(x-min(x))/(max(x)-min(x))))
library(reshape2)
melted_cormat <- melt(t(skin_norm))
melted_cormat$Var2 <- factor(melted_cormat$Var2, levels=rev(c('EUR', 'oth', 'EAS', 'AMR', 'SAS', 'AFR')))

p_skin <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  #scale_fill_brewer(palette='PuBuGn', name='Proportion') +
  scale_fill_gradientn(colors=brewer.pal(6,'YlGnBu'), name='Proportion') +
  labs(x='Skin color', y='Population') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave('skin_pop.pdf', p_skin)

tan <- table(ukbb2$pop, ukbb2$X1727)
colnames(tan) <- c('Prefer no', 'Don\'t know', 'Very tan', 'Moderate tan', 'Mild/occasional tan', 'Only burn')
tan <- tan[,6:1]
tan <- tan[,1:4]
tan_norm <- t(apply(tan, 1, function(x)(x-min(x))/(max(x)-min(x))))
melted_cormat <- melt(t(tan_norm))
melted_cormat$Var2 <- factor(melted_cormat$Var2, levels=rev(c('EUR', 'oth', 'EAS', 'AMR', 'SAS', 'AFR')))

p_tan <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  #scale_fill_brewer(palette='PuBuGn', name='Proportion') +
  scale_fill_gradientn(colors=brewer.pal(6,'YlGnBu'), name='Proportion') +
  labs(x='Tanning', y='Population') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('tan_pop.pdf', p_tan)

hair <- table(ukbb2$pop, ukbb2$X1747)
colnames(hair) <- c('Prefer no', 'Don\'t know', 'Blond', 'Red', 'Light brown', 'Dark brown', 'Black', 'other')
hair <- hair[,3:8]
hair_norm <- t(apply(hair, 1, function(x)(x-min(x))/(max(x)-min(x))))
melted_cormat <- melt(t(hair_norm))
melted_cormat$Var2 <- factor(melted_cormat$Var2, levels=rev(c('EUR', 'oth', 'EAS', 'AMR', 'SAS', 'AFR')))

p_hair <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  #scale_fill_brewer(palette='PuBuGn', name='Proportion') +
  scale_fill_gradientn(colors=brewer.pal(6,'YlGnBu'), name='Proportion') +
  labs(x='Hair color', y='Population') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave('hair_pop.pdf', p_hair)

ggplot(sas)

smoothScatter(ukbb$PC1, ukbb$PC2)

natgeo <- read.csv('~/daly_lab/natgeo/NG_ancestry_proportions.csv', header=T)
natgeo2 <- natgeo %>%
  gather('pop', 'proportion', EUR:SAS)
ggplot(natgeo2, aes(proportion)) +
  geom_density() +
  facet_wrap(~pop)

ukbb_samp <- read.table('sampled_ukbb.txt.bgz', header=T)
ukbb_samp <- merge(ukbb_samp, ukbb, by='s')
ukbb2 <- read.table('../ukb31063.sample_qc.tsv.bgz', header=T)
tgp_proj <- merge(tgp_proj, tgp_pops, by.x='s', by.y='sample')
