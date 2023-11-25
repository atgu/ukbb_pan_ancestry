source('~/CHARR/R/constants.R')
library(magick)

phenotype_levels <- rev(c('Mean sphered cell volume', 'High light scatter reticulocyte count', 'Glycated haemoglobin (HbA1c)', 'Red blood cell (erythrocyte) distribution width', 'Red blood cell (erythrocyte) count'))
phenotype_labels <- rev(c('MCV', 'RET', 'HbA1c', 'RDW', 'RBC'))
OUTPUT <- '~/Dropbox (Partners HealthCare)/analysis/ukb_pan_ancestry/figures/figure5/'

pops = c('AFR', 'AMR', 'EAS', 'MID', 'EUR', 'CSA', 'Meta (raw)', 'Meta (hq)')
names(pops) = c('afr', 'amr', 'eas', 'mid', 'nfe', 'sas', 'meta_raw', 'meta_hq')
pop_colors['mid'] = "#33CC33"
pop_colors['meta_raw'] = "gray40"
pop_colors['meta_hq'] = "black"
pop_colors['meta'] = "black"

figure5a <- function(text_size = 3 , height=85/in2mm, width=100/in2mm, output_path=OUTPUT, save=T){
  g6pd_p <- read_delim('~/Dropbox (Partners HealthCare)/analysis/ukb_pan_ancestry/data/x_153764217_g6pd_all_associations.csv')
  g6pd_corr <- read_delim('~/Dropbox (Partners HealthCare)/analysis/ukb_pan_ancestry/data/rg_rp_5_pheno_g6pd.txt.bgz', delim='\t')
  sub_p <- g6pd_p %>%
    filter(Pvalue < 1e-3) %>% # phewas threshold:
    filter(pop == 'meta') %>%
    select(description, Pvalue) %>%
    distinct()
  rg_missing <- g6pd_corr %>%
    filter(!is.na(rg) & phenotype1 != phenotype2) %>%
    select(phenotype1 = phenotype2, description1=description2, phenotype2=phenotype1, description2=description1, rg=rg)
  g6pd_corr <- g6pd_corr %>%
    merge(., rg_missing, by = colnames(g6pd_corr)[1:4], all.x=T) %>%
    mutate(rg = if_else(is.na(rg.x), rg.y, rg.x)) %>%
    select(-rg.x, -rg.y) %>%
    mutate(r = if_else(phenotype1<phenotype2, rg, rp),
           label = if_else(phenotype1<phenotype2, 'rg', 'rp')) %>%
    merge(., sub_p, by.x='description1', by.y = 'description') %>%
    merge(., sub_p, by.x='description2', by.y = 'description') %>%
    mutate(description1 = factor(description1, levels = phenotype_levels, labels = phenotype_labels),
           description2 = factor(description2, levels = phenotype_levels, labels = phenotype_labels)) %>%
    mutate(description1 = if_else(Pvalue.x < 5e-8, paste0(description1, '*'), description1),
           description2 = if_else(Pvalue.y < 5e-8, paste0(description2, '*'), description2),) %>%
    mutate(description1 = factor(description1, levels = rev(c('MCV*', 'RET*', 'HbA1c*', 'RDW', 'RBC'))),
           description2 = factor(description2, levels = rev(c('MCV*', 'RET*', 'HbA1c*', 'RDW', 'RBC'))))
  p <- g6pd_corr %>%
    mutate(x = 'r[g]', y= 'r[p]') %>%
    filter(!is.na(rg)) %>%
    ggplot() +
    geom_tile(aes(x = description2, y = description1, fill = r)) +
    geom_text(aes(x = description2, y = description1, label = round(r,2)), size=text_size) +
    coord_equal() +
    scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
    labs(x =NULL, y = NULL) +
    themes +
    theme_classic() +
    theme( axis.text = element_text(size = text_size*2.5, color = 'black'),
           strip.text = element_text(size = text_size*2.5, face = 'bold'),
           legend.text = element_text(size = text_size*2),
           legend.title = element_text(size = text_size*2.5, face = 'bold'),
           ) +
    facet_grid(y~x, labeller = label_parsed)

  if(save){
    png(paste0(output_path, 'figure5a_g6pd_phenotype_correlation.png'), height = height, width = width, units = 'in', res = 300)
    print(p)
    dev.off()
  }
  return(p)
}

figure5c <- function(text_size = 6, height = 100/in2mm, width = 174/in2mm, output_path=OUTPUT, save=T, no_amr = F, no_meta_hq = F){
  result_data <- read_delim('/Users/wlu/Dropbox (Partners HealthCare)/analysis/ukb_pan_ancestry/figures/figure5/archived/pop_analysis_x_153764217_g6pd_all.txt.bgz', delim='\t')
  meta_data <- read_delim('/Users/wlu/Dropbox (Partners HealthCare)/analysis/ukb_pan_ancestry/figures/figure5/archived/meta_analysis_x_153764217_g6pd_all.txt.bgz', delim='\t')

  data <- result_data %>%
    mutate(type = 'pop') %>%
    select(phenocode, pop, type, Pvalue, BETA, SE) %>%
    rbind(meta_data %>%
      mutate(pop = 'Meta (raw)',
             type = 'Meta') %>%
      select(phenocode, pop, type, Pvalue, BETA, SE)
    )


  if(!no_meta_hq){
    data <- data %>%
      rbind(meta_data %>%
        mutate(pop = 'Meta (hq)',
               type = 'Meta') %>%
        select(phenocode, pop, type, Pvalue = Pvalue_hq, BETA= BETA_hq, SE=SE_hq)
      )
    n_meta = 2
  }else{
    data <- data %>%
      mutate(pop = if_else(pop == 'Meta (raw)', 'Meta', pop))
    pops <- c(pops[1:(length(pops)-2)], 'Meta')
    names(pops)[length(pops)] <- 'meta'
    n_meta = 1
  }

  data <- data %>%
    filter(complete.cases(.)) %>%
    mutate(pop = factor(pop, levels=rev(pops), labels=rev(names(pops))),
           phenocode = factor(phenocode,
                                 levels = rev(c(30270, 30300, 30750, 30070, 30010)),
                                 labels = rev(c('MCV*', 'RET*', 'HbA1c*', 'RDW', 'RBC'))))

  if(no_amr){
    data <- data %>%
      filter(pop != 'amr')
  }

  figure <- data %>%
    ggplot +
    geom_point(aes(x = BETA, y = phenocode, color = pop, size=pop),position=position_dodge(width = 1)) +
    geom_pointrange(aes(x = BETA, y = phenocode, color = pop, xmax = BETA + 1.96*SE, xmin = BETA-1.96*SE, size=pop),position=position_dodge(width = 1)) + themes +
    geom_hline(yintercept = c(1.5,2.5,3.5,4.5), lty=2, lwd=0.1) +
    labs(color='Ancestry group',  pch='Ancestry group', x = 'Beta', y = NULL) +
    scale_color_manual(name='Ancestry group', values = pop_colors,
                       breaks = rev(names(pops)),
                       labels = rev(pops)) +
    scale_size_manual(name='Ancestry group',
                      values = rev(c(rep(0.2, 6),rep(0.5,n_meta))),
                      breaks = rev(names(pops)),
                      labels = rev(pops)) +
    themes +
    theme_classic() +
    theme( axis.text = element_text(size = text_size*2.5, color = 'black'),
           axis.title = element_text(size = text_size*2.5,  face = 'bold'),
           legend.text = element_text(size = text_size*2),
           legend.title = element_text(size = text_size*2.5, face = 'bold'),
           legend.position = 'top'
           ) +
    guides(colour = guide_legend(nrow = 1),
           size = guide_legend(nrow = 1))

  print(figure)

  if(save){
    png(paste0(output_path, 'figure5c_g6pd_forest_plot',if_else(no_amr, '_no_amr', ''),if_else(no_meta_hq, '_no_meta_hq', ''),'.png'), width = width, height = height, units = 'in', res = 300)
    print(figure)
    dev.off()
  }
  return(figure)
}


figure5 <- function(output_format, output_path=OUTPUT){
  # Layout parameters:
  # image_widths = image_info(figure5b_raw)$width/300 # = 7.23
  # image_heights = image_info(figure5b_raw)$height/300 # = 5.27
  relative_heights = c(1, 1.5)
  factor = 0.9
  width = 7.23 * 2.25 * factor
  height = 5.27 * sum(relative_heights) * factor

  # Prepare panels:
  figure5a <- figure5a(text_size = 7, height=4, width=6, output_path=output_path, save=F)
  figure5b_path <- '~/PycharmProjects/ukbb_pan_ancestry/G6PD_pan_ukb_X_153764217_META_500000_1Mb_continuous-30270-both_sexes--irnt.pdf'
  figure5b_raw <- image_read_pdf(figure5b_path)
  figure5b <- ggplot() +
    annotation_raster(figure5b_raw, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    geom_blank() + theme_nothing()
  figure5c <- figure5c(text_size = 7, height=6, width=10, output_path=output_path, save=F, no_amr=T, no_meta_hq=T)

  # Arrange panels:
  top  = ggpubr::ggarrange(figure5a, figure5b, ncol = 2, labels = c("a", "b"), font.label = list(size = 20), widths = c(1,1.2))
  combined = ggpubr::ggarrange(top, figure5c,  ncol=1, labels = c("", "c"), font.label = list(size = 20))

  # Save the top two panels (sanity check):
  output_type(output_format, paste0(output_path, 'figure5_top.', output_format), height=height/2, width=width)
  print(top)
  dev.off()
  # ggsave(filename = paste0(output_path, 'figure5_top.', output_format),top, height=height/2, width=width)

  # Save the full figure:
  output_type(output_format, paste0(output_path, 'figure5.', output_format), height=height, width=width)
  print(combined)
  dev.off()
  # ggsave(filename = paste0(OUTPUT, 'figure5.', output_format),combined, height=height, width=width)

  return(combined)
}

figure5('pdf', OUTPUT)
figure5('png', OUTPUT)


# Other adjustments
figure5a <- figure5a(text_size = 7, height=4, width=6, output_path=OUTPUT, save=T)
figure5c <- figure5c(text_size = 7, height=6, width=10, output_path=OUTPUT, save=T)
figure5c <- figure5c(text_size = 7, height=6, width=10, output_path=OUTPUT, save=T, no_amr=T)
figure5c <- figure5c(text_size = 7, height=6, width=10, output_path=OUTPUT, save=T, no_meta_hq=T)
figure5c <- figure5c(text_size = 7, height=6, width=10, output_path=OUTPUT, save=T, no_amr=T, no_meta_hq=T)









