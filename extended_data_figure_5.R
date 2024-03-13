source('~/ukbb_pan_ancestry/figure3.R')

known_novel_efo = function() {
  # get_efo_stats(efo)
  # All:
  # Total: 85,960
  # 14,588 (no term)
  # Novelty: novel = 35422, known = 35,950
  # Category: novel = 5,551, known = 65,821
  
  # X chromsosomes:
  # 2,920 total
  # 36 with terms matching
  # 473 (no term)
  # novel terms on X = 2411 (82%)
  # 2920-473-36 = 2411
  # 245 with categories matching on X
  # novel categories on X = 2202 (75%)
  # 2920-473-245 = 2202
  
  efo %>%
    group_by(trait_efo_category, idx) %>%
    summarize(any_same_category=any(trait_efo_category == trait_efo_category_otg)) %>%
    group_by(trait_efo_category) %>%
    summarize(total_count=n(),
              known_count=sum(any_same_category, na.rm=T),
              novel_count=total_count-known_count,
              known_frac=known_count/total_count,
              novel_frac=1-known_frac) %>%
    mutate(trait_efo_category = fct_reorder(trait_efo_category, novel_frac)) %>%
    pivot_longer(-trait_efo_category) %>%
    separate(name, c('novelty', 'scale'), ('_')) %>%
    filter(!is.na(trait_efo_category)) -> data
  
  # data %>%
  #   mutate(value=if_else(scale == 'frac', value, -value)) %>%
  #   ggplot() +
  #   aes(y = value, x = trait_efo_category, fill = novelty) +
  #   geom_bar(stat='identity') +
  #   facet_share(~scale, dir = "h", scales = "free", reverse_num = T,
  #               labeller = labeller(scale = c('count'='Number of hits', 'frac'='Fraction novel'))) + 
  #   coord_flip() +
  #   theme_classic() +
  #   ylab(NULL) + xlab(NULL) +
  #   scale_fill_discrete(labels=c('known' = 'Known', 'novel' = 'Novel'), name=NULL)
  
  fill_colors = c(
    'known' = muted('lightblue'),
    'novel' = muted('orange', l=70, c=100)
  )
  data %>%
    filter(scale == 'count' & novelty == 'novel') %>% 
    arrange(desc(value))
  p_num = data %>%
    filter(scale == 'count' & novelty != 'total') %>%
    ggplot +
    aes(x = value, y = trait_efo_category, fill = novelty, label=comma(value)) +
    geom_bar(stat='identity', position='stack') +
    geom_text(data=data %>% filter(novelty == "total"), hjust=1, nudge_x = -200, size=2.5) +
    theme_classic() + scale_x_reverse(labels=comma, name='Number of significant loci',
                                      expand = expansion(c(0.2, 0))) +
    scale_y_discrete(position='right') +
    theme(axis.text.y=element_blank(),
          legend.position=c(0.01, 0.99),
          legend.justification=c(0, 1)) +
    ylab(NULL) +
    scale_fill_manual(labels=c('known' = 'Known', 'novel' = 'Novel'), name=NULL,
                      values=fill_colors)
  
  data2 = data %>%
    filter(novelty != 'total' & scale == 'count') %>%
    arrange(desc(novelty)) %>%
    group_by(trait_efo_category) %>%
    mutate(frac = value / sum(value),
           x =  cumsum(frac) - frac / 2) %>% ungroup
  
  p_frac = data2 %>%
    ggplot +
    aes(x = frac, y = trait_efo_category, fill = novelty) +
    geom_col() +
    geom_text(aes(x=x, label=percent(frac, accuracy=1), color=novelty),
              data=data2 %>% filter(frac > 0.1),
              hjust=0.5, size=2.5) +
    theme_classic() + scale_x_continuous(labels=percent, name='Fraction novel') +
    theme(axis.text.y=element_blank()) +
    scale_color_manual(values=c('known' = 'white', 'novel' = 'black'), guide=F) +
    ylab(NULL) + guides(fill=F) +
    scale_fill_manual(labels=c('known' = 'Known', 'novel' = 'Novel'), name=NULL,
                      values=fill_colors)
  # p_frac
  p_labels = data %>%
    filter(novelty == 'known' & scale == 'count') %>%
    ggplot + aes(y = trait_efo_category, label = trait_efo_category) +
    geom_text(x = 0.5, size = 3, hjust = 0.5) + 
    theme(axis.text.y=element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank()) + xlab(NULL) + ylab(NULL)
  
  p = egg::ggarrange(p_num + theme(plot.margin = margin(0,0,0,8, "pt")),
                     p_labels + theme(plot.margin = margin(0,0,0,0, "pt")),
                     p_frac + theme(plot.margin = margin(0,8,0,0, "pt")),
                     nrow=1, ncol=3)
  # return(p)
  return(list(p_num + theme(plot.margin = margin(0,0,0,8, "pt")),
              p_labels + theme(plot.margin = margin(0,0,0,0, "pt")),
              p_frac + theme(plot.margin = margin(0,8,0,0, "pt"))))
}

ed5 = known_novel_efo()

extended_data_figure5 = function(output_format = 'png') {
  output_type(output_format, paste0('extended_data_figure5.', output_format), height=3.5, width=7.5)
  print(egg::ggarrange(plots=ed5, ncol = 3))
  dev.off()
}
extended_data_figure5()
