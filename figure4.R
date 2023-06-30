source('~/ukbb_pan_ancestry/constants.R')

efo = read_tsv(gzfile('data/known_ukbb_loci.meta_hq_annotated.txt.bgz'),
               col_types = cols(locus=col_character(),
                                lead_locus=col_character(),
                                locus_otg=col_character()))

novel_hits = efo %>%
  group_by(trait_efo_category, idx) %>% # , trait_type, phenocode, pheno_sex, coding, modifier, lead_locus, nearest_genes) %>%
  mutate(any_same_category=any(trait_efo_category == trait_efo_category_otg, na.rm=T)) %>%
  filter(!any_same_category) %>%
  select(matches("lead_locus") | (!contains("otg") & !matches("locus"))) %>%
  distinct %>% ungroup

novel_hits %>% count(trait_efo_category)

novel_hits %>% 
  mutate(freq_max = pmax(AFR_af, AMR_af, CSA_af, EAS_af, MID_af),
         fc = freq_max / EUR_af) %>% 
  filter(fc > 2) %>%
  filter(!is.na(trait_efo_category) & trait_efo_category != 'Other measurement') %>%
  filter(trait_type %in% c('icd10', 'phecode')) %>%
  arrange(desc(freq_max)) %>% View

novel_hits %>% 
  filter(EUR_af < 0.05 & AFR_af > 0.05) %>%
  filter(!is.na(trait_efo_category) & trait_efo_category != 'Other measurement') %>%
  count(trait_type, phenocode, pheno_sex, coding, modifier) %>% count