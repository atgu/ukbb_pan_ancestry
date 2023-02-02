library(tidyverse)

prefix <- "/Volumes/rahul/"
path_to_pheno_manifest <- paste0(prefix, 'Projects/2020_ukb_diverse_pops/Experiments/',
                                 '200501_ldsc_div_pops_pipeline/Data/',
                                 'phenotype_manifest.tsv')
path_to_round2_h2 <- paste0(prefix,'Projects/2020_ukb_diverse_pops/Experiments/',
                            '200715_compare_h2_results_gnomad_ldsc/Data/',
                            'ukb31063_h2_topline.02Oct2019.tsv')
path_to_round1_h2 <- paste0(prefix,'Projects/2020_ukb_diverse_pops/Experiments/',
                            '200715_compare_h2_results_gnomad_ldsc/Data/',
                            'ukbb_all_h2univar_results.txt')
output <- 'RG_pan_ancestry_round2_round1_map.tsv'
# --------------------------------

paste_noempty <- function(..., sep='_') {
  # Acts like paste(), but assumes that each input is the same length.
  # Less efficient, but excludes any NA values without failure or conversion to string.
  #
  # INPUTS:
  # ...: vectors of equal lengths to be merged across
  # sep: delimiter (default _)
  # 
  # OUTPUTS: merged vector of the same length as the input vectors
  #
  # Example:
  # paste_noempty(c('a',NA), c('b','c'), sep='_') -> c('a_b','c')
  fl = list(...)
  if(length(unique(sapply(fl, length))) != 1) stop('Input vectors must each be the same length.')
  pasted_vec <- sapply(1:length(fl[[1]]), function(idx) {
    vec_for_paste <- sapply(fl, function(x)x[idx])
    return(paste0(vec_for_paste[!is.na(vec_for_paste)],collapse=sep))
  })
  return(pasted_vec)
}

# Combination of phenocode, coding, and modifier is sufficient to match to previous phenotype IDs
# Import phenotype manifest for UKB Pan Ancestry:
phenos <- read_tsv(path_to_pheno_manifest, guess_max = 7000) %>%
  transmute(trait_type, phenocode, pheno_sex, coding, modifier, description,
            description_more, coding_description, category) %>%
  mutate(matching_ID = paste_noempty(phenocode, coding, modifier))

# Import Round 1 and Round 2 phenotype IDs via the heritability tables. Note that
# these are the same identifiers used in the full Round 2 manifest.
ukb_31063_h2 <- read_tsv(path_to_round2_h2) %>%
  transmute(round2_pheno = phenotype, round2_description = description) %>%
  mutate(in_round2 = T)
ukb_h2_round_1 <- read_tsv(path_to_round1_h2) %>%
  transmute(round1_pheno = phenotype, round1_description = description) %>%
  mutate(in_round1 = T)

phenos_augmented <- phenos %>%
  left_join(y = ukb_31063_h2, by = c('matching_ID' = 'round2_pheno')) %>%
  left_join(y = ukb_h2_round_1, by = c('matching_ID' = 'round1_pheno')) %>%
  mutate(in_round2 = ifelse(is.na(in_round2), F, T),
         in_round1 = ifelse(is.na(in_round1), F, T))

write_tsv(phenos_augmented, output)