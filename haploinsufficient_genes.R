source('~/ukbb_pan_ancestry/constants.R')
setwd('~/Dropbox (Partners HealthCare)/analysis/ukb_pan_ancestry/')

efo = read_tsv(gzfile('data/known_ukbb_loci.meta_hq_annotated.txt.bgz'),
               col_types = cols(locus=col_character(),
                                lead_locus=col_character(),
                                locus_otg=col_character()))

novel_hits = efo %>%
  dplyr::group_by(trait_efo_category, idx) %>% #  trait_type, phenocode, pheno_sex, coding, modifier, lead_locus, nearest_genes) %>%
  dplyr::mutate(any_same_category=any(trait_efo_category == trait_efo_category_otg, na.rm=T)) %>%
  dplyr::filter(!any_same_category) %>%
  dplyr::select(matches("lead_locus") | (!contains("otg") & !matches("locus"))) %>%
  distinct %>% ungroup

gene_lists = load_all_gene_list_data()

## Combine gene list information and filter to haploinsufficient gene related hits
x <- novel_hits %>%
  dplyr::rename(gene=nearest_genes) %>%
  separate_rows(gene) %>%
  left_join(gene_lists, multiple = 'all') %>%
  filter(grepl('haploinsufficiency', gene_list))
write_csv(x, 'novel_hits_haploinsufficient.csv')


x <- read_csv('data/novel_hits_haploinsufficient.csv')

## Annotate OMIM and phenotype information
# FROM github repo: https://github.com/macarthur-lab/gene_lists/tree/master/other_data
omim_info <- read_tsv('~/gene_lists/other_data/omim.use.tsv')
colnames(omim_info) <- paste0('omim_', colnames(omim_info))
omim_info <- omim_info %>%
  mutate(gene = strsplit(omim_genes, "\\|")) %>%
  unnest(gene)

# FROM CODE BELOW:
# result_mt = load_final_sumstats_mt(
#     filter_phenos=False,
#     filter_variants=False,
#     filter_sumstats=False,
#     annotate_with_nearest_gene=False,
#     filter_pheno_h2_qc=False,
#     exponentiate_p=True
# )
# pheno_ht = result_mt.select_cols('description', 'description_more', 'coding_description', 'category')
# pheno_ht = result_mt.cols()
# pheno_ht.export(f'{my_bucket}/pan_ukb_phenotype_info.txt.bgz')
pheno_info <- read_delim('~/Downloads/pan_ukb_phenotype_info.txt.bgz', delim = '\t') %>%
  filter(pop_index==0) %>%
  select('trait_type', 'phenocode', 'pheno_sex', 'coding', 'modifier', 'description', 'description_more', 'coding_description')

## Filter to high quality variants
x_hq <- x %>%
  filter(high_quality)

x_hq_omim_annt <- x_hq %>%
  mutate(chrom = str_split(lead_locus, ':') %>% map_chr(., 1))%>%
  select(chrom, lead_locus, rsid, gene, trait_type, phenocode, pheno_sex, coding, modifier, lead_pvalue, AFR_af, AMR_af, CSA_af, EAS_af, EUR_af, MID_af, gene_list) %>%
  merge(., pheno_info , by = c('trait_type', 'phenocode', 'pheno_sex', 'coding', 'modifier'), all.x=T) %>%
  merge(., omim_info, by = 'gene', all.x = T)

write_tsv(x_hq_omim_annt , 'haploinsufficient_selected_hq_omim_info.tsv')