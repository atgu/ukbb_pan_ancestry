__author__ = 'Rahul Gupta'

# Contains relevant constants for running the LDSC pipeline for the Pan Ancestry project.

project_dir = '/Volumes/rahul/Projects/2020_ukb_diverse_pops/Experiments/200501_ldsc_div_pops_pipeline/Data/'
flat_file_location = 'gs://ukb-diverse-pops/sumstats_flat_files/'
bucket = 'rgupta-ldsc'
# ancestries = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']
output_bucket = f'gs://{bucket}/'

#anc_to_ldscore = lambda anc: f'gs://ukb-diverse-pops-public/ld_release/UKBB.{anc}'
anc_to_ldscore = lambda anc: f'gs://rgupta-ldsc/ld/UKBB.{anc}'