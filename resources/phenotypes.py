
bucket = 'gs://ukb-diverse-pops'
pheno_folder = f'{bucket}/Phenotypes/Everyone'


def get_phesant_all_phenos_tsv_paths(sex: str):
    if sex == 'both_sexes':
        return f'{pheno_folder}/uk_round2_allSamples_phenos_phesant_QC.{{}}.tsv.gz'
    else:
        return f'{pheno_folder}/SexSpecific/uk_round2_allSamples_phenos_phesant_QC_Sex{int(sex == "males")}.{{}}.tsv.gz'

phesant_biomarker_phenotypes_tsv_path = f'{pheno_folder}/uk_round2_allSamples_biomarkers_phesant_QC.tsv.gz'

pre_phesant_tsv_path = f'{pheno_folder}/neale_lab_parsed_QC_Oct2019.tsv'
pre_phesant_biomarkers_tsv_path = f'{pheno_folder}/neale_lab_parsed_biomarkers.tsv'


def get_ukb_pheno_ht_path(pop: str = 'all_pops', sex: str = 'both_sexes'):
    return f'{pheno_folder}/ht/{pop}_{sex}.ht'


def get_ukb_pheno_mt_path(data_type: str, pop: str = 'all_pops', sex: str = 'both_sexes'):
    return f'{pheno_folder}/{sex}/{pop}_{data_type}.mt'


def get_biomarker_ht_path(pop: str = 'all_pops', sex: str = 'both_sexes'):
    return f'{pheno_folder}/ht/biomarkers_{pop}_{sex}.ht'


def get_hesin_raw_data_path(data_type: str = None):
    assert data_type in (None, 'diag', 'oper', 'delivery', 'maternity', 'psych')
    data_type = "_" + data_type if data_type else ""
    return f'gs://ukb31063/ukb31063.hesin{data_type}.20191008.txt'

prescription_tsv_path = 'gs://ukb31063/ukb31063.gp_scripts.20191008.txt'
prescription_mapping_path = f'{bucket}/Phenotypes/ukb_prescription_mapping.tsv'


def get_hesin_mt_path(data_type: str):
    assert data_type in ('diag', 'oper')
    return f'{pheno_folder}/hesin/{data_type}.mt'


def get_hesin_delivery_ht_path():
    return f'{pheno_folder}/hesin/delivery.ht'


def get_gp_data_tsv_path(data_type: str = None):
    assert data_type in ('registrations', 'clinical', 'scripts')
    return f'gs://ukb31063/ukb31063.gp_{data_type}.20191008.txt'


def get_phenotype_summary_tsv_path(data_type: str):
    return f'{pheno_folder}/summary/phenos_{data_type}.tsv'


# TODO: add this file
def get_ukb_phesant_summary_tsv_path(pop: str = 'all_pops', sex: str = 'both_sexes'):
    return f'gs://phenotype_pharma/misc_files/phesant_output_combined_{sex}_summary.modified.tsv'

