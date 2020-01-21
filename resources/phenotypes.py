
bucket = 'gs://ukb-diverse-pops'
pheno_folder = f'{bucket}/Phenotypes/Everyone'


def get_phesant_all_phenos_tsv_path(sex: str):
    return f'{pheno_folder}/PHESANT_final_output/January_2020/phesant_output_multi_ancestry_combined_{sex}.tsv.gz'
    # if sex == 'both_sexes':
    #     return f'{pheno_folder}/uk_round2_allSamples_phenos_phesant_QC.{{}}.tsv.gz'
    # else:
    #     return f'{pheno_folder}/SexSpecific/uk_round2_allSamples_phenos_phesant_QC_Sex{int(sex == "males")}.{{}}.tsv.gz'


def get_ukb_phesant_summary_tsv_path(sex: str = 'both_sexes_no_sex_specific'):
    return f'{pheno_folder}/PHESANT_final_output/January_2020/phesant_output_multi_ancestry_combined_{sex}_summary.tsv'


phesant_biomarker_phenotypes_tsv_path = f'{pheno_folder}/uk_round2_allSamples_biomarkers_phesant_QC.tsv.gz'

pre_phesant_tsv_path = f'{pheno_folder}/neale_lab_parsed_QC_Oct2019.tsv'
pre_phesant_biomarkers_tsv_path = f'{pheno_folder}/neale_lab_parsed_biomarkers.tsv'


def get_ukb_pheno_ht_path(sex: str = 'both_sexes'):
    return f'{pheno_folder}/ht/all_pops_{sex}.ht'


def get_ukb_pheno_mt_path(data_type: str = None, sex: str = 'both_sexes_no_sex_specific'):
    if data_type is None:
        data_type, sex = 'full', 'full'
    return f'{pheno_folder}/mt/{sex}/all_pops_{data_type}.mt'


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


def get_phenotype_summary_path(data_type: str, extension = 'ht'):
    return f'{pheno_folder}/summary/phenos_{data_type}.{extension}'
