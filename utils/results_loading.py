import re
from ukbb_common import *
from ukbb_pan_ancestry.resources import *


def get_phenos_to_run(pop: str, limit: int = None, pilot: bool = False, single_sex_only: bool = False,
                      specific_phenos: str = '', skip_case_count_filter: bool = False, first_round_phenos: bool = False,
                      sex_stratified: str = None):
    ht = hl.read_table(get_phenotype_summary_path('full'))
    ht = ht.filter(ht.pop == pop)
    min_cases = MIN_CASES_EUR if pop == 'EUR' else MIN_CASES

    criteria = True
    if first_round_phenos:
        trait_types = hl.literal({'categorical', 'continuous', 'biomarkers', 'icd10', 'prescriptions', 'phecode'})
        criteria = (trait_types.contains(ht.trait_type) &
                    (ht.modifier != 'raw') &
                    ~hl.set({'22601', '22617', '20024', '41230', '41210'}).contains(ht.phenocode) &
                    (ht.n_cases_by_pop >= min_cases))
    if not skip_case_count_filter:
        criteria &= (ht.n_cases_by_pop >= min_cases)

    if single_sex_only:
        prop_female = ht.n_cases_females / (ht.n_cases_males + ht.n_cases_females)
        criteria &= ((prop_female <= 0.1) | (prop_female >= 0.9))

    ht = ht.filter(criteria).key_by()

    if sex_stratified:
        ht_sex_specific = ht.annotate(pheno_sex='males').union(ht.annotate(pheno_sex='females'))
        if sex_stratified == 'all':
            ht = ht.union(ht_sex_specific)
        else:
            ht = ht_sex_specific

    output = set([tuple(x[field] for field in PHENO_KEY_FIELDS) for x in ht.select(*PHENO_KEY_FIELDS).collect()])
    if pilot:
        output = output.intersection(PILOT_PHENOTYPES)
    if specific_phenos:
        specific_phenos = specific_phenos.split(',')
        output = [x for x in output if all(map(lambda y: y is not None, x)) and any([re.match(pcd, '-'.join(x)) for pcd in specific_phenos])]
    if limit:
        output = set(sorted(output)[:limit])

    pheno_key_dict = [dict(zip(PHENO_KEY_FIELDS, x)) for x in output]
    if first_round_phenos:
        pheno_key_dict = recode_pkd_to_legacy(pheno_key_dict)
    return pheno_key_dict


def get_pheno_dict():
    pheno_ht = hl.read_matrix_table(get_ukb_pheno_mt_path()).cols()
    pheno_dict = create_broadcast_dict(pheno_ht.key)
    return pheno_dict


def get_heritability_dict(pop):
    heritability_ht = hl.import_table(get_heritability_txt_path(), impute=True, key=PHENO_KEY_FIELDS)
    heritability_ht = heritability_ht.filter(heritability_ht.pop == pop)
    heritability_dict = create_broadcast_dict(heritability_ht.key)
    return heritability_dict


def get_modified_key(mt):
    key = mt.col_key.annotate(phenocode=format_pheno_dir(mt.phenocode),
                              modifier=hl.case(missing_false=True)
                              .when(mt.trait_type == "biomarkers", "")
                              .when((mt.pop == 'EAS') & (mt.phenocode == '104550'), '104550')
                              .default(mt.modifier))
    return key