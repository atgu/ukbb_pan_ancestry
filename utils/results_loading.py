from ukb_common import *
from ukbb_pan_ancestry.resources import *


def get_phenos_to_run(pop: str, limit: int = None, pilot: bool = False, single_sex_only: bool = False,
                      specific_phenos: str = '', skip_case_count_filter: bool = False, first_round_phenos: bool = False):
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

    ht = ht.filter(criteria)

    output = set([tuple(x[field] for field in PHENO_KEY_FIELDS) for x in ht.key_by().select(*PHENO_KEY_FIELDS).collect()])
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
