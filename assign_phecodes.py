#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import uuid
import hail as hl
from hail.utils import hadoop_open
from ukbb_pan_ancestry import *
from ukb_common import *
hl.init()


def checkpoint_tmp(hail_obj, tmppath='gs://ukbb-diverse-temp-30day/', tmpname=None, overwrite=True):
    if tmpname is None:
        tmpname = str(uuid.uuid4())
    return hail_obj.checkpoint(tmppath + tmpname, overwrite=overwrite)


# Read in the Phecode (v1.2b1) <-> ICD 9/10 codes mapping
with hadoop_open('gs://ukb-diverse-pops/phecode/UKB_Phecode_v1.2b1_ICD_Mapping.txt', 'r') as f:
    df = pd.read_csv(f, delimiter='\t', dtype=str)
list_of_icd_codes_to_include = [row.icd_codes.split(',') for _, row in df.iterrows()]
list_of_phecodes_to_exclude = [row.exclude_phecodes.split(',') for _, row in df.iterrows()]
df['icd_codes'] = list_of_icd_codes_to_include
df['exclude_phecodes'] = list_of_phecodes_to_exclude

# Convert it to HailTable
phecode_ht = hl.Table.from_pandas(df)
phecode_ht = phecode_ht.key_by('icd_codes')
phecode_ht = phecode_ht.checkpoint('gs://ukb-diverse-pops/phecode/UKB_Phecode_v1.2b1_ICD_Mapping.ht', overwrite=True)

# Retreive UKB ICD MatrixTable and combine codes based on Phecode definitions
icd_all = hl.read_matrix_table(get_ukb_pheno_mt_path('icd_all'))
mt = combine_phenotypes(icd_all,
                        icd_all.icd_code,
                        icd_all.any_codes,
                        list_of_icd_codes_to_include,
                        new_col_name='icd_codes',
                        new_entry_name='include_to_cases')
mt = mt.annotate_cols(phecode=phecode_ht[mt.icd_codes].phecode,
                      phecode_sex=phecode_ht[mt.icd_codes].sex,
                      exclude_phecodes=phecode_ht[mt.icd_codes].exclude_phecodes)

# Annotate sex for sex-specific phenotypes
ukb_pheno_ht = hl.read_table(get_ukb_pheno_ht_path())
mt = mt.annotate_rows(isFemale=ukb_pheno_ht[mt.userId].sex == 0)
mt = checkpoint_tmp(mt)

# Compute phecode excluded from controls
mt = mt.key_cols_by()
exclude_mt = combine_phenotypes(mt,
                                mt.phecode,
                                mt.include_to_cases,
                                list_of_phecodes_to_exclude,
                                new_entry_name='exclude_from_controls')
exclude_mt = checkpoint_tmp(exclude_mt)

# Annotate exclusion
mt = mt.key_cols_by('exclude_phecodes')
mt = mt.annotate_entries(exclude_sex=hl.switch(mt.phecode_sex).when("Male",
                                                                    mt.isFemale).when("Female",
                                                                                      ~mt.isFemale).default(False),
                         exclude_from_controls=exclude_mt[mt.userId, mt.exclude_phecodes].exclude_from_controls)

# Compute final case/control status
# `case_control` becomes missing (NA) if a sample 1) is excluded because of sex, 2) is not cases and excluded from controls.
mt = mt.annotate_entries(
    case_control=hl.if_else(mt.exclude_sex |
                            (~mt.include_to_cases & mt.exclude_from_controls), hl.null(hl.tbool), mt.include_to_cases))

mt = mt.key_cols_by('phecode')
mt.describe()

mt.write(get_ukb_pheno_mt_path('phecode'), overwrite=True)
