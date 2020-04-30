
from gnomad.utils import *
from ukb_common import *
from ukbb_pan_ancestry import *

output_file = 'gs://ukb-diverse-pops/misc/l_dopa_melatoma/co_occurrence.ht'

def main(args):
    hl.init(log='/ldopa.log')

    if args.get_n_carriers:
        pheno_mt = hl.read_matrix_table(get_ukb_pheno_mt_path())
        print(pheno_mt.count())
        mt = pheno_mt.filter_cols(
            (pheno_mt.phenocode == "dopamine pro-drug _ decarboxylase inhibitor|Parkinson's") |
            (pheno_mt.phenocode == "C43") | (pheno_mt.phenocode == "G20") | (pheno_mt.phenocode == "1717")
        )
        mt = mt.annotate_entries(both_sexes=hl.case()
                                 .when((mt.pheno == "1717") & (mt.both_sexes <= 2.0), 1.0)  # Coding for "very fair" or "fair":
                                 # http://biobank.ctsu.ox.ac.uk/showcase/coding.cgi?id=100431
                                 .when(mt.pheno == "1717", 0.0)
                                 .when((mt.pheno == "1747") & (mt.both_sexes == 2.0), 1.0)  # Coding for "red":
                                 # http://biobank.ctsu.ox.ac.uk/showcase/coding.cgi?id=100434
                                 .when(mt.pheno == "1747", 0.0)
                                 .default(mt.both_sexes))
        mt = mt.filter_rows(hl.agg.any((mt.both_sexes == 1.0) & (mt.pheno == "1717")))
        print(mt.count())
        ht = mt.localize_entries('entries', 'phenos')
        res = ht.aggregate(hl.agg.array_agg(
            lambda i: hl.agg.array_agg(
                lambda j: hl.agg.count_where((ht.entries[i].both_sexes == 1.0) & (ht.entries[j].both_sexes == 1.0)),
                hl.range(hl.len(ht.phenos))
            ),
            hl.range(hl.len(ht.phenos))), _localize=False)
        phenos = ht.index_globals().phenos
        res_ht = hl.utils.range_table(1).annotate(result=hl.zip(phenos, res)).explode('result').key_by().drop('idx')
        res_ht = res_ht.transmute(pheno1=res_ht.result[0], result=hl.zip(phenos, res_ht.result[1])).explode('result')
        res_ht = res_ht.transmute(pheno2=res_ht.result[0], n_shared=res_ht.result[1])
        res_ht = res_ht.checkpoint(output_file, args.overwrite)
        res_ht.select(pheno1=res_ht.pheno1.meaning, pheno2=res_ht.pheno2.meaning, n_shared=res_ht.n_shared).show(40, width=200)

    def seconds_to_years(seconds, round=True):
        years = seconds / 60 / 60 / 24 / 365.25
        return hl.int(years) if round else years

    if args.get_order:
        prescription_mt = hl.read_matrix_table(get_ukb_pheno_mt_path('prescriptions'))
        prescription_ht = prescription_mt.filter_cols(
            prescription_mt.Drug_Category_and_Indication == "dopamine pro-drug / decarboxylase inhibitor,Parkinson's"
        ).localize_entries('l_dopa', 'phenos')
        prescription_ht = prescription_ht.transmute(l_dopa=hl.or_missing(hl.len(prescription_ht.l_dopa[0].values) > 0,
                                                                         prescription_ht.l_dopa[0].values[0]))
        diag_mt = hl.read_matrix_table(get_hesin_mt_path('diag'))
        diag_ht = diag_mt.filter_cols(diag_mt.diag_icd10.startswith('C43')).localize_entries('melanoma', 'phenos')
        ht = diag_ht.annotate(l_dopa=prescription_ht[diag_ht.key].l_dopa)
        first_diagnosis = hl.sorted(hl.zip(ht.phenos, ht.melanoma), key=lambda x: x[1].admidate)[0]
        ht = ht.annotate(pheno=first_diagnosis[0], melanoma=first_diagnosis[1])
        ht = ht.filter(hl.is_defined(ht.melanoma.ins_index) & hl.is_defined(ht.l_dopa.issue_date)).persist()
        ht = ht.annotate(melanoma_age=seconds_to_years(ht.melanoma.admidate - ht.dob),
                         ldopa_age=seconds_to_years(ht.l_dopa.issue_date - ht.dob))
        ht.flatten().show(width=200)

    # TODO:
    #  red hair? (how many melanoma?), add light olive back in?
    #  compare date of (Parkinson's | l-dopa initiation) vs. melanoma diagnosis
    #  number of melanomas (6 months apart, and/or different site) after the Parkinson's
    #  control: Alzheimer's vs. colon cancer (C18)
    #  add self-reported or family Parkinson's - are we seeing overlap between skin tone and risk?

    #  TODO:
    #   suppressed vitamin D levels will lead to higher propensity of opioid addition
    #   correlate vitamin D levels with: length of taking opioids, or methodone?, cigarette/etoh?
    #   opioid use; depression; vitamin D level, intake, supplement (control by other supplements?)
    #   censor by post-fracture?


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--get_n_carriers', help='Overwrite everything', action='store_true')
    parser.add_argument('--get_order', help='Overwrite everything', action='store_true')
    parser.add_argument('--abo_covid', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)