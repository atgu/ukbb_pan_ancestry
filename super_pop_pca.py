from gnomad_hail import *
from ukbb_pan_ancestry import *
import hail as hl
import argparse
import pandas as pd


def load_ref(dirname, basename):
    """
    Loads a reference plink dataset, writes out a matrix table
    :param dirname: plink file directory name
    :param basename: plink base filename
    :return:
    """
    ref = hl.import_plink(bed=dirname + basename + '.bed',
                          bim=dirname + basename + '.bim',
                          fam=dirname + basename + '.fam',
                          min_partitions=100)
    ref.describe()

    print('sites in ref data: ' + str(ref.count()))  # (639590, 3547)
    ref.write(dirname + basename + '.mt', args.overwrite)


def intersect_ref(dirname, basename, ukbb):
    """
    Intersects reference panel with UKBB data and writes intersections as matrix tables
    :param dirname: directory name to put reference and ukbb file intersections
    :param basename: base filename for reference data
    :param ukbb: ukbb data
    :return:
    """
    print(dirname + basename)
    this_ref = hl.read_matrix_table(dirname + basename + '.mt')

    # filter ukbb to sites in ref & array data
    ukbb_in_ref = ukbb.filter_rows(hl.is_defined(this_ref.rows()[ukbb.row_key]))
    print('sites in ref and UKBB data, inds in UKBB: ' + str(ukbb_in_ref.count()))  # (64233, 488377)

    ##
    ukbb_in_ref.write(dirname + 'intersect_ukbb_' + basename + '.mt', args.overwrite)

    # filter ref to ukbb sites
    ref_in_ukbb = this_ref.filter_rows(hl.is_defined(ukbb.rows()[this_ref.row_key]))
    print('sites in ref and UKBB data, inds in ref: ' + str(ref_in_ukbb.count()))  # (64233, 3547)
    ##
    ref_in_ukbb.write(dirname + 'intersect_' + basename + 'ukbb.mt', args.overwrite)


def run_pca(mt: hl.MatrixTable, out_prefix: str, overwrite: bool = False):
    """
    Run PCA on a dataset
    :param mt: dataset to run PCA on
    :param out_prefix: directory and filename prefix for where to put PCA output
    :return:
    """
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(mt.GT, k=20, compute_loadings=True)
    pca_mt = mt.annotate_rows(pca_af=hl.agg.mean(mt.GT.n_alt_alleles()) / 2)
    pca_loadings = pca_loadings.annotate(pca_af=pca_mt.rows()[pca_loadings.key].pca_af)

    pca_scores.write(out_prefix + 'scores.ht', overwrite)
    pca_scores = hl.read_table(out_prefix + 'scores.ht')
    pca_scores = pca_scores.transmute(**{f'PC{i}': pca_scores.scores[i - 1] for i in range(1, 21)})
    pca_scores.export(out_prefix + 'scores.txt.bgz')  # individual-level PCs

    pca_loadings.write(out_prefix + 'loadings.ht', overwrite)  # PCA loadings


def get_relatedness_path(pop, unrelated: bool = False, extension: str = 'mt'):
    if extension != '': extension = f'.{extension}'
    return f'gs://ukb-diverse-pops/pca/relatedness/{pop}_{"un" if unrelated else ""}rel{extension}'


def project_individuals(pca_loadings, project_mt):
    """
    Project samples into predefined PCA space
    :param pca_loadings: existing PCA space
    :param project_mt: matrix table of data to project
    :param project_prefix: directory and filename prefix for where to put PCA projection output
    :return:
    """
    ht_projections = pc_project(project_mt, pca_loadings)
    ht_projections = ht_projections.transmute(**{f'PC{i}': ht_projections.scores[i - 1] for i in range(1, 21)})
    return ht_projections


def main(args):
    if args.load_ref:
        load_ref(args.dirname, args.basename)

    if args.load_ukbb:
        samples = hl.read_table('gs://ukb-diverse-pops/pigmentation_phenos_covs_pops.ht')
        ukbb = hl.read_matrix_table('gs://ukb31063/ukb31063.genotype.mt')
        ukbb = ukbb.annotate_cols(**samples[ukbb.s])

    if args.intersect_ref:
        intersect_ref(args.dirname, args.basename, ukbb)

    if args.pca_project:
        """
        Compute PCA in global reference panel, project UKBB individuals into PCA space
        """
        ref_in_ukbb = hl.read_matrix_table(args.dirname + 'intersect_' + args.basename + 'ukbb.mt')
        print('Computing reference PCs')
        run_pca(ref_in_ukbb, args.out_prefix + args.basename + '_ukbb_')

        # project ukbb
        pca_loadings = hl.read_table(f'{args.out_prefix}{args.basename}_ukbb_loadings.ht')
        project_mt = hl.read_matrix_table(args.dirname + 'intersect_ukbb_' + args.basename + '.mt')
        ht = project_individuals(pca_loadings, project_mt)
        ht.export(args.out_prefix + 'ukbb_' + args.basename + '_scores.txt.bgz')

    # if args.continental_pca:
    #     """
    #     Compute PCA within reference panel super pops, project UKBB individuals into PCA space
    #     1. Filter UKBB to individuals in continental population
    #     2. Run PCA on continental ref
    #     3. Project UKBB inds
    #     """
    #     pass

    if args.ukbb_pop_pca:
        """
        Compute PCA in each UKBB population (unrelateds), project reference individuals and relateds into PCA space
        1. Filter UKBB to individuals in continental population
        2. Run PC-relate on these individuals
        # New
        2.5 Filter to pruned set of individuals
        #
        3. Filter UKBB population to unrelated individuals
        4. Run PCA on UKBB unrelateds within population
        5. Project relateds
        """

        for pop in POPS:
            mt = hl.read_matrix_table(get_ukb_grm_mt_path(pop))
            pruned_ht = hl.read_table(get_ukb_grm_pruned_ht_path(pop))
            mt = mt.filter_rows(hl.is_defined(pruned_ht[mt.row_key]))

            # run PC-relate
            if args.overwrite or not hl.hadoop_exists(get_relatedness_path(pop, extension='all_scores.ht/_SUCCESS')):
                _, scores, _ = hl.hwe_normalized_pca(mt.GT, k=10, compute_loadings=False)
                scores.write(get_relatedness_path(pop, extension='all_scores.ht'), args.overwrite)
            scores = hl.read_table(get_relatedness_path(pop, extension='all_scores.ht'))
            mt = mt.annotate_cols(scores=scores[mt.col_key].scores)
            # For EUR, required highmem machines with SSDs (Needed ~6T of hdfs space, so 20 workers + 100 pre-emptibles ran in ~7 hours)
            relatedness_ht = hl.pc_relate(mt.GT, min_individual_maf=0.05, scores_expr=mt.scores,
                                          min_kinship=0.05, statistics='kin',
                                          block_size=4096 if pop == 'EUR' else 512).key_by()
            relatedness_ht.write(get_relatedness_path(pop, extension='ht'), args.overwrite)
            relatedness_ht = hl.read_table(get_relatedness_path(pop, extension='ht'))

            # identify individuals in pairs to remove
            related_samples_to_remove = hl.maximal_independent_set(relatedness_ht.i, relatedness_ht.j, False)
            mt_unrel = mt.filter_cols(hl.is_defined(related_samples_to_remove[mt.col_key]), keep=False)
            mt_rel = mt.filter_cols(hl.is_defined(related_samples_to_remove[mt.col_key]), keep=True)

            mt_unrel.write(get_relatedness_path(pop, True, 'mt'), args.overwrite)
            mt_rel.write(get_relatedness_path(pop, extension='mt'), args.overwrite)

    if args.ukb_prune_pca_project:
        for pop in POPS:
            mt_unrel = hl.read_matrix_table(get_relatedness_path(pop, True, 'mt'))
            mt_rel = hl.read_matrix_table(get_relatedness_path(pop, extension='mt'))

            # Removing individuals
            pruned_inds = hl.import_table(get_pruned_tsv_path(), key='s')
            mt_rel = mt_rel.filter_cols(hl.is_defined(pruned_inds[mt_rel.col_key]))
            mt_unrel = mt_unrel.filter_cols(hl.is_defined(pruned_inds[mt_unrel.col_key]))

            # Removing sites
            window = '1e6' if pop != 'EUR' else '1e7'
            pruned_ht = hl.read_table(get_ukb_grm_pruned_ht_path(pop, window))
            mt_unrel = mt_unrel.filter_rows(hl.is_defined(pruned_ht[mt_unrel.row_key]))

            mt_unrel = mt_unrel.repartition(500).checkpoint(hl.utils.new_temp_file())

            pop = pop if window == '1e6' else f'{pop}_{window}'
            run_pca(mt_unrel, get_relatedness_path(pop, unrelated=True, extension='') + '.', args.overwrite)
            pca_loadings = hl.read_table(get_relatedness_path(pop, unrelated=True, extension='loadings.ht'))
            ht = project_individuals(pca_loadings, mt_rel)
            ht.write(get_relatedness_path(pop, extension='scores_projected.ht'), args.overwrite)
            hl.read_table(get_relatedness_path(pop, extension='scores_projected.ht')).export(
                get_relatedness_path(pop, extension='scores_projected.txt.bgz'))

    if args.generate_covariates:
        hts = []
        for pop in POPS:
            pop_path = pop if pop != 'EUR' else f'EUR_1e7'
            ht = hl.read_table(get_relatedness_path(pop_path, extension='scores_projected.ht'))
            hts.append(ht.annotate(pop=pop, related=True))
            ht = hl.read_table(get_relatedness_path(pop_path, True, extension='scores.ht'))
            ht = ht.transmute(**{f'PC{i}': ht.scores[i - 1] for i in range(1, 21)})
            hts.append(ht.annotate(pop=pop, related=False))

        ht = hts[0].union(*hts[1:])
        cov_ht = hl.import_table(get_age_sex_tsv_path(), impute=True, force=True, quote='"', key='userId').select('age', 'sex')
        cov_ht = cov_ht.annotate(age_sex=cov_ht.age * cov_ht.sex,
                                 age2=hl.int32(cov_ht.age ** 2),
                                 age2_sex=hl.int32(cov_ht.age ** 2) * cov_ht.sex)
        ht = ht.annotate(**cov_ht.key_by(userId=hl.str(cov_ht.userId))[ht.key])
        ht.write(get_covariates_ht_path(), args.overwrite)

    get_filtered_mt(imputed=False).cols().export(get_final_sample_set())

    #
    # if args.ukbb_pop_noref:
    #     """
    #     Compute PCA in UKBB population (unrelateds), without reference individuals
    #         Denser SNP set for more precise PC calculation
    #         These will be used as covariates
    #     """
    #     pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--load_ref', action='store_true')
    parser.add_argument('--dirname', default='gs://ukb-diverse-pops/reference_panels/')
    parser.add_argument('--basename', default='HGDP_1kG_maf005_geno05')
    parser.add_argument('--out_prefix', default='gs://ukb-diverse-pops/pca/')
    parser.add_argument('--load_ukbb', action='store_true')
    parser.add_argument('--intersect_ref', action='store_true')
    parser.add_argument('--pca_project', action='store_true')
    parser.add_argument('--ukb_prune_pca_project', action='store_true')

    parser.add_argument('--continental_pca', action='store_true')
    parser.add_argument('--ukbb_pop_pca', action='store_true')
    parser.add_argument('--generate_covariates', action='store_true')

    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    main(args)
