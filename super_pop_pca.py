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


def run_pca(my_data, out_prefix):
    """
    Run PCA on a dataset
    :param my_data: dataset to run PCA on
    :param out_prefix: directory and filename prefix for where to put PCA output
    :return:
    """
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(my_data.GT, k=20, compute_loadings=True)
    pca_mt = my_data.annotate_rows(pca_af=hl.agg.mean(my_data.GT.n_alt_alleles()) / 2)
    pca_loadings = pca_loadings.annotate(pca_af=pca_mt.rows()[pca_loadings.key].pca_af)

    pca_scores.write(out_prefix + 'scores.ht', args.overwrite)
    pca_scores = hl.read_table(out_prefix + 'scores.ht')
    pca_scores = pca_scores.transmute(**{f'PC{i}': pca_scores.scores[i - 1] for i in range(1, 21)})
    pca_scores.export(out_prefix + 'scores.txt.bgz')  # individual-level PCs

    pca_loadings.write(out_prefix + 'loadings.ht', args.overwrite)  # PCA loadings


def project_individuals(ref_pcs, project_data, project_prefix):
    """
    Project samples into predefined PCA space
    :param ref_pcs: existing PCA space
    :param project_data: matrix table of data to project
    :param project_prefix: directory and filename prefix for where to put PCA projection output
    :return:
    """
    pca_scores = hl.read_table(ref_pcs + 'scores.ht')
    pca_loadings = hl.read_table(ref_pcs + 'loadings.ht')

    mt = hl.read_matrix_table(project_data)
    print(f'Projecting population PCs')
    mt_projections = pc_project(mt, pca_loadings)
    mt_projections = mt_projections.transmute(**{f'PC{i}': mt_projections.scores[i - 1] for i in range(1, 21)})
    mt_projections.export(project_prefix + '_scores.txt.bgz')
    # hl.read_table('gs://armartin/pigmentation/pca/ukbb_ref_scores.ht').export(
    #     'gs://armartin/pigmentation/pca/ukbb_ref_scores.txt.bgz')


def main(args):
    if args.load_ref:
        load_ref(args.dirname, args.basename)

    if args.load_ukbb:
        ##
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
        project_individuals(args.out_prefix + args.basename + '_ukbb_',
                            args.dirname + 'intersect_ukbb_' + args.basename + '.mt',
                            args.out_prefix + 'ukbb_' + args.basename)

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
        3. Filter UKBB population to unrelated individuals
        4. Run PCA on UKBB unrelateds within population
        5. Project relateds
        """

        pop_assign = hl.import_table(get_ukb_meta_pop_tsv_path()).drop(*map(lambda x: 'PC' + str(x), range(1,21)))
        pops = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']

        #for pop in ['AFR', 'CSA']:
        for pop in ['EUR']:
            mt = hl.read_matrix_table(get_ukb_grm_mt_path(pop))
            pruned_ht = hl.read_table(get_ukb_grm_pruned_ht_path(pop))
            mt = mt.filter_rows(hl.is_defined(pruned_ht[mt.row_key]))

            # run PC-relate
            relatedness_ht = hl.pc_relate(mt.GT, min_individual_maf=0.05,
                                          k=10,  # block size = 512 for smaller, 4096 for EUR pops
                                          block_size=4096, min_kinship=0.05, statistics='kin2')
            relatedness_ht.write('gs://ukb-diverse-pops/pca/relatedness/' + pop + '_rel.ht', args.overwrite)
            relatedness_ht = hl.read_table('gs://ukb-diverse-pops/pca/relatedness/' + pop + '_rel.ht')

            # identify individuals in pairs to remove
            related_samples_to_remove = hl.maximal_independent_set(relatedness_ht.i, relatedness_ht.j, False)
            mt_unrel = mt.filter_cols(
                hl.is_defined(related_samples_to_remove[mt.col_key]), keep = False)
            mt_rel = mt.filter_cols(
                hl.is_defined(related_samples_to_remove[mt.col_key]), keep=True)
            mt_unrel.write('gs://ukb-diverse-pops/pca/relatedness/' + pop + '_unrel.mt', args.overwrite)
            mt_rel.write('gs://ukb-diverse-pops/pca/relatedness/' + pop + '_rel.mt', args.overwrite)

            mt_rel = hl.read_matrix_table('gs://ukb-diverse-pops/pca/relatedness/' + pop + '_rel.mt')
            mt_unrel = hl.read_matrix_table('gs://ukb-diverse-pops/pca/relatedness/' + pop + '_unrel.mt')

            run_pca(mt_unrel, args.out_prefix + pop + '_')
            # ref_pcs, project_data, project_prefix
            project_individuals(args.out_prefix + pop + '_',
                                'gs://ukb-diverse-pops/pca/relatedness/' + pop + '_rel.mt',
                                'gs://ukb-diverse-pops/pca/relatedness/' + pop + '_rel')

            mt_unrel.count()

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

    parser.add_argument('--continental_pca', action='store_true')
    parser.add_argument('--ukbb_pop_pca', action='store_true')

    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    main(args)
