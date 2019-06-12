from gnomad_hail import *
import hail as hl
import argparse
import pandas as pd


def run_pca(my_data, out_prefix):
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(my_data.GT, k=20, compute_loadings=True)
    pca_mt = my_data.annotate_rows(pca_af=hl.agg.mean(my_data.GT.n_alt_alleles()) / 2)
    pca_loadings = pca_loadings.annotate(pca_af=pca_mt.rows()[pca_loadings.key].pca_af)

    pca_scores.write(out_prefix + 'scores.ht', args.overwrite)
    pca_scores = hl.read_table(out_prefix + 'scores.ht')
    pca_scores = pca_scores.transmute(**{f'PC{i}': pca_scores.scores[i - 1] for i in range(1, 21)})
    pca_scores.export(out_prefix + 'scores.txt.bgz')  # individual-level PCs

    pca_loadings.write(out_prefix + 'loadings.ht', args.overwrite)  # PCA loadings


def project_individuals(ref_pcs, project_prefix):
    pca_scores = hl.read_table(ref_pcs + 'scores.ht')
    pca_loadings = hl.read_table(ref_pcs + 'loadings.ht')
    mt = hl.read_matrix_table(project_prefix + '.mt')
    print(f'Projecting population PCs')
    mt_projections = pc_project(mt, pca_loadings)
    mt_projections = mt_projections.transmute(**{f'PC{i}': mt_projections.scores[i - 1] for i in range(1, 21)})
    mt_projections.export(project_prefix + '_scores.txt.bgz')
    # hl.read_table('gs://armartin/pigmentation/pca/ukbb_ref_scores.ht').export(
    #     'gs://armartin/pigmentation/pca/ukbb_ref_scores.txt.bgz')


def main(args):
    if args.load_ref:
        ref = hl.import_plink(bed='gs://ukb-diverse-pops/pca/data/' + args.pop + 'HGDP_1kG_maf005_geno05.bed',
                              bim='gs://ukb-diverse-pops/pca/data/' + args.pop + 'HGDP_1kG_maf005_geno05.bim',
                              fam='gs://ukb-diverse-pops/pca/data/' + args.pop + 'HGDP_1kG_maf005_geno05.fam',
                              min_partitions=100)
        ref.describe()

        print('sites in ref data: ' + str(ref.count()))  # (639590, 3547)
        ref.write('gs://ukb-diverse-pops/pca/data/HGDP_1kG_maf005_geno05.mt', args.overwrite)

    if args.load_ukbb:
        ref = hl.read_matrix_table('gs://ukb-diverse-pops/pca/data/' + args.pop + 'HGDP_1kG_maf005_geno05.mt')
        samples = hl.read_table('gs://armartin/pigmentation/pigmentation_phenos_covs_pops.ht')
        ukbb = hl.read_matrix_table('gs://phenotype_31063/hail/genotype/ukb31063.genotype.mt')
        ukbb = ukbb.annotate_cols(**samples[ukbb.s])

        # filter ukbb to sites in ref & array data
        ukbb_in_ref = ukbb.filter_rows(hl.is_defined(ref.rows()[ukbb.row_key]))
        print('sites, inds in ref and UKBB data: ' + str(ukbb_in_ref.count()))  # (64233, 488377)

        ukbb_in_ref.write('gs://ukb-diverse-pops/pca/data/' + args.pop + 'ukbb_globalref.mt', args.overwrite)

        # filter ref to ukbb sites
        ref_in_ukbb = ref.filter_rows(hl.is_defined(ukbb.rows()[ref.row_key]))
        print('sites, inds in ref and UKBB data: ' + str(ref_in_ukbb.count()))  # (64233, 3547)
        ref_in_ukbb.write('gs://ukb-diverse-pops/pca/data/globalref_ukbb_intersect.mt', args.overwrite)

        # filter ukbb to unrel individuals
        # ukbb_in_ref_unrel = ukbb_in_ref.filter_cols(ukbb_in_ref.covariates['used_in_pca_calculation'])
    if args.global_pca:
        """
        Compute PCA in global reference panel, project UKBB individuals into PCA space
        """
        ref_in_ukbb = hl.read_matrix_table('gs://ukb-diverse-pops/pca/data/globalref_ukbb_intersect.mt')
        print('Computing reference PCs')
        run_pca(ref_in_ukbb, 'gs://ukb-diverse-pops/pca/data/globalref_ukbb_')

        # project ukbb
        project_individuals('gs://ukb-diverse-pops/pca/data/globalref_ukbb_',
                            'gs://ukb-diverse-pops/pca/data/ukbb_globalref')

    if args.continental_pca:
        """
        Compute PCA within reference panel super pops, project UKBB individuals into PCA space
        1. Filter UKBB to individuals in continental population
        2. Run PCA on continental ref
        3. Project UKBB inds
        """
        pass

    if args.ukbb_pop_pca:
        """
        Compute PCA in each UKBB population (unrelateds), project reference individuals and relateds into PCA space
        1. Filter UKBB to individuals in continental population
        2. Run PC-relate on these individuals
        3. Filter UKBB population to unrelated individuals
        4. Run PCA on UKBB population 
        5. Project reference panel 
        """
        pass

    if args.ukbb_pop_noref:
        """
        Compute PCA in UKBB population (unrelateds), without reference individuals
            Denser SNP set for more precise PC calculation
            These will be used as covariates
        """
        pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--load_ref', action='store_true')
    parser.add_argument('--load_ukbb', action='store_true')
    parser.add_argument('--global_pca', action='store_true')
    parser.add_argument('--continental_pca', action='store_true')
    parser.add_argument('--pop', default='', help='AFR, AMR, CSA, EAS, EUR, MID')
    parser.add_argument('--ukbb_pop_pca', action='store_true')
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    main(args)
