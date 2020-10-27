#!/usr/bin/env python3

__author__ = 'konradk'

import argparse
from hail.linalg import BlockMatrix
# from gnomad.utils import slack
from ukbb_pan_ancestry import *
import math

def get_final_sumstats_mt_for_export():
    mt0 = load_final_sumstats_mt(filter_sumstats=False,
                                 filter_variants=False,
                                 separate_columns_by_pop=False,
                                 annotate_with_nearest_gene=False)
    mt0 = mt0.select_rows()
    return mt0

def resume_write_block_matrices(bms, path_prefix, start, stop):
    r'''
    Used to write the subset of failed intermediate multiply block matrices `bms`
    '''
    # TODO: hl.hadoop_is_file to detect if bms have _SUCCESS file to automatically
    # determine missing bms?
    for i, b in enumerate(bms[start:stop], start):
        b.write(f'{path_prefix}_{i}')

# n_block_cols = 7078
# mul_splits = 200
# inner_brange_size = int(math.ceil(n_block_cols / mul_splits))
# print(f'inner_brange_size: {inner_brange_size}')
# # split_points = list(range(0, n_block_cols, inner_brange_size)) + [n_block_cols]
# split_points = list(range(0, n_block_cols, inner_brange_size)) + [n_block_cols]
# print(len(split_points))
# inner_ranges = list(zip(split_points[:-1], split_points[1:]))
# print(f'len(inner_ranges): {len(inner_ranges)}')
            
def tree_matmul_tree_matsum(bm1, bm2, mul_splits: int, sum_splits: int = None,
                            path_prefix: str = None, read_if_exists=False):
    r'''
    Version of tree_matmul() that allows for intermediate sums of matrix 
    multiplication. `sum_splits` must be a divisor of `mul_splits`
    '''
    # TODO: Make a private function that acts recursively to ensure that the 
    # matrix sums never include more than a maximum number of matrices
    # assert mul_splits%sum_splits==0, '`sum_splits` must be a divisor of `mul_splits'

    if not read_if_exists:
        print(bm1._n_block_cols)
        print(mul_splits)
        inner_brange_size = int(math.ceil(bm1._n_block_cols / mul_splits))
        print(f'inner_brange_size: {inner_brange_size}')
        split_points = list(range(0, bm1._n_block_cols, inner_brange_size)) + [bm1._n_block_cols]
        print(split_points)
        inner_ranges = list(zip(split_points[:-1], split_points[1:]))
        print(f'len(inner_ranges): {len(inner_ranges)}')
        blocks_to_multiply = [(bm1._select_blocks((0, bm1._n_block_rows), (start, stop)),
                                bm2._select_blocks((start, stop), (0, bm2._n_block_cols))) for start, stop in inner_ranges]
    
        intermediate_multiply_exprs = [b1 @ b2 for b1, b2 in blocks_to_multiply]
        print(len(intermediate_multiply_exprs))
        print(f'Writing {mul_splits} intermediate matrices to {path_prefix}')
        hl.experimental.write_block_matrices(intermediate_multiply_exprs, path_prefix)
    
    read_intermediates = [BlockMatrix.read(f"{path_prefix}_{i}") for i in range(0, mul_splits)]
    
    tracked_partial_sums = []

    sum_block_size = math.ceil(mul_splits/sum_splits)
    for i in range(sum_splits): #range(sum_splits):
        partial_sum_path = f"{path_prefix}-partial-{i}"
        # sum(read_intermediates[i*sum_block_size:(i+1)*sum_block_size]).write(partial_sum_path)
        tracked_partial_sums.append(BlockMatrix.read(partial_sum_path))
        
    return sum(tracked_partial_sums)
    

def main(args):
    hl.init(default_reference='GRCh37', log='/prs.log',
            spark_conf={'spark.hadoop.fs.gs.requester.pays.mode': 'AUTO', 'spark.hadoop.fs.gs.requester.pays.project.id': 'ukbb-diversepops-neale'})

    if args.prepare_sumstats_matrix:
        # get meta mt and separate by pop combo
        meta_mt = hl.read_matrix_table(get_meta_analysis_results_path())
        meta_mt = separate_results_mt_by_pop(meta_mt, 'meta_analysis_data', 'meta_analysis')
        meta_mt = meta_mt.annotate_cols(clump_pops=meta_mt.meta_analysis_data.pop)
        meta_mt = meta_mt.key_cols_by('clump_pops', *meta_mt.col_key)
        
        # get sumstats mt and separate by pop combo
        ss_mt = get_final_sumstats_mt_for_export()
        ss_mt = separate_results_mt_by_pop(ss_mt, 'pheno_data', 'summary_stats')
        ss_mt = ss_mt.annotate_cols(clump_pops=hl.array([ss_mt.pheno_data.pop]))
        ss_mt = ss_mt.key_cols_by(*meta_mt.col_key)
        
        # join meta results and sumstats mt
        # NOTE: union_cols() requires the same entry fields schema
        meta_mt = meta_mt.select_entries(BETA = meta_mt.meta_analysis.BETA,
                                         Pvalue = meta_mt.meta_analysis.Pvalue).select_cols().select_rows()
        ss_mt = ss_mt.select_entries(BETA = ss_mt.summary_stats.BETA,
                                     Pvalue = ss_mt.summary_stats.Pvalue).select_cols().select_rows()
        mt = meta_mt.union_cols(ss_mt)
        
        # filter to distinct cols
        # NOTE: distinct_by_col() does not allow a col key of type `list`
        mt = mt.annotate_cols(clump_pops_str = hl.delimit(mt.clump_pops)).key_cols_by('clump_pops_str', *[k for k in mt.col_key if k!='clump_pops']).distinct_by_col()
        mt = mt.distinct_by_col()
        
        # ensure that betas are not missing
        ss_mt = ss_mt.annotate_cols(clump_pops_str = hl.delimit(ss_mt.clump_pops)).key_cols_by('clump_pops_str', *[k for k in ss_mt.col_key if k!='clump_pops'])
        mt = mt.annotate_entries(BETA = hl.or_else(mt.BETA, ss_mt[mt.row_key, mt.col_key].BETA),
                                 Pvalue = hl.or_else(mt.Pvalue, ss_mt[mt.row_key, mt.col_key].Pvalue))
        
        # read clump mt and separate by pop combo
        clump_mt = hl.read_matrix_table(get_clumping_results_path(high_quality=args.high_quality, 
                                                                  max_pops=args.max_pops))
        if args.max_pops:
            # if max_pops=True, the clump_mt is already separated by pop
            # these steps are necessary to make downstream code usable for both max_pops=True/False
            clump_mt = clump_mt.annotate_entries(plink_clump = hl.struct(TOTAL = clump_mt.TOTAL))
            clump_mt = clump_mt.annotate_cols(pop_index = 0)
        else:
            clump_mt = separate_results_mt_by_pop(clump_mt, 'clump_pops', 'plink_clump', skip_drop=True)
        
        clump_mt = clump_mt.annotate_cols(clump_pops_str = hl.delimit(clump_mt.clump_pops))
        clump_mt = clump_mt.drop('clump_pops').key_cols_by(*mt.col_key)
        
        # join sumstats/meta-analysis with clump mt
        mt = all_axis_join(mt, clump_mt)
        
        mt = mt.filter_cols(hl.is_defined(mt.pop_index))
        
        print(f'\n\nMatrix dimensions (before explode by p-threshold): {mt.count()}\n')
        mt = explode_by_p_threshold(mt).unfilter_entries()
        # Write pheno data for later use
        mt.add_col_index('idx').key_cols_by('idx').cols().write(
            get_clump_sumstats_col_ht_path(high_quality=args.high_quality,
                                           max_pops=args.max_pops), 
            args.overwrite)
        BlockMatrix.write_from_entry_expr(
            hl.or_else(mt.BETA * hl.is_defined(mt.plink_clump.TOTAL) * hl.int(mt.Pvalue < mt.p_threshold), 0.0),
            get_clump_sumstats_bm_path(high_quality=args.high_quality,
                                        max_pops=args.max_pops), 
            args.overwrite)
        # 2020-06-25 01:49:32 Hail: INFO: Wrote all 7078 blocks of 28987534 x 3530 matrix with block size 4096.
        # If clump_mt is significantly smaller than meta_mt, consider putting that on the left of the join,
        # then filter the genotype matrix to only those SNPs (pilot would go from 28.9M -> 21.2M)

    if args.prepare_genotype_matrix:
        meta_mt = hl.read_matrix_table(get_meta_analysis_results_path())
        mt = get_filtered_mt_with_x()
        mt = mt.filter_rows(hl.is_defined(meta_mt.rows()[mt.row_key]))
        # Write sample data for later use
        mt = mt.key_cols_by(userId=hl.int32(mt.s))
        mt.cols().add_index().write(genotype_samples_ht_path, args.overwrite)
        BlockMatrix.write_from_entry_expr(mt.dosage, genotype_bm_path, args.overwrite)
        # 2020-06-25 19:18:14 Hail: INFO: Wrote all 764424 blocks of 28987534 x 441345 matrix with block size 4096.

    if args.compute_prs:
        sumstats_bm = BlockMatrix.read(get_clump_sumstats_bm_path(high_quality=args.high_quality, 
                                                                  max_pops=args.max_pops))
        genotype_bm = BlockMatrix.read(genotype_bm_path)
        mul_splits = 197 # sumstats_bm.shape[1]//10000*10
        sum_splits = 20 #int(mul_splits/10)
        assert mul_splits>10 # if not more than 10, sum_splits is not necessary
        prs_bm = tree_matmul_tree_matsum(genotype_bm.T, sumstats_bm, mul_splits=mul_splits, 
                                         sum_splits=sum_splits, path_prefix = f'{temp_bucket}/prs/tree_matmul{"_max_pops" if args.max_pops else ""}',
                                         read_if_exists = True)
        prs_bm.write(get_prs_bm_path(high_quality=args.high_quality,
                                     max_pops=args.max_pops), args.overwrite)

    if args.create_prs_mt:
        prs_bm = BlockMatrix.read(get_prs_bm_path(high_quality=args.high_quality,
                                                  max_pops=args.max_pops))
        pheno_ht = hl.read_table(get_clump_sumstats_col_ht_path(high_quality=args.high_quality,
                                                                max_pops=args.max_pops)).key_by('idx')
        samples_ht = hl.read_table(genotype_samples_ht_path).key_by('idx')
        # 10k partitions for 370 GB table (441k x 108k) = 37 MB/partition
        # 5014 partitions for 240 GB table (441k x 72k) = 48 MB/partition (max_pops)
        n_partitions = 15000 #int(1000*(pheno_ht.count()/72*5)//1000) # or hard code
        mt = BlockMatrix.to_matrix_table_row_major(prs_bm, n_partitions=n_partitions).rename({'element': 'score'}) 
        mt = mt.annotate_cols(**pheno_ht[mt.col_key]).key_cols_by(*PHENO_KEY_FIELDS)
        mt = mt.annotate_rows(**samples_ht[mt.row_key]).key_rows_by('userId')
        mt.write(get_prs_mt_path(high_quality=args.high_quality, 
                                 max_pops=args.max_pops), 
                 args.overwrite)

    if args.assess_prs:
        prs_mt = hl.read_matrix_table(get_prs_mt_path(high_quality=args.high_quality, 
                                                      max_pops=args.max_pops))
        pheno_mt = get_ukb_pheno_mt()  # TODO: fix all phenos to new keying scheme
        pheno_mt = pheno_mt.key_cols_by(
            **pheno_mt.col_key.annotate(modifier=hl.if_else(pheno_mt.trait_type == "biomarkers", "irnt", pheno_mt.modifier)))
        mt = prs_mt.annotate_entries(**pheno_mt[prs_mt.row_key, prs_mt.col_key])
        mt = mt.annotate_cols(description = pheno_mt.cols()[mt.col_key].description)
        for pop in POPS:
            mt_pop = mt.filter_rows(mt.pop==pop)
            mt_pop = mt_pop.annotate_cols(prs_corr=hl.agg.linreg(mt_pop.both_sexes, [1.0, mt_pop.score]))
            cols = mt_pop.cols()
            cols.select('description', 
                        'p_threshold',
                        clump_pops_str=hl.delimit(cols.clump_pops,'-'),
                        prs_corr_r2=cols.prs_corr.multiple_r_squared, 
                        prs_corr_pval=cols.prs_corr.p_value[1], 
                        prs_corr_n=cols.prs_corr.n).export(f'gs://ukbb-diverse-temp-30day/prs/assess_prs{"_max_pops" if args.max_pops else ""}.{pop}.tsv.gz')
        # mt = mt.filter_cols(mt.description == 'Type 2 diabetes')
        # import gnomad.utils.file_utils
        # gnomad.utils.file_utils.select_primitives_from_ht(mt.entries()).export(f'{temp_bucket}/prs/prs_test.txt.bgz')
        # mt = mt.annotate_cols(prs_corr=hl.agg.linreg(mt.both_sexes, [1.0, mt.score]))
        # ht = mt.cols().checkpoint(get_prs_assess_ht_path(args.high_quality), args.overwrite, _read_if_exists=not args.overwrite)
        # print(ht.aggregate(hl.agg.fraction(ht.prs_corr.p_value[1] < 0.05)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--prepare_sumstats_matrix', help='Prepare summary stats blockmatrix', action='store_true')
    parser.add_argument('--prepare_genotype_matrix', help='Prepare genotype blockmatrix', action='store_true')
    parser.add_argument('--compute_prs', help='Compute PRS', action='store_true')
    parser.add_argument('--create_prs_mt', help='Convert PRS blockmatrix to MT', action='store_true')
    parser.add_argument('--assess_prs', help='Assess PRS performance', action='store_true')
    parser.add_argument('--high_quality', help='Overwrite everything', action='store_true')
    parser.add_argument('--max_pops', help='Whether to use "max_pops" results', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    # if args.slack_channel:
    #     from slack_token_pkg.slack_creds import slack_token
    #     with slack.slack_notifications(slack_token, args.slack_channel):
    main(args)
