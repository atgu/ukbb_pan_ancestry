from ukbb_pan_ancestry.resources.results import *
from ukbb_pan_ancestry.heritability.import_heritability import get_h2_ht
P_THRESHOLDS = {'s1': 5e-8, 's2': 1e-6, 's3': 1e-4, 's4': 1e-3, 's5': 1e-2, 's6': .05, 's7': .1, 's8': .2, 's9': .5, 's10': 1.}


def annotate_nearest_gene(t, add_contig: bool = False, add_only_gene_symbols_as_str: bool = False, loc: str = 'nearest_genes'):
    intervals_ht = hl.read_table(get_gene_intervals_path())
    if add_contig:
        intervals_ht = intervals_ht.annotate(contig=intervals_ht.interval.start.contig)
    annotation = intervals_ht.index(t.locus, all_matches=True)
    if add_only_gene_symbols_as_str:
        annotation = hl.delimit(annotation.gene_name)
    if loc: annotation = {loc: annotation}
    return t.annotate_rows(**annotation) if isinstance(t, hl.MatrixTable) else t.annotate(**annotation)


def filter_lambda_gc(lambda_gc):
    return (lambda_gc > 0.5) & (lambda_gc < 2)


def load_final_sumstats_mt(filter_phenos: bool = True, filter_variants: bool = True,
                           filter_sumstats: bool = True, separate_columns_by_pop: bool = True,
                           annotate_with_nearest_gene: bool = True, add_only_gene_symbols_as_str: bool = False,
                           load_contig: str = None, filter_pheno_h2_qc: bool = True):
    mt = hl.read_matrix_table(get_variant_results_path('full', 'mt')).drop('gene', 'annotation')
    if load_contig:
        mt = mt.filter_rows(mt.locus.contig == load_contig)
    variant_qual_ht = hl.read_table(get_variant_results_qc_path())
    mt = mt.annotate_rows(**variant_qual_ht[mt.row_key])
    pheno_qual_ht = hl.read_table(get_analysis_data_path('lambda', 'lambdas', 'full', 'ht'))
    mt = mt.annotate_cols(**pheno_qual_ht[mt.col_key])
    h2_qc_ht = hl.read_table(get_h2_ht())
    mt = mt.annotate_cols(heritability = h2_qc_ht[mt.col_key].heritability)


    def update_pheno_struct(pheno_struct, mt):
        pheno_struct = pheno_struct.annotate(saige_heritability = pheno_struct.heritability,
                                             heritability = mt.paired_pop_h2.get(pheno_struct.pop, 
                                                                                 hl.missing(hl.tstruct(**mt.heritability[0].dtype))))
        pheno_struct = pheno_struct.annotate(heritability = pheno_struct.heritability.drop('pop'))
        return pheno_struct


    mt = mt.annotate_cols(paired_pop_h2 = hl.dict(hl.map(lambda x: (x.pop,x), mt.heritability)))
    mt = mt.annotate_cols(pheno_data = mt.pheno_data.map(lambda x: update_pheno_struct(x, mt)))
    mt = mt.drop('heritability', 'paired_pop_h2')

    if filter_phenos:
        keep_phenos = hl.enumerate(mt.pheno_data).filter(
            lambda x: filter_lambda_gc(x[1].lambda_gc))

        mt = mt.annotate_cols(
            pheno_indices=keep_phenos.map(lambda x: x[0]),
            pheno_data=keep_phenos.map(lambda x: x[1]))
        mt = mt.annotate_entries(
            summary_stats=hl.zip_with_index(mt.summary_stats).filter(
                lambda x: mt.pheno_indices.contains(x[0])).map(lambda x: x[1])
        ).drop('pheno_indices')
        mt = mt.filter_cols(hl.len(mt.pheno_data) > 0)

    if filter_pheno_h2_qc:
        mt = mt.filter_cols(mt.pheno_data.any(lambda x: (hl.is_defined(x.heritability.qcflags.pass_all)) & \
                                              (x.heritability.qcflags.pass_all)))
        # filter arrays
        def filter_pop_array(mt, fields):
            expr_dict = {}
            for x in fields:
                this_expr = hl.zip(mt[x], mt.tf_filt_h2).filter(lambda x: x[1]).map(lambda x: x[0])
                expr_dict.update({x: this_expr})
            return expr_dict
        
        mt = mt.annotate_cols(tf_filt_h2 = mt.pheno_data.map(lambda x: x.heritability.qcflags.pass_all))
        # column arrays
        mt = mt.annotate_cols(**filter_pop_array(mt, ['pheno_data']))
        # entry arrays
        mt = mt.annotate_entries(**filter_pop_array(mt, ['summary_stats']))
        mt = mt.drop('tf_filt_h2')

    if filter_sumstats:
        mt = mt.annotate_entries(summary_stats=mt.summary_stats.map(
            lambda x: hl.or_missing(~x.low_confidence, x)
        ))
        mt = mt.filter_entries(~mt.summary_stats.all(lambda x: hl.is_missing(x.Pvalue)))

    if filter_variants:
        mt = mt.filter_rows(mt.high_quality)

    if annotate_with_nearest_gene:
        mt = annotate_nearest_gene(mt, add_only_gene_symbols_as_str=add_only_gene_symbols_as_str)

    if separate_columns_by_pop:
        mt = separate_results_mt_by_pop(mt)

    return mt


def separate_results_mt_by_pop(mt, col_field = 'pheno_data', entry_field = 'summary_stats', skip_drop: bool = False):
    mt = mt.annotate_cols(col_array=hl.enumerate(mt[col_field])).explode_cols('col_array')
    mt = mt.transmute_cols(pop_index=mt.col_array[0], **{col_field: mt.col_array[1]})
    mt = mt.annotate_entries(**{entry_field: mt[entry_field][mt.pop_index]})
    if not skip_drop:
        mt = mt.drop('pop_index')
    return mt


def explode_by_p_threshold(mt):
    mt = mt.annotate_cols(p_threshold=hl.literal(list(P_THRESHOLDS.items()))).explode_cols('p_threshold')
    mt = mt.transmute_cols(p_threshold_name=mt.p_threshold[0], p_threshold=mt.p_threshold[1])
    return mt

