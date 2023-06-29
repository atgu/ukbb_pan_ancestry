#!/usr/bin/env python3
# coding: utf-8

__author__ = "mkanai"

import argparse
import atexit
from hail.typecheck import check
from hail.utils.hadoop_utils import hadoop_exists
import numpy as np
import hail as hl
import uuid
from hail.linalg import BlockMatrix
from hail.utils import new_temp_file, hadoop_open, timestamp_path
from ukbb_pan_ancestry.resources import get_filtered_mt, get_variant_results_qc_path, POPS, temp_bucket_7day
from ukbb_pan_ancestry.resources.ld import (
    get_hm3_snplist_path,
    get_ld_score_flat_file_path,
    get_ld_variant_index_path,
    get_ld_matrix_path,
    get_ld_score_ht_path,
)


def new_gs_temp_path():
    return f"{temp_bucket_7day}/{str(uuid.uuid4())}"


def checkpoint_tmp(hail_obj, path=None, overwrite=False, force_row_major=True):
    if path is None:
        path = new_gs_temp_path()

    if isinstance(hail_obj, BlockMatrix):
        return hail_obj.checkpoint(path, overwrite=overwrite, force_row_major=force_row_major)
    else:
        return hail_obj.checkpoint(path, overwrite=overwrite)


def normalize_bm(bm):
    n = bm.shape[1]
    m1 = checkpoint_tmp(bm.sum(axis=1))
    m2 = checkpoint_tmp((bm ** 2).sum(axis=1))
    mean = m1 / n
    # biased is n; unbiased is n - 1
    stdev = ((m2 - m1 ** 2 / n) / n).sqrt()
    # add a min float value to prevent zero division due to machine precision; 3.35e-14
    bm_norm = (bm - mean) / (stdev + 1.18e-38)
    return bm_norm


# cf. https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/0.2/create_ldsc_hm3_table.py
def write_ldsc_hm3_snplist(info_threshold=0.9, maf_threshold=0.01, overwrite=False):
    # Filter variants
    ht = hl.read_table(get_variant_results_qc_path())
    # in autosomes
    ht = ht.filter(ht.locus.in_autosome())
    # no MHC
    ht = ht.filter(~hl.parse_locus_interval("6:28477797-33448354").contains(ht.locus))
    # info > 0.9
    ht = ht.filter(ht.info > info_threshold)
    # SNP only
    ht = ht.filter(hl.is_snp(ht.alleles[0], ht.alleles[1]))
    # no multi-allelic sites
    loc_count = ht.group_by(ht.locus).aggregate(nloc=hl.agg.count())
    loc_count = loc_count.filter(loc_count.nloc > 1)
    multi_sites = loc_count.aggregate(hl.agg.collect_as_set(loc_count.locus), _localize=False)
    ht = ht.filter(~multi_sites.contains(ht.locus))

    # in HM3
    hm3_snps = hl.read_table("gs://ukbb-ldsc-dev/ukb_hm3_snplist/hm3.r3.b37.auto_bi_af.ht")
    hm3_snps = hm3_snps.select()
    ht = ht.join(hm3_snps, "right")
    # no strand ambiguity
    ht = ht.filter(~hl.is_strand_ambiguous(ht.alleles[0], ht.alleles[1]))

    ht = checkpoint_tmp(ht)

    def get_maf(af):
        return 0.5 - hl.abs(0.5 - af)

    # MAF > 1% in UKB & gnomad genome/exome (if defined) for each population
    for pop in POPS:
        snplist = ht.filter(
            hl.rbind(
                ht.freq[ht.freq.index(lambda x: x.pop == pop)],
                lambda y: (get_maf(y.af) > maf_threshold)
                & (hl.is_missing(y.gnomad_genomes_af) | (get_maf(y.gnomad_genomes_af) > maf_threshold))
                & (hl.is_missing(y.gnomad_exomes_af) | (get_maf(y.gnomad_exomes_af) > maf_threshold)),
            )
        )
        snplist = snplist.select("rsid")
        snplist.write(get_hm3_snplist_path(pop), overwrite=overwrite)


# cf: https://github.com/broadinstitute/gnomad_qc/blob/master/gnomad_qc/v2/annotations/generate_ld_data.py
def copmute_ldscore(ht, r2_adj, out_name, overwrite):
    # Note that the original ld matrix is triangular
    l2row = checkpoint_tmp(r2_adj.sum(axis=0)).T
    l2col = checkpoint_tmp(r2_adj.sum(axis=1))
    r2_diag = checkpoint_tmp(r2_adj.diagonal()).T
    l2 = l2row + l2col - r2_diag
    l2_bm_tmp = new_temp_file()
    l2_tsv_tmp = new_gs_temp_path()

    l2.write(l2_bm_tmp, force_row_major=True)
    BlockMatrix.export(l2_bm_tmp, l2_tsv_tmp)

    ht_scores = hl.import_table(l2_tsv_tmp, no_header=True, impute=True)
    ht_scores = ht_scores.add_index().rename({"f0": "ld_score"})
    ht_scores = ht_scores.key_by("idx")
    ht = ht.add_index()
    ht = ht.annotate(**ht_scores[ht.idx]).drop("idx")
    ht = ht.checkpoint(out_name, overwrite)
    return ht


def copmute_statified_ldscore(ht, r2_adj, out_name, overwrite, ht_annot=None):
    if ht_annot is not None:
        # Merge with LD matrix indices
        ht_annot = ht.select().join(ht_annot, "left")
        # Replace NA with zero (assuming binary annotations)
        cols = list(ht_annot.row.keys() - ht_annot.key.keys())
        ht_annot = ht_annot.annotate(
            **{key: hl.if_else(hl.is_defined(ht_annot[key]), ht_annot[key], 0) for key in cols}
        )
    else:
        ht_annot = ht.select()

    # Add a base annotation (includes all variants)
    ht_annot = ht_annot.annotate(base=1)
    # Create an annotation BlockMatrix
    cols = sorted(list(ht_annot.row.keys() - ht_annot.key.keys()))
    mt_annot = ht_annot.to_matrix_table_row_major(columns=cols, entry_field_name="value", col_field_name="annotation")
    tmp_bm_path = new_gs_temp_path()
    BlockMatrix.write_from_entry_expr(mt_annot.value, tmp_bm_path, mean_impute=False, center=False, normalize=False)
    bm_annot = BlockMatrix.read(tmp_bm_path)

    # Note that the original ld matrix is triangular
    l2row = r2_adj @ bm_annot
    l2col = r2_adj.T @ bm_annot
    l2_diag = r2_adj.diagonal().T * bm_annot
    l2 = l2row + l2col - l2_diag
    l2_bm_tmp = new_gs_temp_path()
    l2_tsv_tmp = new_gs_temp_path()

    l2.write(l2_bm_tmp, force_row_major=True)
    BlockMatrix.export(l2_bm_tmp, l2_tsv_tmp)

    ht_scores = hl.import_table(l2_tsv_tmp, no_header=True, impute=True)
    ht_scores = ht_scores.add_index().rename({f"f{i}": v for i, v in enumerate(cols)})
    ht_scores = ht_scores.key_by("idx")
    ht = ht.add_index()
    ht = ht.annotate(**ht_scores[ht.idx]).drop("idx")
    ht = ht.checkpoint(out_name, overwrite)
    return ht


def export_ldscore(ht, pop):
    hm3_snps = hl.read_table(get_hm3_snplist_path(pop))

    ht = ht.select(
        CHR=ht.locus.contig,
        SNP=hl.variant_str(ht.locus, ht.alleles),
        RSID=ht.rsid,
        BP=ht.locus.position,
        L2=ht.ld_score,
        MAF=0.5 - hl.abs(0.5 - ht.AF),
    )
    count = ht.aggregate(hl.struct(M=hl.agg.count(), M_5_50=hl.agg.sum(ht.MAF > 0.05)))
    ht = ht.filter(hl.is_defined(hm3_snps[ht.locus, ht.alleles]))
    ht = ht.key_by().drop("locus", "alleles", "MAF")

    with hadoop_open(get_ld_score_flat_file_path(pop, extension="M"), "w") as f:
        f.write(f"{count.M}\n")
    with hadoop_open(get_ld_score_flat_file_path(pop, extension="M_5_50"), "w") as f:
        f.write(f"{count.M_5_50}\n")

    # LD score with variant ids
    ht.drop("RSID").export(get_ld_score_flat_file_path(pop))
    # with rsids
    ht.transmute(SNP=ht.RSID).export(get_ld_score_flat_file_path(pop, rsid=True))


def export_stratified_ldscore(ht, ht_annot, pop, annot):
    hm3_snps = hl.read_table(get_hm3_snplist_path(pop))

    cols = ['base'] + sorted(list(ht_annot.row.keys() - ht_annot.key.keys()))

    ht_annot = ht.select().join(ht_annot, 'left')
    ht_annot = ht_annot.annotate(base=1)
    ht_annot = ht_annot.join(ht.annotate(MAF=0.5 - hl.abs(0.5 - ht.AF)).select("MAF"), "left")
    mt_annot = ht_annot.to_matrix_table_row_major(columns=cols, entry_field_name="value", col_field_name="annotation")
    mt_annot = mt_annot.annotate_cols(
        M=hl.agg.count_where(mt_annot.value == 1),
        M_5_50=hl.agg.count_where((mt_annot.value == 1) & (mt_annot.MAF > 0.05)),
    )
    count = mt_annot.cols().to_pandas().set_index("annotation")

    with hadoop_open(get_ld_score_flat_file_path(pop, extension="M", annot=annot), "w") as f:
        count[["M"]].T.to_csv(f, sep="\t", header=False, index=False)
    with hadoop_open(get_ld_score_flat_file_path(pop, extension="M_5_50", annot=annot), "w") as f:
        count[["M_5_50"]].T.to_csv(f, sep="\t", header=False, index=False)

    # export .annot.gz
    ht_annot = ht_annot.join(hm3_snps, "inner")
    ht_annot = ht_annot.annotate(
        CHR=ht_annot.locus.contig,
        SNP=hl.variant_str(ht_annot.locus, ht_annot.alleles),
        RSID=ht_annot.rsid,
        BP=ht_annot.locus.position,
    )
    ht_annot = ht_annot.key_by().select("CHR", "SNP", "RSID", "BP", *cols)
    ht_annot = checkpoint_tmp(ht_annot)
    ht_annot.drop("RSID").export(get_ld_score_flat_file_path(pop, extension="annot.gz", annot=annot))
    ht_annot.transmute(SNP=ht_annot.RSID).export(
        get_ld_score_flat_file_path(pop, extension="annot.gz", annot=annot, rsid=True)
    )

    # export .l2.ldscore.gz
    ht = ht.join(hm3_snps, "inner")
    ht = ht.annotate(CHR=ht.locus.contig, SNP=hl.variant_str(ht.locus, ht.alleles), RSID=ht.rsid, BP=ht.locus.position)
    ht = ht.key_by().select("CHR", "SNP", "RSID", "BP", *cols)
    ht = ht.rename({key: f"{key}L2" for key in cols})
    ht = checkpoint_tmp(ht)
    # LD score with variant ids
    ht.drop("RSID").export(get_ld_score_flat_file_path(pop, annot=annot))
    # with rsids
    ht.transmute(SNP=ht.RSID).export(get_ld_score_flat_file_path(pop, annot=annot, rsid=True))


def main(args):
    pop = args.pop
    num_pcs = 10
    basic_covars = ["sex", "age", "age2", "age_sex", "age2_sex"]
    covariates = basic_covars + [f"PC{x}" for x in range(1, num_pcs + 1)]

    tmp_mt_path = f"{temp_bucket_7day}/{pop}.mt"
    tmp_bm_path = f"{temp_bucket_7day}/{pop}.bm"
    tmp_r2_adj_path = f"{temp_bucket_7day}/{pop}.r2_adj.bm"

    if args.write_mt:
        mt = get_filtered_mt(chrom="all", pop=pop, entry_fields=["dosage"], min_mac=19, filter_mac_instead_of_ac=True)
        mt_x = get_filtered_mt(chrom="X", pop=pop, entry_fields=["dosage"], min_mac=19, filter_mac_instead_of_ac=True)
        mt = mt.union_rows(mt_x)
        mt = mt.annotate_rows(AF=hl.agg.mean(mt.dosage) / 2)
        # mt = mt.checkpoint(tmp_mt_path, overwrite=args.overwrite)
        n = mt.count()[1]

        # write variant indexes
        ht = mt.rows().select("rsid", "AF").add_index()
        ht = ht.annotate_globals(n_samples=n, pop=pop)
        ht.write(get_ld_variant_index_path(pop), overwrite=args.overwrite)

    if args.write_bm:
        mt = hl.read_matrix_table(tmp_mt_path)
        # convert mt to bm
        BlockMatrix.write_from_entry_expr(
            mt.dosage, tmp_bm_path, mean_impute=True, center=False, normalize=False, overwrite=args.overwrite
        )

    if args.compute_ld_matrix:
        mt = hl.read_matrix_table(tmp_mt_path)
        n = mt.count()[1]
        bm = BlockMatrix.read(tmp_bm_path)
        print(f"BlockMatrix shape: {bm.shape}")

        # mean-center and normalize bm
        bm_norm = normalize_bm(bm)
        bm_norm = checkpoint_tmp(bm_norm)

        # take covariates (with intercept), make hat bms for FWL projection
        cov = mt.cols().select(*covariates).to_pandas().drop(["s"], axis=1)
        cov["Intercept"] = 1.0
        hat1 = cov.values
        hat2 = np.dot(np.linalg.inv(np.dot(cov.transpose(), cov)), cov.transpose())
        bm_hat1 = checkpoint_tmp(BlockMatrix.from_numpy(hat1))
        bm_hat2 = checkpoint_tmp(BlockMatrix.from_numpy(hat2))

        # Cov-adjustement; conducting in three steps due to huge matrix operation
        bm_Z = checkpoint_tmp(bm_norm @ bm_hat1)
        bm_Z = checkpoint_tmp(bm_Z @ bm_hat2)
        bm_Z = checkpoint_tmp(bm_norm - bm_Z)

        # compute ld matrix with a specified radius
        bm_ldadj = (bm_Z @ bm_Z.T) / n
        starts_and_stops = hl.linalg.utils.locus_windows(mt.locus, radius=args.radius, _localize=False)
        bm_ldadj = bm_ldadj._sparsify_row_intervals_expr(starts_and_stops, blocks_only=False)

        # sparcify to a triangle matrix
        bm_ldadj = bm_ldadj.sparsify_triangle()
        bm_ldadj = bm_ldadj.checkpoint(get_ld_matrix_path(pop), overwrite=args.overwrite, force_row_major=True)

    if args.write_ldsc_hm3_snplist:
        # Note: currently, this writes snplists for all the populations at once
        write_ldsc_hm3_snplist(overwrite=args.overwrite)

    if args.compute_ldscore or args.compute_stratified_ldscore:
        ht = hl.read_table(get_ld_variant_index_path(pop))
        n = ht.n_samples.collect()[0]

        if args.overwrite or (not hadoop_exists(tmp_r2_adj_path)):
            bm_ldadj = BlockMatrix.read(get_ld_matrix_path(pop))
            r2 = bm_ldadj ** 2
            r2_adj = ((n - 1.0) / (n - 2.0)) * r2 - (1.0 / (n - 2.0))

            # This is required, as the squaring/multiplication densifies, so this re-sparsifies.
            starts_and_stops = hl.linalg.utils.locus_windows(ht.locus, args.ld_score_radius, _localize=False)
            r2_adj = r2_adj._sparsify_row_intervals_expr(starts_and_stops, blocks_only=False)
            r2_adj = r2_adj.sparsify_triangle()
            r2_adj = r2_adj.checkpoint(tmp_r2_adj_path, overwrite=args.overwrite)
        else:
            r2_adj = BlockMatrix.read(tmp_r2_adj_path)

        if args.compute_ldscore:
            ht_ldscore = copmute_ldscore(
                ht, r2_adj, out_name=get_ld_score_ht_path(pop), overwrite=args.overwrite,
            )
            export_ldscore(ht_ldscore, pop)

        if args.compute_stratified_ldscore:
            annot = args.annot
            ht_annot = hl.read_table(
                f"gs://ukb-diverse-pops/rg-pcgc/annot/ht/{pop}_pananc_31063_pcgc_{annot}_annotations.ht"
            )
            ht_sldscore = copmute_statified_ldscore(
                ht,
                r2_adj,
                out_name=get_ld_score_ht_path(pop, annot=annot),
                overwrite=args.overwrite,
                ht_annot=ht_annot,
            )
            ht_sldscore = hl.read_table(get_ld_score_ht_path(pop, annot=annot))
            export_stratified_ldscore(ht_sldscore, ht_annot, pop, annot)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pop", type=str, required=True, help="Population to compute a LD matrix")
    parser.add_argument("--radius", type=int, default=1e7, help="Radius of window for LD matrix")
    parser.add_argument("--ld-score-radius", type=int, default=1e6, help="Radius of window for LD score")
    parser.add_argument("--write-mt", action="store_true", help="Write MatrixTable from bgen")
    parser.add_argument("--write-bm", action="store_true", help="Write BlockMatrix from MatrixTable")
    parser.add_argument("--compute-ld-matrix", action="store_true", help="Compute LD matrix")
    parser.add_argument("--compute-ldscore", action="store_true", help="Compute LD score")
    parser.add_argument("--compute-stratified-ldscore", action="store_true", help="Compute stratified LD score")
    parser.add_argument("--annot", type=str, help="Annotations for stratified LD score")
    parser.add_argument("--write-ldsc-hm3-snplist", action="store_true", help="Write QCed HM3 snplist for ldsc")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite data")
    args = parser.parse_args()

    atexit.register(lambda: hl.copy_log(timestamp_path(f"gs://ukbb-diverse-temp-30day/ld/{args.pop}/ld", suffix=".log")))

    main(args)
