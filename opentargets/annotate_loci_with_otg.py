import argparse
import requests
import numpy as np
import pandas as pd
import hail as hl
from hail.expr.functions import _sort_by
from hail.utils import hadoop_open, new_temp_file
import pyranges as pr
from ukbb_pan_ancestry import (
    get_meta_analysis_results_path,
    get_pheno_manifest_path,
    get_distance_clumping_results_path,
    get_ukb_pheno_efo_mapping_path,
    get_munged_otg_v2d_path,
    get_known_ukbb_loci_path,
    load_final_sumstats_mt,
    otg_release,
)

OPENTARGET_PREFIX = "gs://ukbb-diverse-open-targets-genetics-releases/releases"

def get_efo_otg_mappings(release: str = otg_release):
    # EFO mappings from OpenTargets Genetics
    ht = hl.read_table(f"{OPENTARGET_PREFIX}/{release}/lut/study-index.ht")

    ht_neale = ht.filter(ht.study_id.startswith("NEALE2_"))
    ht_neale = ht_neale.annotate(study_id=ht_neale.study_id.replace("^NEALE2_", "").replace("_raw$", ""))
    ht_neale = ht_neale.annotate(x=ht_neale.study_id.split("_", 2))
    ht_neale = ht_neale.select(
        phenocode=ht_neale.x[0],
        coding=hl.if_else(hl.len(ht_neale.x) == 2, ht_neale.x[1], ""),
        trait_efos=ht_neale.trait_efos,
    )
    # trait type unknown
    ht_neale = ht_neale.key_by("phenocode", "coding").distinct()

    ht_phecode = ht.filter(ht.study_id.startswith("SAIGE_"))
    ht_phecode = ht_phecode.select(
        trait_type="phecode",
        phenocode=ht_phecode.study_id.replace("^SAIGE_", "").replace("_", "."),
        coding="",
        trait_efos=ht_phecode.trait_efos,
    )
    ht_phecode = ht_phecode.key_by("trait_type", "phenocode", "coding").distinct()

    return ht_phecode, ht_neale


def get_efo_ukbb_mappings():
    # EFO/UKBB mappings from https://github.com/EBISPOT/EFO-UKB-mappings
    r = requests.get(
        "https://raw.githubusercontent.com/EBISPOT/EFO-UKB-mappings/4ed870a3fc127622b14ead663fcba1a1ff9a8231/UK_Biobank_master_file.tsv"
    )
    lines = r.text.split("\n")
    # fix offending lines
    lines[1134] = lines[1134].replace(" ", "\t")
    for i in [1135, 1136, 1137, 1312]:
        lines[i] = lines[i].rstrip() + "\t\r"

    tmpfile = new_temp_file()
    with hadoop_open(tmpfile, "w") as f:
        f.write("\n".join(lines))

    ht_efo = hl.import_table(tmpfile, impute=True)
    ht_efo = ht_efo.select(
        trait_type="",
        phenocode=ht_efo["ICD10_CODE/SELF_REPORTED_TRAIT_FIELD_CODE"],
        coding="",
        trait_efos=ht_efo["MAPPED_TERM_URI"].split("[,| ]+"),
    )
    ht_efo = ht_efo.annotate(
        phenocode=hl.case()
        .when(ht_efo.phenocode.contains("_"), ht_efo.phenocode.split("_")[0])
        .default(ht_efo.phenocode),
        coding=hl.case().when(ht_efo.phenocode.contains("_"), ht_efo.phenocode.split("_")[1]).default(""),
    )

    ht_icd10 = ht_efo.filter(ht_efo.phenocode.matches("^[A-Z]"))
    ht_icd10 = ht_icd10.annotate(trait_type="icd10")
    ht_icd10 = ht_icd10.key_by("trait_type", "phenocode", "coding").distinct()
    # trait type unknown
    ht_efo = ht_efo.filter(~ht_efo.phenocode.matches("^[A-Z]"))
    ht_efo = ht_efo.key_by("phenocode", "coding").distinct()

    return ht_icd10, ht_efo


def get_efo_manual_mappings():
    # manual mappings for biomarkers
    ht_biomarker = hl.import_table("gs://ukb-diverse-pops/Phenotypes/biomarkers_efo.txt")
    ht_biomarker = ht_biomarker.select(
        trait_type=ht_biomarker.trait_type,
        phenocode=ht_biomarker.phenocode,
        coding="",
        trait_efos=[ht_biomarker.trait_efos],
    )
    ht_biomarker = ht_biomarker.key_by("trait_type", "phenocode", "coding").distinct()
    ht_continuous = ht_biomarker.filter(ht_biomarker.trait_type == "continuous")
    ht_biomarker = ht_biomarker.filter(ht_biomarker.trait_type == "biomarkers")

    return ht_continuous, ht_biomarker


def get_obsolete_efo_dict():
    # see https://github.com/EBISPOT/efo/issues/1381
    with hadoop_open("gs://ukb-diverse-pops/Phenotypes/obsoletes_replacement.simplified.tsv", "r") as f:
        df = pd.read_csv(f, delimiter="\t")
    df["cls"] = df["cls"].str.split("/").str[-1]
    df["replCls"] = df["replCls"].str.split("/").str[-1]
    obsolete_efo_dict = hl.literal(
        df[["cls", "replCls"]].drop_duplicates().set_index("cls").T.to_dict("records")[0]
    )
    return obsolete_efo_dict


def get_efo_dicts():
    # GWAS catalog EFO Parent category
    df = pd.read_csv("https://www.ebi.ac.uk/gwas/api/search/downloads/trait_mappings", delimiter="\t")
    df["efo"] = df["EFO URI"].str.split("/").str[-1]
    efo_category = hl.literal(df[["efo", "Parent term"]].drop_duplicates().set_index("efo").T.to_dict("records")[0])
    efo_term = hl.literal(df[["efo", "EFO term"]].drop_duplicates().set_index("efo").T.to_dict("records")[0])

    return efo_category, efo_term


def munge_otg_variant_index(release: str = otg_release):
    ht = hl.read_table(f"{OPENTARGET_PREFIX}/{release}/variant-index.ht")
    ht = ht.annotate(
        locus_b38=hl.locus("chr" + ht.chr_id, ht.position, reference_genome="GRCh38"),
        locus=hl.locus(ht.chr_id_b37, ht.position_b37, reference_genome="GRCh37"),
    )
    ht = ht.filter(hl.is_defined(ht.locus))
    ht = ht.key_by("locus_b38")
    return ht


def munge_otg_v2d(efo_category, efo_term, release: str = otg_release):
    ukbb_sources = hl.set(["NEALE", "SAIGE"])

    ht = hl.read_table(f"{OPENTARGET_PREFIX}/{release}/v2d.ht")
    ht = ht.filter(~ukbb_sources.contains(ht.source) & (hl.len(ht.trait_efos) > 0))
    ht = ht.key_by(locus_b38=hl.locus("chr" + ht.lead_chrom, hl.int32(ht.lead_pos), reference_genome="GRCh38"))

    ht_variant = munge_otg_variant_index()
    ht = ht.join(ht_variant, "inner")

    ht = ht.select("locus", "study_id", "trait_efos", "trait_reported")
    ht = ht.explode("trait_efos", name="trait_efo")

    obsolete_efo_dict = get_obsolete_efo_dict()
    ht = ht.annotate(trait_efo=obsolete_efo_dict.get(ht.trait_efo, ht.trait_efo))

    ht = ht.group_by("locus", "trait_efo").aggregate(study_id=hl.agg.collect_as_set(ht.study_id))
    ht = ht.annotate(
        trait_efo_term=hl.or_missing(efo_term.contains(ht.trait_efo), efo_term[ht.trait_efo]),
        trait_efo_category=hl.or_missing(efo_category.contains(ht.trait_efo), efo_category[ht.trait_efo]),
    )
    return ht


def agg_distance_clumping(
    locus_expr, alleles_expr, pvalue_expr, radius=500000, nlog_p=False, merge_overlapping_loci=True
):
    def _distance_clumping(f, lead_st_arr, remaining_st_arr):
        new_lead_st = remaining_st_arr[0]
        new_region = new_lead_st.region
        remaining_st_arr = remaining_st_arr.filter(lambda x: ~new_region.contains(x.lead_locus))

        # check whether a new region overlaps with previous ones
        overlapping_st = (
            lead_st_arr.find(lambda x: new_region.overlaps(x.region))
            if merge_overlapping_loci
            else hl.missing(new_lead_st.dtype)
        )
        lead_st_arr = hl.if_else(
            hl.is_defined(overlapping_st),
            lead_st_arr.map(
                lambda x: hl.if_else(
                    x == overlapping_st,
                    x.annotate(
                        region=hl.rbind(
                            overlapping_st.region,
                            lambda old_region: hl.interval(
                                hl.if_else(
                                    old_region.start < new_region.start,
                                    old_region.start,
                                    new_region.start,
                                ),
                                hl.if_else(old_region.end > new_region.end, old_region.end, new_region.end),
                            ),
                        )
                    ),
                    x,
                )
            ),
            lead_st_arr.append(new_lead_st),
        )
        converged = hl.len(remaining_st_arr) == 0

        return hl.if_else(converged, lead_st_arr, f(lead_st_arr, remaining_st_arr))

    # Modified from hl.expr.functions.sorted -- key is now another collection, not callable
    def _sorted(collection, key, reverse=False):
        def comp(left, right):
            return (
                hl.case()
                .when(hl.is_missing(left), False)
                .when(hl.is_missing(right), True)
                .when(reverse, hl._compare(right, left) < 0)
                .default(hl._compare(left, right) < 0)
            )

        return _sort_by(hl.zip(key, collection), lambda l, r: comp(l[0], r[0])).map(lambda elt: elt[1])

    # lead locus / region pair is necessary to check overlapping regions during the loop
    x = hl.agg.collect(
        hl.struct(
            st=hl.struct(
                lead_locus=locus_expr,
                lead_alleles=alleles_expr,
                region=locus_expr.window(radius, radius),
                lead_pvalue=pvalue_expr,
            )
        )
    )
    sorted_st = _sorted(x.st, x.st.lead_pvalue, reverse=nlog_p)
    t_struct = sorted_st[0].dtype

    run_distance_clumping = hl.experimental.define_function(
        lambda sorted_st: hl.experimental.loop(
            _distance_clumping, hl.tarray(t_struct), hl.empty_array(t_struct), sorted_st
        ),
        hl.tarray(t_struct),
    )
    return run_distance_clumping(sorted_st)


def main(args):
    efo_category, efo_term = get_efo_dicts()
    obsolete_efo_dict = get_obsolete_efo_dict()

    if args.annotate_efo_mappings:
        ht_phecode, ht_neale = get_efo_otg_mappings()
        ht_icd10, ht_efo = get_efo_ukbb_mappings()
        ht_continuous, ht_biomarker = get_efo_manual_mappings()
        ht_union = hl.Table.union(ht_phecode, ht_icd10, ht_biomarker, ht_continuous).cache()

        neale_trait_types = hl.set(["categorical", "continuous"])

        ht = hl.import_table(
            # get_pheno_manifest_path(),
            "gs://ukb-diverse-pops/combined_results/221215_phenotype_manifest.tsv.bgz",
            key=["trait_type", "phenocode", "pheno_sex", "coding", "modifier"], impute=True
        )
        ht = ht.annotate(
            trait_efos=hl.case()
            .when(
                hl.is_defined(ht_union[ht.trait_type, ht.phenocode, ht.coding]),
                ht_union[ht.trait_type, ht.phenocode, ht.coding].trait_efos,
            )
            .when(
                neale_trait_types.contains(ht.trait_type) & hl.is_defined(ht_neale[ht.phenocode, ht.coding]),
                ht_neale[ht.phenocode, ht.coding].trait_efos,
            )
            .when(
                neale_trait_types.contains(ht.trait_type) & hl.is_defined(ht_efo[ht.phenocode, ht.coding]),
                ht_efo[ht.phenocode, ht.coding].trait_efos,
            )
            .when((ht.phenocode == "Smoking") & ht.modifier.startswith("CPD_combined"), ["EFO_0006525"])
            .when((ht.phenocode == "Smoking") & (ht.modifier == "Ever_Never"), ["EFO_0006527"])
            .default(hl.missing("array<str>"))
        )
        # update obsolete EFO terms
        ht = ht.annotate(trait_efos=hl.map(lambda x: obsolete_efo_dict.get(x, x), ht.trait_efos))
        # annotate EFO term and categories
        ht = ht.annotate(
            trait_efo_terms=hl.or_missing(
                hl.is_defined(ht.trait_efos),
                ht.trait_efos.map(lambda x: hl.or_missing(efo_term.contains(x), efo_term[x])),
            ),
            trait_efo_categories=hl.or_missing(
                hl.is_defined(ht.trait_efos),
                hl.set(ht.trait_efos.map(lambda x: hl.or_missing(efo_category.contains(x), efo_category[x]))).filter(
                    hl.is_defined
                ),
            ),
        )
        ht_efo = ht.checkpoint(get_ukb_pheno_efo_mapping_path(), overwrite=args.overwrite)

        ht_efo.group_by(ht_efo.trait_type).aggregate(
            n_annotated=hl.agg.sum(hl.is_defined(ht_efo.trait_efos)),
            n_efo_defined=hl.agg.sum(hl.len(ht_efo.trait_efos) > 0),
            n_efo_missing=hl.agg.sum(hl.is_missing(ht_efo.trait_efos)),
            n_category_defined=hl.agg.sum(hl.len(ht_efo.trait_efo_categories) > 0),
            n_category_multi=hl.agg.sum(hl.len(ht_efo.trait_efo_categories) > 1),
        ).show()
    else:
        ht_efo = hl.read_table(get_ukb_pheno_efo_mapping_path())

    if args.munge_otg_v2d:
        ht_v2d = munge_otg_v2d(efo_category, efo_term)
        ht_v2d = ht_v2d.checkpoint(get_munged_otg_v2d_path(), overwrite=args.overwrite)
    else:
        ht_v2d = hl.read_table(get_munged_otg_v2d_path())

    if args.distance_clump:
        if args.sumstats_pop == "meta_hq":
            mt = hl.read_matrix_table(get_meta_analysis_results_path(filter_pheno_h2_qc=True))
            ht = mt.entries()
            ht = ht.annotate(Pvalue=ht.meta_analysis.Pvalue[0])
        elif args.sumstats_pop == "meta_raw":
            mt = hl.read_matrix_table(get_meta_analysis_results_path(filter_pheno_h2_qc=False))
            mt_hq = hl.read_matrix_table(get_meta_analysis_results_path(filter_pheno_h2_qc=True))
            mt = mt.semi_join_cols(mt_hq.cols())
            ht = mt.entries()
            ht = ht.annotate(Pvalue=ht.meta_analysis.Pvalue[0])
        elif args.sumstats_pop == "EUR":
            mt = load_final_sumstats_mt(
                filter_phenos=True,
                filter_variants=False,
                filter_sumstats=True,
                separate_columns_by_pop=True,
                annotate_with_nearest_gene=False,
                filter_pheno_h2_qc=True,
            )
            mt = mt.filter_cols(mt.pheno_data.pop == "EUR")
            ht = mt.entries()
            ht = ht.annotate(Pvalue=ht.summary_stats.Pvalue)
        elif args.sumstats_pop == "EUR_round2":
            # TODO: match phenotype keys
            mt = hl.experimental.load_dataset(
                name="UK_Biobank_Rapid_GWAS_both_sexes", version="v2", reference_genome="GRCh37"
            )
            mt = mt.filter_cols(mt.variable_type != "continuous_raw")
            ht = mt.entries()
            ht = ht.annotate(Pvalue=ht.pval)
        else:
            raise ValueError()

        ht = ht.filter(ht.Pvalue < np.log(args.p_threshold))
        ht = ht.group_by(*list(mt.col_key)).aggregate(
            distance_clumps=agg_distance_clumping(
                ht.locus,
                ht.alleles,
                ht.Pvalue,
                radius=args.radius,
                nlog_p=False,
                merge_overlapping_loci=(not args.not_merge_overlapping_loci),
            )
        )
        ht = ht.annotate(n_lead_loci=hl.len(ht.distance_clumps))
        ht = ht.checkpoint(
            get_distance_clumping_results_path(
                pop=args.sumstats_pop, radius=args.radius, merged=(not args.not_merge_overlapping_loci)
            ),
            overwrite=args.overwrite,
        )
        ht.drop("distance_clumps").export(
            get_distance_clumping_results_path(
                pop=args.sumstats_pop,
                radius=args.radius,
                merged=(not args.not_merge_overlapping_loci),
                extension="n_lead_loci.tsv.bgz",
            )
        )
    else:
        ht = hl.read_table(
            get_distance_clumping_results_path(
                pop=args.sumstats_pop, radius=args.radius, merged=(not args.not_merge_overlapping_loci)
            )
        )

    if args.annotate_known_loci:

        def _to_pandas(ht):
            df = ht.to_pandas()
            # a temporary fix: cf. https://github.com/biocore-ntnu/pyranges/pull/264
            df = df.astype({k: "object" if v == "string[python]" else str(v).lower() for k, v in df.dtypes.items()})
            return df

        # munge v2d for pyranges
        ht_v2d = ht_v2d.annotate(
            study_id=hl.delimit(ht_v2d.study_id),
            Chromosome=ht_v2d.locus.contig,
            Start=ht_v2d.locus.position,
            End=ht_v2d.locus.position,
        )
        df_v2d = _to_pandas(ht_v2d)

        # munge efo
        ht_efo = ht_efo.explode("trait_efos", name="trait_efo")
        ht_efo = ht_efo.annotate(
            trait_efo_term=hl.or_missing(efo_term.contains(ht_efo.trait_efo), efo_term[ht_efo.trait_efo]),
            trait_efo_category=hl.or_missing(efo_category.contains(ht_efo.trait_efo), efo_category[ht_efo.trait_efo]),
        )
        ht_efo = ht_efo.select(
            "description", "description_more", "coding_description", "trait_efo", "trait_efo_term", "trait_efo_category"
        )

        ht = ht.drop("n_lead_loci").explode("distance_clumps")
        ht = ht.transmute(**ht.distance_clumps)
        # unique index for trait-locus pairs
        ht = ht.add_index()
        # this join will duplicate trait-locus pairs if a trait matches multilple trait_efo
        ht = ht.join(ht_efo, "left")
        ht2 = ht.annotate(Chromosome=ht.region.start.contig, Start=ht.region.start.position, End=ht.region.end.position)
        df = _to_pandas(ht2)

        gr = pr.PyRanges(df)
        gr_v2d = pr.PyRanges(df_v2d)
        df_annotated = gr.join(gr_v2d, how="left", suffix="_otg").as_df()
        df_annotated = df_annotated.drop(columns=["Chromosome", "Start", "End", "Start_otg", "End_otg"])

        # write out and read back b/c hl.Table.from_pandas is extreamly slow
        tmpfile = new_temp_file()
        with hadoop_open(tmpfile, "wb") as f:
            df_annotated.to_csv(f, sep="\t", na_rep="NA", index=False)

        ht = hl.import_table(tmpfile, impute=True, min_partitions=1000)
        ht = ht.rename({"locus": "locus_otg", "study_id": "study_id_otg"})
        ht = ht.annotate(
            **{
                key: hl.if_else(ht[key] == "-1", hl.missing(hl.tstr), ht[key])
                for key in filter(lambda x: x.endswith("_otg"), list(ht.row))
            }
        )
        ht = ht.annotate(
            lead_locus=hl.parse_locus(ht.lead_locus),
            region=hl.parse_locus_interval(ht.region),
            locus=hl.parse_locus(ht.locus_otg),
            study_id_otg=ht.study_id_otg.split(","),
        )

        for field in ["trait_efo", "trait_efo_category"]:
            x = df_annotated[~df_annotated[field].isna()]
            known_idx = set(x[x[field] == x[f"{field}_otg"]].idx)
            novel_idx = set(x.idx) - set(known_idx)
            na_idx = set(df_annotated.idx) - set(x.idx)
            na_match = na_idx.intersection(
                set(df_annotated[df_annotated[field].isna() & ~df_annotated.locus.isna()].idx)
            )

            k = "efo" if field == "trait_efo" else "efo_category"
            ht = ht.annotate(
                **{
                    f"is_known_{k}": hl.case()
                    .when(hl.set(known_idx).contains(ht.idx), True)
                    .when(hl.set(novel_idx).contains(ht.idx), False)
                    .default(hl.missing(hl.tbool))
                }
            )

            print(field)
            # category known
            print(len(known_idx))
            # category novel
            print(len(novel_idx))
            # category NA
            print(len(na_idx))
            # category NA match
            print(len(na_match))

        ht = ht.key_by("trait_type", "phenocode", "pheno_sex", "coding", "modifier", "idx")
        ht = ht.checkpoint(get_known_ukbb_loci_path(pop=args.sumstats_pop), overwrite=args.overwrite)
        ht.export(get_known_ukbb_loci_path(pop=args.sumstats_pop, extension="tsv.bgz"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--annotate-efo-mappings",
        help="XXXX",
        action="store_true",
    )
    parser.add_argument(
        "--munge-otg-v2d",
        help="XXXX",
        action="store_true",
    )
    parser.add_argument(
        "--keep-low-heritability-pairs",
        help="Keep trait-ancestry pairs that failed at heritability QC",
        action="store_true",
    )
    parser.add_argument(
        "--distance-clump",
        help="XXXX",
        action="store_true",
    )
    parser.add_argument(
        "--sumstats-pop", default="meta_hq", type=str, choices=["meta_hq", "meta_raw", "EUR", "EUR_round2"], help="XXXX"
    )
    parser.add_argument(
        "--not-merge-overlapping-loci",
        help="XXXX",
        action="store_true",
    )
    parser.add_argument(
        "--annotate-known-loci",
        help="XXXX",
        action="store_true",
    )
    parser.add_argument(
        "--otg-release",
        default=otg_release,
        help="Release version of OpenTargets Genetics",
    )
    parser.add_argument(
        "--p-threshold",
        type=float,
        default=5e-8,
        help="P-value threshold for genome-wide significane",
    )
    parser.add_argument(
        "--radius",
        type=int,
        default=500000,
        help="Window radius for distance-based clumping",
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    args = parser.parse_args()

    main(args)
