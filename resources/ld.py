from .generic import *


def get_ld_matrix_path(pop: str, extension: str = "bm"):
    return f"{bucket}/ld/{pop}/UKBB.{pop}.ldadj.{extension}"


def get_ld_score_ht_path(pop: str, annot: str = None):
    annot = f"{annot}." if annot is not None else ""
    return f"{public_bucket}/ld_release/UKBB.{pop}.{annot}ldscore.ht"


def get_ld_score_flat_file_path(pop: str, extension: str = "ldscore.gz", annot: str = None, rsid: bool = False):
    annot = f"{annot}." if annot is not None else ""
    rsid = "rsid." if rsid else ""
    return f"{public_bucket}/ld_release/UKBB.{pop}.{annot}{rsid}l2.{extension}"


def get_ld_variant_index_path(pop: str, extension: str = "ht"):
    return f"{bucket}/ld/{pop}/UKBB.{pop}.ldadj.variant.{extension}"


def get_hm3_snplist_path(pop: str, extension: str = "ht"):
    return f"{bucket}/ld/{pop}/HM3.UKBB.{pop}.qc.snplist.{extension}"
