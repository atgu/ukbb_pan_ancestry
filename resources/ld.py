from .generic import *


def get_ld_matrix_path(pop: str, extension: str = 'bm'):
    return f'{bucket}/ld/{pop}/UKBB.{pop}.ldadj.{extension}'


def get_ld_score_ht_path(pop: str):
    return f'{public_bucket}/ld_release/UKBB.{pop}.ldscore.ht'


def get_ld_score_flat_file_path(pop: str, extension: str = 'ldscore.gz', rsid: bool = False):
    if extension == 'ldscore.gz' and rsid:
        return f'{public_bucket}/ld_release/UKBB.{pop}.rsid.l2.{extension}'
    return f'{public_bucket}/ld_release/UKBB.{pop}.l2.{extension}'


def get_ld_variant_index_path(pop: str, extension: str = 'ht'):
    return f'{bucket}/ld/{pop}/UKBB.{pop}.ldadj.variant.{extension}'


def get_hm3_snplist_path(pop: str, extension: str = 'ht'):
    return f'{bucket}/ld/{pop}/HM3.UKBB.{pop}.qc.snplist.{extension}'