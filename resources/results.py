from .generic import *


def get_variant_results_path(pop: str, extension: str = 'mt'):
    return f'{bucket}/combined_results/results_{pop}.{extension}'

