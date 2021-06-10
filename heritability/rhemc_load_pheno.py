__author__ = 'Rahul Gupta'

import hail as hl
import argparse
import logging
import re

from pcgc_pipeline import saige_phenos, get_pheno_split_names, get_famfiles, \
                          parse_ancestries, _read_pheno_data, get_pheno_filename, \
                          convert_pheno_id_to_potential_saige, htcheckpoint


def compatiblify_phenotype_id(phenotype_id):
    """PCGC throws errors on reading weird phenotype
    headers. This function resolves these characters.
    """
    to_replace = [' ', '|', '>', '<', '/', '\\']
    new_id = phenotype_id
    for i in to_replace:
        new_id = new_id.replace(i, '_')
    return new_id


def rename_phenotype_id(phenotype_id, all_phenotypes):
    """ For some reason, there are a subset of biomarker, continuous, and prescription
    phentoypes that don't match between the manifest and the phenotype MatrixTable. This
    function manually converts these phenotypes to ensure that mapping is successful.
    """
    if phenotype_id in all_phenotypes:
        return phenotype_id
    else:
        if re.search(r'^biomarkers.+-irnt$', phenotype_id):
            new_id = re.sub(r'-irnt$', '', phenotype_id)
        elif re.search(r'^continuous-.+-both_sexes$', phenotype_id):
            middle_term = phenotype_id.split('-')[1]
            new_id = phenotype_id + '-' + middle_term
        elif re.search(r'^prescriptions.+', phenotype_id):
            new_id = re.sub(r'/', '_', phenotype_id)
        else:
            new_id = phenotype_id

        return new_id


def generate_indiv_pheno_file(phenotype_id, checkpoint, log, random=False):
    """ Create individual phenotype file. This does not merge with the fam file.
    """
    pheno_mt_string = _read_pheno_data(checkpoint)
    this_pheno_loc = htcheckpoint + phenotype_id + '.ht'

    if log:
        logging.info(f'Phenotype MatrixTable successfully imported.')

    if random:
        pheno_tab = pheno_mt_string.rows().annotate(value = hl.rand_norm())
        pheno_tab = pheno_tab.rename({'value':compatiblify_phenotype_id(phenotype_id)})
        pheno_tab = pheno_tab.cache()
        if log:
            logging.info(f'Random phenotype table for {phenotype_id} created successfully.')
    else:
        if hl.hadoop_exists(this_pheno_loc):
            pheno_tab = hl.read_table(this_pheno_loc)
        else:
            all_phenos = pheno_mt_string.aggregate_cols(hl.agg.collect_as_set(pheno_mt_string.phenotype_id))
            nmatches = len([x for x in all_phenos if x == phenotype_id])

            if nmatches <= 1:
                new_id = rename_phenotype_id(phenotype_id, all_phenos)
            else:
                raise ValueError('Filtering to phenotype ' + phenotype_id + ' resulted in ' + \
                                str(nmatches) + ' matches. Only one match is allowed.')

            # convert to tab and output
            pheno_mt_filt = pheno_mt_string.filter_cols(pheno_mt_string.phenotype_id == new_id
                                                        ).key_cols_by()
            pheno_tab = pheno_mt_filt.entries().drop('phenotype_id'
                                                    ).rename({'value':compatiblify_phenotype_id(phenotype_id)})
            pheno_tab = pheno_tab.cache()
        
        if log:
            logging.info(f'Phenotype table for {phenotype_id} created successfully.')

    return pheno_tab


def generate_final_pheno_file(phenotype_id, ancestries, checkpoint, override_check, log=False, random=False):
    """ Generate a single phenotype file for each of a set of ancestries. These are in the proper
    format for PCGC. Fam files must exist for this to work, so if they do not exist
    use generate_geno_annot_split.
    """
    if log:
        logging.info(f'Generating phenotype file for {args.phenotype_id}.')

    # load phenotype data
    pheno_tab = generate_indiv_pheno_file(phenotype_id, checkpoint, log=log, random=random)
    famfiles = get_famfiles(ancestries, checkpoint, override_check=override_check)

    # join with fam files and output
    pheno_dirs = get_pheno_split_names(ancestries, 
                                       phenotype_id=compatiblify_phenotype_id(phenotype_id), 
                                       dictout=True)
    for anc in ancestries:
        famfile_this = famfiles[anc]
        famfile_this = famfile_this.add_index('idx').key_by('IID')
        famfile_pheno = famfile_this.join(pheno_tab, how='left')
        famfile_pheno = famfile_pheno.key_by('idx').order_by('idx').drop('idx')
        if log:
            logging.info(f'{anc} .fam file merged with phenotype for {args.phenotype_id}.')
        famfile_pheno.export(output=pheno_dirs[anc], header=True, delimiter='\t')
        if log:
            logging.info(f'{pheno_dirs[anc]} written.')


def format_saige_phenotype(phenotype_id, anc, log):
    """ Imports a Saige phenotype flat file as a HailTable and formats to a similar schema
    as generate_indiv_pheno_file.
    """
    saige_file = saige_phenos + anc + '/' + convert_pheno_id_to_potential_saige(phenotype_id)
    tab = hl.import_table(saige_file, impute=True).repartition(n=10).key_by('userId')
    if log:
        logging.info(f'Phenotype table imported from {saige_file}')
    tab = tab.select('value')
    tab = tab.annotate(value=hl.float64(tab.value))
    tab = tab.rename({'value': compatiblify_phenotype_id(phenotype_id)})
    return tab


def convert_saige_to_pcgc(phenotype_id, ancestries, checkpoint, override_check, log=False):
    """ Converts a pre-existing filtered phenotype flat file (constructed for Saige)
    to the proper format for PCGC. The results are uploaded in the usual format dictated
    by `get_pheno_split_names`, allowing for use by ancestry-level workers.
    """
    if log:
        logging.info(f'Generating phenotype file for {args.phenotype_id} from Saige files.')
    
    famfiles = get_famfiles(ancestries, checkpoint, override_check=override_check)   
    pheno_dirs = get_pheno_split_names(ancestries, 
                                       phenotype_id=compatiblify_phenotype_id(phenotype_id), 
                                       dictout=True)
    for anc in ancestries:
        # load phenotype data
        pheno_tab = format_saige_phenotype(phenotype_id, anc, log=log)

        # join with fam files
        famfile_this = famfiles[anc]
        famfile_this = famfile_this.add_index('idx').key_by('IID')
        famfile_pheno = famfile_this.join(pheno_tab, how='left')
        famfile_pheno = famfile_pheno.key_by('idx').order_by('idx').drop('idx')
        if log:
            logging.info(f'{anc} .fam file merged with phenotype for {args.phenotype_id}.')
        
        famfile_pheno.export(output=pheno_dirs[anc], header=True, delimiter='\t')
        if log:
            logging.info(f'{pheno_dirs[anc]} written.')


parser = argparse.ArgumentParser()
parser.add_argument('--phenotype-id', type=str,
                    help='Phenotype identifier. See pcgc_pipeline.construct_phenotype_id ' + \
                         'to see how these are constructed.')
parser.add_argument('--logging', action='store_true',
                    help='If enabled, outputs a text log')
parser.add_argument('--ancestries', type=str,
                    help='Comma-delimited set of ancestries to include.')
parser.add_argument('--checkpoint', action='store_true',
                    help='If enabled, will checkpoint all individual level data in '+ \
                         'Hail format prior to outputting flat/binary files in check_indiv_files.')
parser.add_argument('--override-check', action='store_true',
                    help='If enabled, will assume that fam files in Hail Table format have been created already.')
parser.add_argument('--pull-from-saige', action='store_true',
                    help='If used, will assume that the file exists as a .tsv from the Saige run.  ' + \
                         'This tool will NOT check for this, so do ensure that the file exists.  ' + \
                         'A try catch sequence will be used; if the import fails, this will fall back ' + \
                         'to stock import.')
parser.add_argument('--random', action='store_true',
                    help='If used, will generate a phenotype at random form a Normal(0,1). Draws ' + \
                         'will be independent across ancestries.')


if __name__ == '__main__':
    args = parser.parse_args()
    ancestries = parse_ancestries(args)
    if args.logging:
        log_name = get_pheno_filename(compatiblify_phenotype_id(args.phenotype_id), enable_suffix=False) + '.log'
        logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s", level='INFO', filename=log_name)

    if args.pull_from_saige:
        try:
            convert_saige_to_pcgc(args.phenotype_id, ancestries,
                                  checkpoint=args.checkpoint, 
                                  override_check=args.override_check, 
                                  log=args.logging)
        except Exception as e:
            if args.logging:
                logging.info(f'Though Saige phenotype import was used for {args.phenotype_id}, an excpetion was raised: {e}.' + \
                            f'Trying to generate this phenotype from the original MatrixTable.')
            generate_final_pheno_file(args.phenotype_id, ancestries, 
                                      checkpoint=args.checkpoint, 
                                      override_check=args.override_check, 
                                      log=args.logging, random=args.random)
    else:
        generate_final_pheno_file(args.phenotype_id, ancestries, 
                                  checkpoint=args.checkpoint, 
                                  override_check=args.override_check,
                                  log=args.logging, random=args.random)
