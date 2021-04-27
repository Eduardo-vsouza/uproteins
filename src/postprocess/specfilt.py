import os
import pandas as pd

from pyteogenomics.postsearch import SequenceFinder, LinkData
from pyteogenomics.utilities import PercolatorConverter
from pyteogenomics import StringTieGFF, GenomeCoordinates, RefSeqGFF, GenomeCoordinatesRNA, PercolatorUTP, StillCounting, Enrichment


class PostPercolator(object):
    def __init__(self, args, folder, filetype):
        self.args = args
        self.folder = folder
        self.percDir = f'{self.folder}/post_perc'
        self.filetype = filetype
        self.__check_dir()

    def __check_dir(self):
        if not os.path.exists(f'{self.folder}/post_perc'):
            os.system(f'mkdir {self.folder}/post_perc')

    def convert_output(self):
        """ Fix some inconsistencies. Pre-processing step. """
        pout = PercolatorConverter(pout=f'{self.percDir}/{self.filetype}_results_psm.txt', conversion_file='.',
                                   handle='psm', gff=self.args.gff)
        pout.convert_entries(pattern='WP', output=f"{self.percDir}/{self.filetype}_converted_psm.txt")

    def get_coordinates_rna(self):
        """ Finds the genome coordinates of each transcript. Must be used on transcriptome database exclusively. """
        str_gff = StringTieGFF(gff='assembled.gtf')
        orf_dict = str_gff.get_dict()
        ref_gff = RefSeqGFF(gff=self.args.gff)
        ref_dict = ref_gff.get_dict()
        coordinates = GenomeCoordinatesRNA(f'{self.percDir}/{self.filetype}_converted_psm.txt', orf_dict,
                                           ref_dict)
        coordinates.get_proteins().save_table(output=f'{self.percDir}/{self.filetype}_psm_coords')
        return self

    def get_coordinates_genome(self):
        """ Finds the genome coordinates of each ORF from the genome database. """
        ref_gff = RefSeqGFF(gff=self.args.gff)
        ref_dict = ref_gff.get_dict()
        coordinates = GenomeCoordinates(f'{self.percDir}/{self.filetype}_converted_psm.txt', ref_dict)
        coordinates.get_coords().save_table(output=f'{self.percDir}/{self.filetype}_psm_coords')
        return self

    def filter_novel(self):
        anno = AnnoFilter(f'{self.percDir}/{self.filetype}_psm_coords.txt')
        anno.remove_annotated(f'{self.percDir}/{self.filetype}_no_anno.txt')

    def unique_peptides(self):
        """ Remove non-unique peptides. Check uProteInS methods for unique peptide classification. """
        unique = PercolatorUTP(coord_df=f'{self.percDir}/{self.filetype}_no_anno.txt', pep=self.args.pep,
                               qvalue=self.args.qvalue)
        unique.get_utps().save_table(output=f'{self.percDir}/{self.filetype}_utps', keep='all')
        return self

    def msgf_info(self):
        """ Adds MSGF info to percolator's output. """
        os.system(f'cat {self.folder}/tsv_msgf/*.tsv > {self.percDir}/{self.filetype}_search.tsv')
        link = LinkData(f'{self.percDir}/{self.filetype}_search.tsv', f'{self.percDir}/{self.filetype}_utps.txt')
        link.filter_msgf(f'{self.percDir}/{self.filetype}_linked')

    def protein_seqs(self):
        seq = SequenceFinder(f'{self.percDir}/{self.filetype}_linked.tsv', f'{self.filetype}_database.fasta')
        seq.df_proteins().save(f'{self.percDir}/{self.filetype}_proteined')



class AnnoFilter(object):
    def __init__(self, df):
        self.df = pd.read_csv(df, sep='\t')

    def remove_annotated(self, output):
        df = self.df[self.df["proteinIds"].str.contains('ANNO') == False]
        df.to_csv(output, sep='\t', index=False)

