import os
import pandas as pd
from Bio import SeqIO
import numpy as np

from ..sequtils.postsearch import SequenceFinder, LinkData, TSVChunks
from ..sequtils.utilities import PercolatorConverter
from ..sequtils import StringTieGFF, GenomeCoordinates, RefSeqGFF, GenomeCoordinatesRNA, PercolatorUTP, StillCounting, Enrichment


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
        print('Converting output\n')
        pout = PercolatorConverter(pout=f'{self.folder}/Percolator/{self.filetype}_results_psm.txt', conversion_file='.',
                                   handle='psm', gff=self.args.gff)
        pout.convert_entries(pattern='WP', output=f"{self.percDir}/{self.filetype}_converted_psm")

    def get_coordinates_rna(self):
        """ Finds the genome coordinates of each transcript. Must be used on transcriptome database exclusively. """
        print('Getting coordinates\n')
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
        print('Getting coordinates\n')
        ref_gff = RefSeqGFF(gff=self.args.gff)
        ref_dict = ref_gff.get_dict()
        coordinates = GenomeCoordinates(f'{self.percDir}/{self.filetype}_converted_psm.txt', ref_dict)
        coordinates.get_coords().save_table(output=f'{self.percDir}/{self.filetype}_psm_coords')
        return self

    def filter_novel(self):
        print('Filtering novel peptides\n')
        anno = AnnoFilter(f'{self.percDir}/{self.filetype}_psm_coords.txt')
        anno.remove_annotated(f'{self.percDir}/{self.filetype}_no_anno.txt')

    def unique_peptides(self):
        """ Remove non-unique peptides. Check uProteInS methods for unique peptide classification. """
        print("Removing non-unique peptides\n")
        os.system(f'cp {self.percDir}/{self.filetype}_no_anno.txt {self.percDir}/{self.filetype}_utps.txt')
        # unique = PercolatorUTP(coord_df=f'{self.percDir}/{self.filetype}_no_anno.txt', pep=self.args.pep,
        #                        qvalue=self.args.qvalue)
        # unique.get_utps().save_table(output=f'{self.percDir}/{self.filetype}_utps', keep='all')
        return self

    def msgf_info(self):
        """ Adds MSGF info to percolator's output. """
        print('Adding MSGF info\n')
        os.system(f'cat {self.folder}/tsv_msgf/*.tsv > {self.percDir}/{self.filetype}_search.tsv')
        if self.args.Transcriptome == "YES" and self.filetype == 'transcriptome':
            chunks = TSVChunks(folder='Transcriptome', filetype='transcriptome')
            chunks.filter_search()
        if self.filetype == 'genome':
            chunks = TSVChunks(folder=self.folder, filetype=self.filetype)
            chunks.filter_search()
        link = LinkData(f'{self.percDir}/{self.filetype}_chunk_search.tsv', f'{self.percDir}/{self.filetype}_utps.txt')
        link.filter_msgf(f'{self.percDir}/{self.filetype}_linked')

    def protein_seqs(self):
        """
        Adds protein seqs to the output.
        """
        print('Adding protein sequences\n')
        seq = SequenceFinder(f'{self.percDir}/{self.filetype}_linked.tsv', f'{self.filetype}_database.fasta')
        seq.df_proteins().save(f'{self.percDir}/{self.filetype}_proteined')

    def protein_threshold(self):
        """
        Applies a protein FDR to the results.
        """
        print('Protein threshold\n')
        protein = ProteinFDR(folder=self.folder, filetype=self.filetype)
        protein.protein_cutoff()
        protein.apply_to_psm()
        protein.filter_from_utp()
        protein.add_proteins(f'{self.percDir}/{self.filetype}_results_02.txt')

    def add_coordinates(self, qvalue=0.01):
        coords = Coordinator(proteined=f'{self.percDir}/{self.filetype}_proteined.tsv', utps=f'{self.percDir}/{self.filetype}_utps.txt', qvalue=qvalue)
        coords.add_information(f'{self.percDir}/{self.filetype}_results_02.txt')


class AnnoFilter(object):
    def __init__(self, df):
        self.df = pd.read_csv(df, sep='\t')

    def remove_annotated(self, output):
        df = self.df[self.df["proteinIds"].str.contains('ANNO') == False]
        df.to_csv(output, sep='\t', index=False)


class Coordinator(object):
    def __init__(self, utps, proteined, qvalue=0.01):
        self.UTPs = pd.read_csv(utps, sep='\t')
        print(self.UTPs)
        self.UTPs = self.UTPs[self.UTPs["q-value"] != "q-value"]
        print(self.UTPs)
        # print(self.UTPs)
        self.UTPs["q-value"] = pd.to_numeric(self.UTPs["q-value"], downcast='float')
        print(self.UTPs)
        # self.UTPs["q-value"] = self.UTPs["q-value"].astype(float)
        # qvs = self.UTPs["q-value"].tolist()
        # for i in qvs:
        #     if type(i) != float:
        #         print(i)
        self.UTPs = self.UTPs[self.UTPs["q-value"] <= qvalue]
        print(self.UTPs)
        self.proteined = pd.read_csv(proteined, sep='\t')
        self.coordinates = self._get_coordinates()

    def _get_coordinates(self):
        coordinates = self.UTPs["Genome Coordinates"].tolist()
        orfs = self.UTPs["proteinIds"].tolist()
        coordict = {}
        for coords, orf in zip(coordinates, orfs):
            coord_list = coords.split(",")
            orf_list = orf.split(",")
            for i in range(len(coord_list)):
                if orf_list[i] not in coordict:
                    coordict[orf_list[i]] = coord_list[i]
        return coordict

    def add_information(self, output):
        # proteins = self.proteined["Protein"].tolist()
        ndf = pd.DataFrame(columns=self.proteined.columns)

        for protein in self.coordinates:
            df = self.proteined[self.proteined["Protein"].str.contains(protein)]
            proteins = df["Protein"].tolist()
            new_proteins = []
            coordinates = []
            for entry in proteins:
                new_proteins.append(protein)
                coordinates.append(self.coordinates[protein])
            df.insert(8, "entry", new_proteins)
            df.insert(5, "Genome Coordinates", coordinates)
            ndf = ndf.append(df)
        ndf = ndf.drop_duplicates()
        ndf = ndf.drop(columns='Protein')
        ndf = ndf.rename(columns={'entry': 'Protein'})
        ndf.to_csv(output, sep='\t', index=False)


class ProteinFDR(object):
    def __init__(self, folder, filetype, fdr=0.5):
        self.folder = folder
        self.filetype = filetype
        self.fdr = fdr
        self.percDir = f'{self.folder}/post_perc'
        self.filteredProtein = None
        self.__cat_protein_results()

    def __cat_protein_results(self):
        self.catProteinFile = f'{self.percDir}/{self.filetype}_cat_protein_results.txt'
        os.system(f'cat {self.folder}/Percolator/*protein_results* > {self.catProteinFile}')
        return self

    def protein_cutoff(self):
        df = pd.read_csv(self.catProteinFile, sep='\t')
        df = df[df["q-value"] != "q-value"]
        # df["q-value"] = df["q-value"].astype(float).fillna(1, 1)
        new_q = [float(i) for i in df["q-value"].tolist()]
        df = df.drop(columns="q-value")
        df.insert(3, "q-value", new_q)
        df = df[df['q-value'] <= self.fdr]
        self.filteredProtein = df

    def apply_to_psm(self):
        """ filters {filetype}_proteined based on proteins present in self.filteredProtein. Must define this attribute
        with self.protein_cutoff() before calling this function. """
        protein_df = self.filteredProtein[self.filteredProtein["ProteinId"].str.contains('ANNO') == False]
        names = protein_df["ProteinId"].tolist()
        results = pd.read_csv(f'{self.percDir}/{self.filetype}_proteined.tsv', sep='\t')
        filtered_results = pd.DataFrame(columns=results.columns)
        for i in names:
            # print(i)
            df = results[results["Protein"].str.contains(i) == True]
            filtered_results = filtered_results.append(df)
        filtered_results = filtered_results.drop_duplicates()
        self.ProteinPSMFiltered = f'{self.percDir}/{self.filetype}_psm_protein_filtered.txt'
        filtered_results.to_csv(self.ProteinPSMFiltered, sep='\t', index=False)

    def __get_utp_prots_and_coords(self):
        utps = f'{self.percDir}/{self.filetype}_utps.txt'
        utp_df = pd.read_csv(utps, sep='\t')
        prots_utp = utp_df["proteinIds"].tolist()
        coords_utp = utp_df["Genome Coordinates"].tolist()
        prots = []
        coords = []
        for prot_list, coord_list in zip(prots_utp, coords_utp):
            splat_prot = prot_list.split(",")
            splat_coords = coord_list.split(",")
            # if len(splat_prot) > 1:
            for prot, coord in zip(splat_prot, splat_coords):
                prots.append(prot)
                coords.append(coord)
            # else:
            #     prots.append(splat_prot)
            #     coords.append(splat_coords)
        return prots, coords

    def __rename_orfs(self):
        """ Removes the (...) from ORFs entries inside psm_protein_filtered.txt"""
        df = pd.read_csv(f'{self.percDir}/{self.filetype}_psm_protein_filtered.txt', sep='\t')
        names = df["Protein"].tolist()
        renamed = []
        for name in names:
            print(name)
            renamed.append(name.split("(")[0])
        df = df.drop(columns="Protein")
        df.insert(8, "Protein", renamed)
        return df

    def filter_from_utp(self):
        """
        Filters protein FDR results, checking if they are inside the filetype_utps.txt. Also adds coordinates to the
        filtered protein results.
        """
        prots, coords = self.__get_utp_prots_and_coords()
        renamed_protein_filtered_df = self.__rename_orfs()
        # df = renamed_protein_filtered_df[renamed_protein_filtered_df["Protein"].isin(prots)] # it's this piece of shit
        df = renamed_protein_filtered_df
        filtered_df = pd.DataFrame(columns=df.columns)
        filtered_df.insert(5, "Genome Coordinates", [])
        for prot, coord in zip(prots, coords):
            ndf = df[df["Protein"] == prot]
            ncoords = []
            for i in range(len(ndf["Protein"].tolist())):
                ncoords.append(coord)
            ndf.insert(5, "Genome Coordinates", ncoords)
            filtered_df = filtered_df.append(ndf)
        filtered_df = filtered_df.drop_duplicates()
        filtered_df.to_csv(f"{self.percDir}/{self.filetype}_results_01.txt", sep='\t', index=False)
        # print(filtered_df)
        return self

    def __read_db(self):
        protein_dict = {}
        db = f'{self.filetype}_database.fasta'
        records = SeqIO.parse(db, 'fasta')
        for record in records:
            if str(record.description) not in protein_dict:
                protein_dict[str(record.description)] = str(record.seq)
        return protein_dict

    def add_proteins(self, output):
        db_proteins = self.__read_db()
        results = pd.read_csv(f"{self.percDir}/{self.filetype}_results_01.txt", sep='\t')
        results = results.drop(columns="ORF Sequence")
        fixed_seqs = []
        entries = results["Protein"].tolist()
        for entry in entries:
            fixed_seqs.append(db_proteins[entry])
        results.insert(9, "db entry", fixed_seqs)
        results.to_csv(output, sep='\t', index=False)
        return self























