# Copyright © 2021-2025 Eduardo Vieira de Souza
# Copyright © 2021-2025 Adriana Canedo
# Copyright © 2021-2025 Cristiano Valim Bizarro
#
# This file is part of uProteInS.
#
# uProteInS is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# uProteInS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# uProteInS. If not, see <https://www.gnu.org/licenses/>.


import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib_venn import venn2


class DecoyVoid(object):
    def __init__(self, genome_linked, transcriptome_linked, genome_db, transcriptome_db):
        self.genomeDataFrame = pd.read_csv(genome_linked, sep='\t')
        self.genomeDataFrame = self.genomeDataFrame[self.genomeDataFrame["ORF Sequence"].str.len() <= 100]
        self.transcriptomeDataFrame = pd.read_csv(transcriptome_linked, sep='\t')
        self.transcriptomeDataFrame = self.transcriptomeDataFrame[self.transcriptomeDataFrame["ORF Sequence"].str.len() <= 100]

        self.genomeORFs = self.genomeDataFrame["ORF Sequence"].tolist()

        self.transcriptomeORFs = self.transcriptomeDataFrame["ORF Sequence"].tolist()

        self.genomeDB = genome_db
        self.transcriptomeDB = transcriptome_db

        self.genomeDict = {}
        self.transcriptomeDict = {}
        self.__add_orfs()

        self.tORFs = []
        self.gORFs = []
        self.bORFs = []
        self.__add_subsets()

        self.data = {'SpecFile': [], 'SpecID': [], 'ScanNum': [], 'FragMethod': [],	'Precursor': [], 'IsotopeError': [],
                                     'PrecursorError(ppm)': [], 'Charge': [], 'Peptide': [], 'Protein': [], 'DeNovoScore': [], 'MSGFScore': [],
                                     'SpecEValue': [], 'EValue': []}


        self.genomeUnique = pd.DataFrame(self.data)
        self.trancriptomeUnique = pd.DataFrame(self.data)
        self.both = pd.DataFrame(self.data)

    def count_peptides(self):
        peps = []
        gpeps = self.genomeDataFrame["Peptide"].tolist()
        rpeps = self.transcriptomeDataFrame["Peptide"].tolist()
        for pep in gpeps:
            if pep not in peps:
                peps.append(pep)
        for pep in rpeps:
            if pep not in peps:
                peps.append(pep)
        print(len(peps))

    def __add_orfs(self):
        g_records = SeqIO.parse(self.genomeDB, 'fasta')
        t_records = SeqIO.parse(self.transcriptomeDB, 'fasta')

        for record in g_records:
            if record.seq not in self.genomeDict:
                self.genomeDict[record.seq] = record.id
        for record in t_records:
            if record.seq not in self.transcriptomeDict:
                self.transcriptomeDict[record.seq] = record.id
        return self

    def __add_subsets(self):
        for orf in self.genomeORFs:
            if orf not in self.gORFs:
                if orf not in self.transcriptomeORFs:
                    self.gORFs.append(orf)

                else:
                    if orf not in self.bORFs:
                        self.bORFs.append(orf)
        for orf in self.transcriptomeORFs:
            if orf not in self.tORFs:
                if orf not in self.genomeORFs:
                    self.tORFs.append(orf)
                else:
                    if orf not in self.bORFs:
                        self.bORFs.append(orf)

    def check_decoy(self):
        for orf in self.gORFs:
            if orf in self.transcriptomeDict:
                if orf not in self.bORFs:
                    self.bORFs.append(orf)
        for orf in self.tORFs:
            if orf in self.genomeDict:
                if orf not in self.bORFs:
                    self.bORFs.append(orf)
        self.__remove_orfs()

        output = f'bORFs: {len(self.bORFs)}\n' \
            f'gORFs: {len(self.gORFs)}\n' \
            f'tORFs: {len(self.tORFs)}'

        print(output)
        return self

    def separate_subsets(self, folder):
        gorfs = self.genomeDataFrame[self.genomeDataFrame["ORF Sequence"].isin(self.gORFs)]
        torfs = self.transcriptomeDataFrame[self.transcriptomeDataFrame["ORF Sequence"].isin(self.tORFs)]
        borfs = self.genomeDataFrame[self.genomeDataFrame["ORF Sequence"].isin(self.bORFs)]
        t_to_borfs = self.transcriptomeDataFrame[self.transcriptomeDataFrame["ORF Sequence"].isin(self.bORFs)]
        borfs = borfs.append(t_to_borfs)

        gorfs.to_csv(f'{folder}genome_unique.tsv', sep='\t')
        torfs.to_csv(f'{folder}transcriptome_unique.tsv', sep='\t')
        borfs.to_csv(f'{folder}both.tsv', sep='\t')

    def subset_single_df(self, output):
        gorfs = self.genomeDataFrame[self.genomeDataFrame["ORF Sequence"].isin(self.gORFs)]
        gorf_type = ["Genome" for i in range(len(gorfs["ORF Sequence"].tolist()))]

        torfs = self.transcriptomeDataFrame[self.transcriptomeDataFrame["ORF Sequence"].isin(self.tORFs)]
        torf_type = ["Transcriptome" for i in range(len(torfs["ORF Sequence"].tolist()))]

        borfs = self.genomeDataFrame[self.genomeDataFrame["ORF Sequence"].isin(self.bORFs)]

        t_to_borfs = self.transcriptomeDataFrame[self.transcriptomeDataFrame["ORF Sequence"].isin(self.bORFs)]
        borfs = borfs.append(t_to_borfs)
        borf_type = ["Both" for i in range(len(borfs["ORF Sequence"].tolist()))]

        gorfs.insert(10, "Subset", gorf_type)
        torfs.insert(10, "Subset", torf_type)
        borfs.insert(10, "Subset", borf_type)
        df = gorfs.append(torfs)
        df = df.append(borfs)
        df.to_csv(f'{output}.xls', sep='\t', index=False)

    def __remove_orfs(self):
        """ called after adding from decoy to both"""
        for orf in self.bORFs:
            if orf in self.gORFs:
                self.gORFs.remove(orf)
            if orf in self.tORFs:
                self.tORFs.remove(orf)

    def create_venn(self):
        venn2(subsets=(len(self.tORFs), len(self.gORFs), len(self.bORFs)), set_labels=('tORFs', 'gORFs'))
        plt.title('ORF subsets after applying decoy method')
        plt.show()

    def adapt_manual_inspect(self, inspected_df):
        ''' only for our smeg analysis '''
        df = pd.read_csv(inspected_df, sep='\t')
        df = df[(df["Result"] == "High") | (df["Result"] == "high")]
        df = df[df["ORF Sequence"].str.len() <= 100]
        seqs = df["ORF Sequence"].tolist()
        total = []
        gorfs = []
        torfs = []
        borfs = []

        for orf in seqs:
            if orf not in total:
                total.append(orf)
            # if orf in self.gORFs and orf not in gorfs:
            #     gorfs.append(orf)
            # elif orf in self.tORFs and orf not in torfs:
            #     torfs.append(orf)
            # elif orf in self.bORFs and orf not in borfs:
            #     borfs.append(orf)
        print(len(total))
        # venn2(subsets=(len(torfs), len(gorfs), len(borfs)), set_labels=('tORFs', 'gORFs'))
        # plt.title('ORF subsets after manual inspection')
        # plt.show()

    def save_manual_inspect(self, inspected_df, output):
        ''' only for our smeg analysis '''
        df = pd.read_csv(inspected_df, sep='\t')
        df = df[(df["Result"] == "High") | (df["Result"] == "high")]
        torf = df[df["ORF Sequence"].isin(self.tORFs)]
        gorf = df[df["ORF Sequence"].isin(self.gORFs)]
        borf = df[df["ORF Sequence"].isin(self.bORFs)]
        torf.to_csv(f'{output}_transcriptome_unique.xls', sep='\t', index=False)
        gorf.to_csv(f'{output}_genome_unique.xls', sep='\t', index=False)
        borf.to_csv(f'{output}_both.xls', sep='\t', index=False)

