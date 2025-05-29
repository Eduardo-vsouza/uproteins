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


class SequenceFinder(object):
    def __init__(self, df, fasta_db):
        self.df = pd.read_csv(df, sep="\t")
        self.df = self.df[self.df["Protein"].str.contains("contaminant", regex=False) == False]
        self.df = self.df[self.df["Protein"].str.contains("lcl|", regex=False) == False]
        self.df = self.df[self.df["Protein"].str.contains("decoy", regex=False) == False]
        self.df = self.df[self.df["Protein"].str.contains("sp|", regex=False) == False]
        self.proteins = self.df["Protein"].tolist()
        self.fasta = fasta_db
        self.proteinDict = self.__get_db_proteins()

    def __get_db_proteins(self):
        protein_dict = {}
        records = SeqIO.parse(self.fasta, 'fasta')
        for record in records:
            if record.description not in protein_dict:
                protein_dict[record.description.split(" ")[0]] = record.seq
        return protein_dict

    def df_proteins(self):
        protein_list = []
        for i in self.proteins:
            protein_set = i.split(";")
            seq_set = ""
            for protein in protein_set:
                if 'decoy' not in protein:
                    # print(protein)
                    # if ')' not in protein:
                    pos = protein.rfind("(pre")
                    fixed = protein[:pos]
                    # fixed = protein
                    seq = self.proteinDict[fixed]
                    # print(seq)
                    if len(seq_set) > 0:
                        seq_set += f",{seq}"
                    else:
                        seq_set += seq
            protein_list.append(seq_set)
        self.df.insert(5, "ORF Sequence", protein_list)
        return self

    def save(self, output):
        self.df.to_csv(f'{output}.tsv', sep='\t', index=False)
        return self


class Subsets(object):
    def __init__(self, genome, transcriptome):
        self.genomeDataFrame = pd.read_csv(genome, sep='\t')
        self.transcriptomeDataFrame = pd.read_csv(transcriptome, sep='\t')

        self.genomeORFs = self.genomeDataFrame["ORF Sequence"].tolist()
        self.transcriptomeORFs = self.transcriptomeDataFrame["ORF Sequence"].tolist()

        self.bORFs = []
        self.gORFs = []
        self.tORFs = []

    def get_orfs(self):
        genome = []
        transcriptome = []
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

        # self.gORFs = [orf for orf in genome if orf not in transcriptome]
        # self.bORFs = [orf for orf in genome if orf in transcriptome]
        # self.tORFs = [orf for orf in transcriptome if orf not in genome]
        return self

    def count_orfs(self):
        output = f'bORFs: {len(self.bORFs)}\n' \
            f'gORFs: {len(self.gORFs)}\n' \
            f'tORFs: {len(self.tORFs)}'
        print(output)
        return self

    def create_venn(self):
        venn2(subsets=(len(self.tORFs), len(self.gORFs), len(self.bORFs)), set_labels=('tORFs', 'gORFs'))
        plt.title('ORF subsets')
        plt.show()

    def adapt_manual_inspect(self, inspected_df):
        ''' only for our smeg analysis '''
        df = pd.read_csv(inspected_df, sep='\t')
        seqs = df["ORF Sequence"].tolist()
        gorfs = []
        torfs = []
        borfs = []

        for orf in seqs:
            if orf in self.gORFs:
                gorfs.append(orf)
            elif orf in self.tORFs:
                torfs.append(orf)
            elif orf in self.bORFs:
                borfs.append(orf)
        venn2(subsets=(len(torfs), len(gorfs), len(borfs)), set_labels=('tORFs', 'gORFs'))
        plt.title('ORF subsets after manual inspection')
        plt.show()


class FastaSubsetter(object):
    def __init__(self, subset_df):
        self.df = pd.read_csv(subset_df, sep='\t')
        self.df = self.df[self.df["ORF Sequence"].str.len() <= 100]
        self.proteins = self.df["Protein"].tolist()
        self.seqs = self.df["ORF Sequence"].tolist()

    def __get_seqs(self):
        proteins = []
        to_write = []
        for i in range(len(self.proteins)):
            pro_set = self.proteins[i].split(";")
            for protein in pro_set:
                if protein not in proteins:
                    proteins.append(protein)
                    to_write.append(f'>{protein}\n{self.seqs[i]}\n')
        return to_write

    def create_fasta(self, output):
        fasta = self.__get_seqs()
        with open(f'{output}.fasta', 'w') as fa:
            fa.writelines(fasta)
        return self

    def filter_smorfs(self, output):
        df = self.df[self.df["ORF Sequence"].str.len() <= 100]
        df.to_csv(f'{output}_smorfs.tsv', sep='\t', index=False)


class PeptideSubsets(object):
    def __init__(self, genome, transcriptome):
        self.genomeDataFrame = pd.read_csv(genome, sep='\t')
        self.transcriptomeDataFrame = pd.read_csv(transcriptome, sep='\t')
        # self.bothDataFrame = pd.read_csv(both, sep='\t')
        self.genomePeptides = self.genomeDataFrame["peptide"].tolist()
        self.transcriptomePeptides = self.transcriptomeDataFrame["peptide"].tolist()

        self.fixedGenomePeptides = self.__get_peptides(self.genomePeptides)
        self.fixedTranscriptomePeptides = self.__get_peptides(self.transcriptomePeptides)

        self.gPeptides = []
        self.tPeptides = []
        self.bPeptides = []

    def __get_peptides(self, pepset):
        def reformat_peptide(peptide):
            pep = peptide.replace(".", "")
            pep = pep.replace("-", "")
            while '[' in pep:
                pos1 = pep.find("[")
                pos2 = pep.find("]") + 1
                to_replace = pep[pos1:pos2]
                pep = pep.replace(to_replace, "")
            return pep

        peptides = []
        for pep in pepset:
            peptide = reformat_peptide(pep)
            if peptide not in peptides:
                peptides.append(peptide)
        return peptides

    def shared(self):
        for gpep in self.fixedGenomePeptides:
            if gpep in self.fixedTranscriptomePeptides and gpep not in self.bPeptides:
                self.bPeptides.append(gpep)
            elif gpep not in self.fixedTranscriptomePeptides and gpep not in self.gPeptides:
                self.gPeptides.append(gpep)
        for tpep in self.fixedTranscriptomePeptides:
            # if tpep in self.fixedGenomePeptides and tpep not in self.bPeptides:
            #     self.bPeptides.append(tpep)
            if tpep not in self.fixedGenomePeptides and tpep not in self.tPeptides:
                self.tPeptides.append(tpep)
        print(self.tPeptides)
        print(self.gPeptides)
        print(self.bPeptides)
        return self

    def venn(self):
        venn2(subsets=(len(self.tPeptides), len(self.gPeptides), len(self.gPeptides)), set_labels=('tPeptides', 'gPeptides'))
        plt.show()


class CollectionSubsets(object):
    def __init__(self, ncbi, uniprot, mycobrowser):
        self.ncbi = self.__get_records(ncbi)
        self.uniprot = self.__get_records(uniprot)
        self.myco = self.__get_records(mycobrowser)

        self.ncbiSeqs = self.__get_seqs(self.ncbi)
        self.uniprotSeqs = self.__get_seqs(self.uniprot)
        self.mycoSeqs = self.__get_seqs(self.myco)

        self.ncbiUnique = []
        self.uniprotUnique = []
        self.mycoUnique = []

    @staticmethod
    def __get_records(db):
        records = SeqIO.parse(db, 'fasta')
        return records

    @staticmethod
    def __get_seqs(db):
        seqs = []
        for record in db:
            if str(record.seq) not in seqs:
                seqs.append(str(record.seq))
        return seqs

    def unique_sequences(self):
        for seq in self.ncbiSeqs:
            if seq not in self.uniprotSeqs and seq not in self.mycoSeqs:
                self.ncbiUnique.append(seq)
        for seq in self.uniprotSeqs:
            if seq not in self.ncbiSeqs and seq not in self.mycoSeqs:
                self.uniprotUnique.append(seq)
        for seq in self.mycoSeqs:
            if seq not in self.ncbiSeqs and seq not in self.uniprotSeqs:
                self.mycoUnique.append(seq)

        print(f'Unique sequences in each database\nUniprot: {len(set(self.uniprotUnique))}\nNCBI: {len(set(self.ncbiUnique))}\n'
              f'Mycobrowser: {len(set(self.mycoUnique))}')
        return self

    def join_collections(self):
        joint = self.uniprotUnique + self.ncbiUnique + self.mycoUnique
        return joint

    def filter_data(self, proteined_df, **kwargs):
        joint = self.join_collections()
        output = ""
        if kwargs.get("output"):
            output = kwargs.get("output")
        df = pd.read_csv(proteined_df, sep='\t')
        print(df.shape)
        df = df[df["ORF Sequence"].isin(joint) == False]
        print(df.shape)
        if kwargs.get("save"):
            df.to_csv(f'{output}.xls', sep='\t', index=False)
        return self

    def filter_inspected(self, inspected, filtered_df, **kwargs):
        output = ""
        if kwargs.get("output"):
            output = kwargs.get("output")
        insp = pd.read_csv(inspected, sep='\t')
        filtered = pd.read_csv(filtered_df, sep='\t')
        scans = set(insp["Peptide"].tolist())
        filtered = filtered[filtered["Peptide"].isin(scans) == False]
        if kwargs.get('save'):
            filtered.to_csv(f'{output}.xls', sep='\t', index=False)
        return self



