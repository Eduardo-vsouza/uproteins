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
from ..percolator import PercolatorData
from ..orflib import AltORF, ORF, reformat_peptide


class AltStart(PercolatorData):
    def __init__(self, pout, database):
        super().__init__(pout)
        """ 'pout' accepts a Percolator output data frame with UTPs and Genome Coordinates. """
        self.coordinates = self.get_coordinates()
        self.peptides = self.get_peptides()
        self.ids = self.get_ids()
        self.altORFs = {}
        self.totalORFs = self.__get_proteins()
        self.db = database
        self.ORFSequences = self.__get_db_proteins()
        self.includedORFs = []

    def __get_proteins(self):
        total_orfs = []
        for orfs in self.ids:
            orf_set = orfs.split(",")
            for orf in orf_set:
                total_orfs.append(orf)
        return total_orfs

    def __get_db_proteins(self):
        records = SeqIO.parse(self.db, 'fasta')
        orf_sequences = {}
        for record in records:
            if record.id in self.totalORFs:
                orf_sequences[record.id] = str(record.seq)
        self.totalORFs = []
        return orf_sequences

    def create_fake_db(self, output):
        fasta = []
        for entry in self.ORFSequences:
            fasta.append(f'>{entry}\n{self.ORFSequences[entry]}\n')
        with open(f'{output}.fasta', 'w') as fa:
            fa.writelines(fasta)

    def get_alternatives(self):
        for i in range(len(self.coordinates)):
            coord_set = self.coordinates[i].split(",")
            id_set = self.ids[i].split(",")
            peptide = self.peptides[i]
            peptide = reformat_peptide(peptide)
            if len(coord_set) > 1:
                for coord, orf in zip(coord_set, id_set):
                    coords = coord.split("-")
                    start = coords[0]
                    stop = coords[1]
                    if 'reverse' in orf:
                        start = coords[1]
                        stop = coords[0]
                        strand = 'reverse'
                    else:
                        start = coords[0]
                        stop = coords[1]
                        strand = 'forward'
                    orf_obj = ORF(seq=self.ORFSequences[orf], name=orf, start=start, end=stop, strand=strand)
                    orf_obj.MSPeptides.append(peptide)
                    if stop not in self.altORFs:
                        altorf = AltORF(strand=strand)
                        altorf.add_info(stop=stop, start=start, peptide=peptide, entry=orf)
                        altorf.add_orfs(orf_obj)
                        self.altORFs[stop] = altorf
                    else:
                        altorf = self.altORFs[stop]
                        altorf.add_orfs(orf_obj)
                        altorf.add_info(start=start, peptide=peptide, entry=orf)
                        self.altORFs[stop] = altorf
        for bla in self.altORFs:
            alto = self.altORFs[bla]
            for orfe in alto.ORFs:
                orf = alto.ORFs[orfe]
                orf.find_ms_peptides()
            included_orfs = alto.check_starts()
            for i in range(len(included_orfs)):
                self.includedORFs.append(included_orfs[i])
        return self.includedORFs


class SubsetFilter(object):
    def __init__(self, subset_df):
        self.df = pd.read_csv(subset_df, sep='\t')
        self.df = self.df.drop(self.df.columns[0], axis=1)

        self.altORFs = []

    def get_alt_info(self, genome_alts, transcriptome_alts):
        for orf in genome_alts:
            if orf not in self.altORFs:
                self.altORFs.append(orf)
        for orf in transcriptome_alts:
            if orf not in self.altORFs:
                self.altORFs.append(orf)
        return self

    def filter_alternatives(self, output):
        data = {'SpecFile': [], 'SpecID': [], 'ScanNum': [], 'FragMethod': [], 'Protein': [], 'Genome Coordinates': [],
                'Precursor': [], 'ORF Sequence': [], 'IsotopeError': [], 'PrecursorError(ppm)': [], 'Charge': [],
                'Peptide': [], 'DeNovoScore': [], 'MSGFScore': [], 'SpecEValue': [], 'EValue': []}
        new_df = pd.DataFrame(data)
        for orf in self.altORFs:
            df = self.df[self.df["Protein"].str.contains(orf.name)]
            coords = []
            names = []
            seqs = []
            if df.shape[0] != 0:
                for i in range(len(df["Protein"].tolist())):
                    coords.append(f'{orf.start}-{orf.end}')
                    name = fix_name(orf.name)
                    names.append(name)
                    seqs.append(orf.seq)
                df = df.drop(columns=['Protein', 'ORF Sequence'], axis=1)
                df.insert(4, 'Protein', names)
                df.insert(5, 'Genome Coordinates', coords)
                df.insert(7, "ORF Sequence", seqs)
                new_df = new_df.append(df)
        new_df.to_csv(f'{output}.txt', sep='\t', index=False)


class AltInspected(object):
    def __init__(self, inspected_df_subset):
        """ Filters the inspected subsets ('inspected_genome_unique', for instance) using the df generated by
        filter_alternatives() method from SubsetFilter class. """
        self.df = pd.read_csv(inspected_df_subset, sep='\t')
        self.proteins = self.df["names"].tolist()
        self.seqs = self.df["ORF Sequence"].tolist()
        self.names = []

    def get_alternative(self, filtered_orfs_df):
        df = pd.read_csv(filtered_orfs_df, sep='\t')
        names = df["Protein"].tolist()
        # names = [fix_name(name) for name in names]
        self.names = names

    def filter_inspected(self, output):
        to_write = []
        protein_list = []
        protein_seqs = []
        for i in range(len(self.proteins)):
            pro_set = self.proteins[i].split(";")
            for proteina in pro_set:
                protein=proteina
                # protein = replace_prepost(proteina)
                # protein = fix_name(protein)
                if protein not in protein_list:
                    protein_list.append(protein)
                    protein_seqs.append(self.seqs[i])
        print(protein_list)
        for i in range(len(protein_list)):
            if protein_list[i] in self.names:
                to_write.append(f'>{protein_list[i]}\n{protein_seqs[i]}\n')
        with open(f'{output}_altfiltered.fasta', 'w') as fa:
            fa.writelines(to_write)


# def replace_prepost(entry):
#     pos1 = entry.find("(")
#     protein = entry[:pos1]
#     return protein

def replace_prepost(entry):
    return entry[:-14]


def fix_name(entry):
    if 'forward' not in entry and 'reverse' not in entry:
        start = entry.rfind('_', 0, len(entry)-1)
    else:
        start = entry.find('_', 6)
    mid = entry[:start].rfind('_')
    mid = entry[:mid]
    end = entry[start:]
    final = mid[1:] + end
    return final









