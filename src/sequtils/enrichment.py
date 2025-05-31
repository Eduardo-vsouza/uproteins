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
import matplotlib.pyplot as plt
import seaborn as sns


class Enrichment(object):
    def __init__(self, rna_df, genome_df, both_df):
        self.RNADataFrame = pd.read_csv(rna_df, sep='\t')
        self.RNADataFrame = self.RNADataFrame.drop(self.RNADataFrame.columns[0], axis=1)
        self.RNAFiles = self.RNADataFrame["SpecFile"].tolist()
        self.RNASeqs = self.RNADataFrame["ORF Sequence"].tolist()

        self.DNADataFrame = pd.read_csv(genome_df, sep='\t')
        self.DNADataFrame = self.DNADataFrame.drop(self.DNADataFrame.columns[0], axis=1)
        self.DNAFiles = self.DNADataFrame["SpecFile"].tolist()
        self.DNASeqs = self.DNADataFrame["ORF Sequence"].tolist()

        self.bothDataFrame = pd.read_csv(both_df, sep='\t')
        self.bothDataFrame = self.bothDataFrame.drop(self.bothDataFrame.columns[0], axis=1)
        self.bothFiles = self.bothDataFrame["SpecFile"].tolist()
        self.bothSeqs = self.bothDataFrame["ORF Sequence"].tolist()

    def methods_df(self):
        rna = self.__orfs_by_method(self.RNAFiles, self.RNASeqs)
        dna = self.__orfs_by_method(self.DNAFiles, self.DNASeqs)
        both = self.__orfs_by_method(self.bothFiles, self.bothSeqs)

        def count_orfs(subset):
            sub_dict = {}
            for meth in subset:
                orfs = subset[meth].split(",")
                orf_numb = len(orfs)
                sub_dict[meth] = orf_numb
            return sub_dict

        rna_dict = count_orfs(rna)
        dna_dict = count_orfs(dna)
        both_dict = count_orfs(both)

        orfs = []
        meth = []
        subset = []

        def df_info(sub_dict, subset_name):
            for enrich in sub_dict:
                orfs.append(sub_dict[enrich])
                meth.append(enrich)
                subset.append(subset_name)

        df_info(rna_dict, 'transcriptome')
        df_info(dna_dict, 'genome')
        df_info(both_dict, 'both')

        data = {'ORFs': orfs, 'Enrichment Method': meth, 'Subset': subset}
        self.methodsDataFrame = pd.DataFrame(data)
        return self

    def plot_bars(self):
        x = "Enrichment Method"
        y = 'ORFs'
        hue = 'Subset'
        ax = sns.barplot(x=x, y=y, hue=hue, data=self.methodsDataFrame)
        plt.title('ORFs identified using different enrichment methods')
        plt.show()

    def save_df(self, output):
        rna = self.methodsDataFrame[self.methodsDataFrame["Subset"] == 'transcriptome']
        dna = self.methodsDataFrame[self.methodsDataFrame["Subset"] == "genome"]
        both = self.methodsDataFrame[self.methodsDataFrame["Subset"] == "both"]
        total = self.methodsDataFrame[""]
        rna.to_csv(f'{output}_rna.txt', sep='\t', index=False)
        dna.to_csv(f'{output}_dna.txt', sep='\t', index=False)
        both.to_csv(f'{output}_both.txt', sep='\t', index=False)
        return self


    @staticmethod
    def __orfs_by_method(files, seqs):
        methods = {}
        for i in range(len(files)):
            pos1 = files[i].find("_") + 1
            pos2 = files[i].find(".") - 2
            file = files[i][pos1:pos2]
            if file not in methods:
                methods[file] = seqs[i]
            else:
                if seqs[i] not in methods[file].split(","):
                    methods[file] += f',{seqs[i]}'
        return methods

