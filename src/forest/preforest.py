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


import os
import pandas as pd


class PreFiltering(object):
    def __init__(self, pin_folder, results_04, testing=False):
        self.pinFolder = pin_folder
        self.pinFiles = os.listdir(self.pinFolder)
        self.testing = testing
        self.results = pd.read_csv(results_04, sep='\t')
        self.__check_testing()
        self.percolatorProteins = self.results["Protein"].tolist()
        self.proteins = []

    def __check_testing(self):
        if self.testing:
            self.results = self.results.head(20)

    def filter_proteins(self):
        for prot in self.percolatorProteins:
            if prot not in self.proteins:
                self.proteins.append(prot)

    def filter_pin_files(self, outdir):
        if not os.path.exists(outdir):
            os.system(f'mkdir {outdir}')
        i = 0
        for file in self.pinFiles:
            if 'fixed' in file and 'pin' in file:
                i += 1
                print(i, len(self.pinFiles), end="\r")
                df = pd.read_csv(f'{self.pinFolder}/{file}', sep='\t')
                # ndf = pd.DataFrame(columns=df.columns)
                df = df[df["Proteins"].str.contains('|'.join(self.proteins)).any(level=0)]
                # for prot in self.proteins:
                #     prot_df = df[df["Proteins"].str.contains(prot)]
                #     ndf = ndf.append(prot_df)
                df.to_csv(f'{outdir}/filtered_{file}', sep='\t', index=False)


class FeatureFishing(object):
    def __init__(self, results, pin_folder, testing=False):
        """
        :param results: should be the filetype_results_02.txt from uproteins.
        :param pin_folder: the directory containing the split pin files.
        """
        self.results = pd.read_csv(results, sep='\t')
        self.testing = testing
        self.__check_testing()
        self.specFiles = self.results["SpecFile"].tolist()
        self.scanNum = self.results["ScanNum"].tolist()
        # self.specTags = self.__add_spec_tag()  # should be added afterwards to easily identify the predictions
        self.pinFolder = pin_folder
        self.pinFiles = os.listdir(pin_folder)

    def __check_testing(self):
        if self.testing:
            # self.results = self.results.head(20)
            ...

    def __add_spec_tag(self):
        tags = []
        for i in range(len(self.scanNum)):
            tags.append(f'Spec {i+1}')
        return tags

    def add_features(self):
        chunk_size = 10**6
        i = 0
        j = 0
        for file in self.pinFiles:
            print(j, len(self.pinFiles), end='\r')
            if 'fixed' in file and 'pin' in file:
                # if 'pin' in file:
                chunk = pd.read_csv(f'{self.pinFolder}/{file}', sep='\t')
                chunk = chunk[chunk["ScanNr"].isin(self.scanNum)]
                # for chunk in pd.read_csv(f'{self.pinFolder}/{file}', sep='\t', chunksize=chunk_size):
                if i == 0:
                    self.dataWithFeatures = pd.DataFrame(columns=chunk.columns)
                    i += 1
                for spec, scan in zip(self.specFiles, self.scanNum):
                    self.__filter_by_scan(chunk, spec, scan)
                # print(j)
                j += 1
        return self

    def __filter_by_scan(self, chunk, spec, scan):
        df = chunk[(chunk["SpecId"].str.contains(spec[:-5])) & (chunk["ScanNr"] == scan)]
        self.dataWithFeatures = self.dataWithFeatures.append(df)
        return self

    def save_table(self, output):
        self.dataWithFeatures.to_csv(output, sep='\t', index=False)
        return self


if __name__ == '__main__':
    # data = FeatureFishing(results='for_predicting/genome_results_02.txt', pin_folder='pin_files/mtb')  ## the ol' working code
    # data.add_features()
    # data.save_table('for_predicting/Transcriptome/transcriptome_results_with_features.txt')
    # genome_prefiltering = PreFiltering(pin_folder='genome_fixed_pin_files',
    #                                    results_04='uproteins_results_1608/genome_results_04.txt')
    # genome_prefiltering.filter_proteins()
    # genome_prefiltering.filter_pin_files('genome_fixed_filtered_pin_files')


    # data = FeatureFishing(results='uproteins_results_1608/genome_results_04.txt', pin_folder='genome_fixed_filtered_pin_files')
    # data.add_features()
    # data.save_table('for_predicting_1608/genome_results_with_features.txt')

    # data_3008 = FeatureFishing(results='uproteins_results_3008/genome_results_04.txt', pin_folder='genome_fixed_filtered_pin_files')
    # data_3008.add_features()
    # data_3008.save_table('uproteins_results_3008/for_predicting/genome_results_with_features.txt')

    transcriptome_3008 = FeatureFishing(results='uproteins_results_3008/transcriptome_results_04.txt',
                                        pin_folder='transcriptome_fixed_pin_files')
    transcriptome_3008.add_features()
    transcriptome_3008.save_table('uproteins_results_3008/for_predicting/transcriptome_results_with_features.txt')


