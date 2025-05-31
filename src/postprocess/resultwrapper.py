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

from ..sequtils.utilities import check_dir


class ResultsWrapper(object):
    def __init__(self, df, folder, filetype):
        self.df = pd.read_csv(df, sep='\t')
        self.folder = folder
        self.filetype = filetype
        self.resultsFolder = f'{folder}/Results'
        check_dir(self.resultsFolder)

    def __get_orf_names(self):
        new_names = []
        not_extended = self.df["Protein"].tolist()
        for orf in not_extended:
            new = '_'.join(orf.split("_")[:-2]).replace('__', '_')
            new_names.append(new)
        return new_names

    def reformat(self, pre_validation=True):
        new_names = self.__get_orf_names()
        to_remove = ['ORF Sequence', 'Protein', 'Fixed Peptides', "Extended ORF", "Genome Coordinates"]
        # if not pre_validation:
        #     to_remove.append("Renamed Files")
        self.df = self.df.drop(columns=to_remove)
        self.df = self.df.rename(columns={'Extended Sequence': 'Protein sequence', "Extended Coordinates":
                                          "Genome Coordinates"})
        self.df.insert(1, "ORF name", new_names)
        if pre_validation:
            pattern = "pre_validation"
        else:
            pattern = "post_validation"
        self.df.to_csv(f'{self.resultsFolder}/{self.filetype}_{pattern}_results.txt', sep='\t', index=False)
