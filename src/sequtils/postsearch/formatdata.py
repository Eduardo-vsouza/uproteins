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


class LinkData(object):
    def __init__(self, cat_msgf, peptide):
        """ cat_msgf is a table cointaining the concatenated MSGF plus results, after converting them from mzID to
        tsv. Peptide is the table containing the percolator output, after aplying the UTP identification method. """
        pd.set_option('display.max_columns', None)

        self.catDataFrame = pd.read_csv(cat_msgf, sep="\t")
        self.catDataFrame.columns = ['SpecFile', 'SpecID', 'ScanNum', 'FragMethod',	'Precursor', 'IsotopeError',
                                     'PrecursorError(ppm)', 'Charge', 'Peptide', 'Protein', 'DeNovoScore', 'MSGFScore',
                                     'SpecEValue', 'EValue']

        self.peptideDataFrame = pd.read_csv(peptide, sep="\t")
        self.peptideDataFrame = self.peptideDataFrame[self.peptideDataFrame["PSMId"] != "PSMId"]
        self.pepIds = self.peptideDataFrame["PSMId"].tolist()
        self.__get_peptide_scans()

        data = {'SpecFile': [], 'SpecID': [], 'ScanNum': [], 'FragMethod': [],	'Precursor': [], 'IsotopeError': [],
                                     'PrecursorError(ppm)': [], 'Charge': [], 'Peptide': [], 'Protein': [], 'DeNovoScore': [], 'MSGFScore': [],
                                     'SpecEValue': [], 'EValue': []}
        self.joinedDataFrame = pd.DataFrame(data)


    def __get_peptide_scans(self):
        """ retrieves the scan number for each peptide in the peptide data frame. """
        scans = []
        files = []
        for psm in self.pepIds:
            pos1 = psm.find("SII") + 4
            pos0 = pos1 - 5
            pos2 = psm[pos1:].find("_")
            scan = psm[pos1:][:pos2]
            scans.append(scan)
            file = f'{psm[:pos0]}.mzML'
            files.append(file)
        self.peptideDataFrame.insert(2, "scanNumber", scans)
        self.peptideDataFrame.insert(0, "file", files)


    def filter_msgf(self, output):
        """ filters MSGF results based on percolator output. """
        scans = self.peptideDataFrame["scanNumber"].tolist()
        files = self.peptideDataFrame["file"].tolist()

        for i in range(len(scans)):
            df = self.catDataFrame[(self.catDataFrame["ScanNum"] == int(scans[i])) & (self.catDataFrame["SpecFile"] == files[i])]
            df = df[df["Protein"].str.contains("Decoy") == False]
            df = df[df["SpecFile"] == files[i]]
            self.joinedDataFrame = self.joinedDataFrame.append(df)
        self.joinedDataFrame.to_csv(f'{output}.tsv', sep="\t", index=False)