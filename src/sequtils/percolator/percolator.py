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


class PercolatorData(object):
    def __init__(self, pout):
        self.df = pd.read_csv(pout, sep='\t')

    def get_coordinates(self):
        coords = self.df["Genome Coordinates"].tolist()
        return coords

    def get_peptides(self):
        peps = self.df["peptide"].tolist()
        return peps

    def get_ids(self):
        entries = self.df["proteinIds"].tolist()
        return entries
