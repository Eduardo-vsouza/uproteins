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

from .__helpers import HandlerError, SourceError
from .pepsearch import UProteInS


class SpectralCounting(object):
    def __init__(self, df=None, handler='file', source='MSGF'):
        """ df must be the path to a tab-separated data frame, or a pandas DataFrame instance. Handler specifies
        whether it is a 'file' or a 'object'. The 'source' attribute specifies the software used for Peptide Search.
        Currently supported ones are MSGF, uproteins. """
        self.handler = handler
        self.df = self.__check_df(df)
        self.source = source
        self.data = self.__check_source()
        self.orfs = self.data.orfs

    def __check_df(self, df):
        if self.handler == "file":
            df = pd.read_csv(df, sep="\t")
            return df
        elif self.handler == "object":
            return df
        else:
            raise HandlerError

    def __check_source(self):
        if self.source == "uproteins":
            data = UProteInS(self.df)
            return data
        else:
            raise SourceError

    def get_nsaf(self):
        orfs = self.data.get_nsaf()
        return orfs