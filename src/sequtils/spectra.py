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