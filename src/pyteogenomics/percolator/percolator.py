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
