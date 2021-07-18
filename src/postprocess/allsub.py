import os
import pandas as pd


class AllSub(object):
    def __init__(self, folder, filetype):
        self.folder = folder
        self.filetype = filetype
        self.pinFile = pd.read_csv(f'{self.folder}/Percolator/{self.folder}_pin.txt', sep='\t')

    def remove_irrelevant(self):
        """ Removes from the pin files any peptide that matches an annotated protein, leaving only those relevant
        to the hypothesis of the pipeline. """
        self.pinFile = self.pinFile[self.pinFile["Proteins"].str.contains('ANNO') == False]
        os.system(f'mv {self.folder}/Percolator/{self.folder}_pin.txt {self.folder}/Percolator/{self.folder}_all_pin.txt')

    def save(self):
        self.pinFile.to_csv(f'{self.folder}/Percolator/{self.folder}_pin.txt', sep='\t', index=False)


