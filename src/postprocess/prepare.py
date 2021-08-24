import os
import pandas as pd


class ForestCleaning(object):
    def __init__(self, folder, filetype):
        """
        Class for cleaning data before using the random forest model for filtering spectra.
        """
        self.folder = folder
        self.filetype = filetype
        self.percDir = f'{self.folder}/post_perc'
        self.input = pd.read_csv(f'{self.percDir}/{self.filetype}_psm_protein_filtered.txt', sep='\t')
        self.pinFile = f'{folder}/{self.folder}_pin.txt'
        self.output = None
        self.MLFolder = f'{self.folder}/rf_predictions'

    def fetch_spectra(self):
        """
        Grabs the necessary information to run the model from the pin files.
        """
        files = self.input["SpecFile"].tolist()
        spec_num = self.input["SpecID"].tolist()
        chunk_size = 10**6
        j = 0
        identifiers = []  # holds spectrum information across multiple files
        for i in range(len(files)):
            for chunk in pd.read_csv(self.pinFile, sep='\t', chunksize=chunk_size):
                df = chunk[chunk["SpecId"].str.contains(files[i][:-5])]  # ignores the '.mzML' suffix when comparing
                df = df[df["ScanNr"] == spec_num[i]]
                if j == 0:
                    new_df = pd.DataFrame(columns=chunk.columns)
                    new_df.insert(2, "Identifier", [])
                identifiers.append(f'spec {j}')
                ids, j = self.__add_identifiers(df, j)
                df.insert(2, "Identifier", ids)
                new_df = new_df.append(df)
        self.output = new_df

    def __create_directory(self):
        if not os.path.exists(self.MLFolder):
            os.system(f'mkdir {self.MLFolder}')
        return self

    def save_table(self):
        self.__create_directory()
        if self.output is not None:
            self.output.to_csv(f'{self.MLFolder}/{self.filetype}_model_input.txt', sep='\t', index=False)
        return self

    @staticmethod
    def __add_identifiers(df, j):
        scans = df["ScanNr"].tolist()
        ids = []
        for i in range(len(scans)):
            ids.append(f'Spec {j}')
            j += 1
        return ids, j

