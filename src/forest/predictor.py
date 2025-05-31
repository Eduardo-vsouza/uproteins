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


import pickle
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import MinMaxScaler
import os


class SpectralForest(object):
    def __init__(self, model_pickle, results):
        self.df = pd.read_csv(results, sep='\t')
        # self.df = self.df[self.df["SpecFile"] != "DefaultDirection"]
        self.df = self.df[self.df["ExpMass"] != "ExpMass"]
        self.df = self.df[self.df["Proteins"].str.contains("Decoy") == False]
        self.pickleFile = open(model_pickle, 'rb')
        self.model = pickle.load(self.pickleFile)

        self.properties = list(self.df.columns)
        # self.properties.remove('label_spectra')
        self.__remove_features()
        self.x = self.df[self.properties]
        # self.y = self.df["label_spectra"]
        self.__scale_data()

    def __remove_features(self):
        to_rem = ['SpecFile', 'SpecId', 'Label', 'ScanNr', 'Peptide', 'Proteins', 'proteinIds',
                  'enzN',	'enzC', 'enzInt', 'StdevErrorTop7', 'MeanErrorTop7']
        for i in range(20):
            to_rem.append(f'Charge{i}')
        for i in to_rem:
            try:
                self.properties.remove(i)
            except:
                pass

    def __scale_data(self):
        scaler = MinMaxScaler()
        scaler.fit(self.x)
        self.x = scaler.fit_transform(self.x)

    def predict(self):
        predicted = self.model.predict(self.x)
        self.__insert_predicted(predicted)
        return predicted

    def __insert_predicted(self, predicted):
        self.df.insert(2, 'Prediction', predicted)
        self.df = self.df[self.df["Prediction"] == 1]
        return self

    def save(self, output):
        self.df.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    # model = SpectralForest(model_pickle='semisupervised_rf_30_04.pickle', results='for_predicting/genome_results_with_features.txt')
    # model.predict()
    # model.save('for_predicting/genome_predicted.txt')

    # direc = 'pin_files/mtb/more_splits'
    # direc = 'pin_files/mtb/transcriptome'
    # files = os.listdir(direc)
    # for file in files:
    #     if 'specfiles' in file:
    #         print(file)
    #         model = SpectralForest(model_pickle='semisupervised_rf_30_04.pickle', results=f'{direc}/{file}')
    #         # model = SpectralForest(model_pickle='rf_proteins_30_04.pickle', results=f'{direc}/{file}')
    #
    #         model.predict()
    #         model.save(f'for_predicting/Transcriptome/{file}_predicted_semi_supervised.txt')
    # model = SpectralForest(model_pickle='semisupervised_rf_30_04.pickle', results='for_predicting/Genome/genome_results_with_features.txt')
    # model.predict()
    # model.save('predicted/genome_predicted_2906.txt')
    exo_model = SpectralForest(model_pickle='semisupervised_rf_30_04.pickle', results='pin_files/exosomes_mtb/mtb_pin.txt_specfiles_full.txt')
    exo_model.predict()
    exo_model.save('predicted/exosomes_mtb_predicted.txt')
    """  python3 ~/Programas/uproteins_0908/ProteInS.py testing --skip_assembly TRUE --skip_db TRUE --skip_ms TRUE --skip_postms TRUE --outdir testing_validate/
"""
