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
import sys

from ..forest import ProteinFixer, PreFiltering, FeatureFishing, SpectralForest, SpectrumMiner
from ..postprocess import ResultsWrapper


class ValidatePipeline(object):
    def __init__(self, args, testing=False):
        self.args = args
        self.pypath = sys.path[0]
        self.testing = testing

    def validate_genome(self):
        if not os.path.exists("Genome/Percolator/for_predicting"):
            os.mkdir('Genome/Percolator/for_predicting')

        genome_pin = ProteinFixer(pin_folder='Genome/Percolator')
        genome_pin.fix_files(outdir='Genome/Percolator')
        genome_prefiltering = PreFiltering(pin_folder='Genome/Percolator',
                                           results_04='Genome/post_perc/genome_results_04.txt', testing=self.testing)
        genome_prefiltering.filter_proteins()
        genome_prefiltering.filter_pin_files('Genome/Percolator/for_predicting')
        genome_fishing = FeatureFishing(results='Genome/post_perc/genome_results_04.txt',
                                        pin_folder='Genome/Percolator/for_predicting',
                                        testing=self.testing)
        genome_fishing.add_features()
        genome_fishing.save_table('Genome/Percolator/for_predicting/genome_results_with_features.txt')

        genome_sf = SpectralForest(model_pickle=f'{self.pypath}/src/model/model.pkl',
                                   results='Genome/Percolator/for_predicting/genome_results_with_features.txt')
        genome_sf.predict()
        if not os.path.exists('Genome/predicted'):
            os.system('mkdir Genome/predicted')
        genome_sf.save('Genome/predicted/genome_predicted.txt')
        genome_mining = SpectrumMiner(predicted='Genome/predicted/genome_predicted.txt',
                                      results='Genome/post_perc/genome_results_04.txt', testing=self.testing)
        genome_mining.filter_results('Genome/post_perc/genome_results_05.txt')
        genome_mining.add_lc_spectra_for_hc_orf(filetype='genome', folder='Genome')
        genome_results = ResultsWrapper(df='Genome/post_perc/genome_results_05.txt', folder='Genome', filetype='genome')
        genome_results.reformat(pre_validation=False)

    def validate_transcriptome(self):
        if self.args.transcriptome:
            if not os.path.exists('Transcriptome/Percolator/for_predicting'):
                os.system('mkdir Transcriptome/Percolator/for_predicting')
            transcriptome_pin = ProteinFixer(pin_folder='Transcriptome/Percolator')
            transcriptome_pin.fix_files(outdir='Transcriptome/Percolator')
            transcriptome_prefiltering = PreFiltering(pin_folder='Transcriptome/Percolator',
                                                      results_04='Transcriptome/post_perc/transcriptome_results_04.txt',
                                                      testing=self.testing)
            transcriptome_prefiltering.filter_proteins()

            transcriptome_prefiltering.filter_pin_files('Transcriptome/Percolator/for_predicting')
            transcriptome_fishing = FeatureFishing(results='Transcriptome/post_perc/transcriptome_results_04.txt',
                                                   pin_folder='Transcriptome/Percolator/for_predicting', testing=self.testing)
            transcriptome_fishing.add_features()
            transcriptome_fishing.save_table(
                'Transcriptome/Percolator/for_predicting/transcriptome_results_with_features.txt')

            transcriptome_sf = SpectralForest(model_pickle=f'{self.pypath}/src/model/model.pkl',
                                              results='Transcriptome/Percolator/for_predicting/transcriptome_results_with_features.txt')
            transcriptome_sf.predict()
            if not os.path.exists('Transcriptome/predicted'):
                os.system('mkdir Transcriptome/predicted')
            transcriptome_sf.save('Transcriptome/predicted/transcriptome_predicted.txt')
            transcriptome_mining = SpectrumMiner(predicted='Transcriptome/predicted/transcriptome_predicted.txt',
                                                 results='Transcriptome/post_perc/transcriptome_results_04.txt',
                                                 testing=self.testing)
            transcriptome_mining.filter_results('Transcriptome/post_perc/transcriptome_results_05.txt')
            transcriptome_mining.add_lc_spectra_for_hc_orf(filetype='transcriptome', folder='Transcriptome')
            transcriptome_results = ResultsWrapper(df='Transcriptome/post_perc/transcriptome_results_05.txt',
                                                   folder='Transcriptome', filetype='transcriptome')
            transcriptome_results.reformat(pre_validation=False)
