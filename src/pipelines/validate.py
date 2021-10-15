from ..forest import ProteinFixer, PreFiltering, FeatureFishing, SpectralForest
import os
import sys


class ValidatePipeline(object):
    def __init__(self, args):
        self.args = args
        self.pypath = sys.path[0]

    def validate_genome(self):
        genome_pin = ProteinFixer(pin_folder='Genome/Percolator')
        genome_pin.fix_files(outdir='Genome/Percolator')
        genome_prefiltering = PreFiltering(pin_folder='Genome/Percolator',
                                           results_04='Genome/post_perc/genome_results_04.txt')
        genome_prefiltering.filter_proteins()
        if not os.path.exists("Genome/Percolator/for_predicting"):
            os.system('mkdir Genome/Percolator/for_predicting')
        genome_prefiltering.filter_pin_files('Genome/Percolator/for_predicting')
        genome_fishing = FeatureFishing(results='Genome/post_perc/genome_results_04.txt',
                                        pin_folder='Genome/Percolator/for_predicting')
        genome_fishing.add_features()
        genome_fishing.save_table('Genome/Percolator/for_predicting/genome_results_with_features.txt')

        genome_sf = SpectralForest(model_pickle=f'{self.pypath}/src/model/model.pkl',
                                   results='Genome/Percolator/for_predicting/genome_results_with_features.txt')
        genome_sf.predict()
        if not os.path.exists('Genome/predicted'):
            os.system('mkdir Genome/predicted')
        genome_sf.save('Genome/predicted/genome_predicted.txt')

    def validate_transcriptome(self):
        if self.args.Transcriptome == "YES":
            # transcriptome_pin = ProteinFixer(pin_folder='Transcriptome/Percolator')
            # transcriptome_pin.fix_files(outdir='Transcriptome/Percolator')
            # transcriptome_prefiltering = PreFiltering(pin_folder='Transcriptome/Percolator',
            #                                           results_04='Transcriptome/post_perc/transcriptome_results_04.txt')
            # transcriptome_prefiltering.filter_proteins()
            # if not os.path.exists('Transcriptome/Percolator/for_predicting'):
            #     os.system('mkdir Transcriptome/Percolator/for_predicting')
            # transcriptome_prefiltering.filter_pin_files('Transcriptome/Percolator/for_predicting')
            # transcriptome_fishing = FeatureFishing(results='Transcriptome/post_perc/transcriptome_results_04.txt',
            #                                        pin_folder='Transcriptome/Percolator/for_predicting')
            # transcriptome_fishing.add_features()
            # transcriptome_fishing.save_table(
            #     'Transcriptome/Percolator/for_predicting/transcriptome_results_with_features.txt')

            transcriptome_sf = SpectralForest(model_pickle=f'{self.pypath}/src/model/model.pkl',
                                              results='Transcriptome/Percolator/for_predicting/transcriptome_results_with_features.txt')
            transcriptome_sf.predict()
            if not os.path.exists('Transcriptome/predicted'):
                os.system('mkdir Transcriptome/predicted')
            transcriptome_sf.save('Transcriptome/predicted/genome_predicted.txt')
