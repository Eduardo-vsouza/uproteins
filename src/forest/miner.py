import pandas as pd
import os


class SpecMiner(object):
    def __init__(self, results, pin_dir):
        self.df = pd.read_csv(results, sep='\t')
        self.pinFolder = pin_dir
        self.df = self.df.drop_duplicates()
        # self.fixedFileNames = self.__fix_filenames()
        self.fixedFileNames = self.df["SpecId"]
        self.scanNums = self.df["ScanNum"].tolist()
        self.extendedORFs = self.df["Extended ORF"].tolist()
        self.extendedSequences = self.df["Extended Sequence"].tolist()

    def __fix_filenames(self):
        """ legacy code """
        names = self.df["SpecFile"].tolist()
        fixed = [name[:-5] for name in names]
        self.df.insert(0, "FixedFiles", fixed)
        return fixed

    def mine_spectra(self, output):
        files = os.listdir(self.pinFolder)
        i = 0
        extended_orfs = []
        extended_seqs = []
        for scan, file in zip(self.scanNums, self.fixedFileNames):
            for file in files:
                full_df = pd.read_csv(f'{self.pinFolder}/{file}', sep='\t')
                if i == 0:
                    df = pd.DataFrame(columns=full_df.columns)
                    # df.columns = full_df.columns
                intersection = full_df[(full_df["ScanNr"] == scan) & (full_df["SpecId"].str.contains(file))]
                df = df.append(intersection)
                scans = intersection["ScanNr"].tolist()
                for j in range(len(scans)):
                    extended_orfs.append(self.extendedORFs[i])
                    extended_seqs.append(self.extendedSequences[i])
            print(i)
            i += 1

            if i == 5:
                break

        df.insert(3, "Extended ORF", extended_orfs)
        df.insert(4, "Extended Sequence", extended_seqs)
        df.to_csv(output, sep='\t', index=False)


class MinePredicted(object):
    def __init__(self, predicted, results):
        self.predicted = pd.read_csv(predicted, sep='\t')
        self.df = pd.read_csv(results, sep='\t')

        self.specFiles = self.df["SpecFile"].tolist()
        self.scans = self.df["ScanNum"].tolist()

        self.highConfidenceScans = self.predicted["ScanNr"].tolist()
        self.highConfidenceFiles = self.predicted["SpecId"].tolist()

    def filter(self, output):
        high_scans = []
        high_files = []
        for file, scan in zip(self.specFiles, self.scans):
            for i in range(len(self.highConfidenceFiles)):
                if file[:-5] in self.highConfidenceFiles[i] and self.highConfidenceScans[i] == scan:
                    high_scans.append(scan)
                    high_files.append(file)
        df = pd.DataFrame(columns=self.df.columns)
        for file, scan in zip(high_files, high_scans):
            ndf = self.df[(self.df["SpecFile"] == file) & (self.df["ScanNum"] == scan)]
            df = df.append(ndf)
        df = df.drop_duplicates()
        df.to_csv(output, sep='\t', index=False)

    def filter_v2(self):
        hcs = self.__get_hc()



    def __get_hc(self):
        hcs = {}  # spectra predicted to have high confidence by the ML model
        for file, scan in zip(self.highConfidenceFiles, self.highConfidenceScans):
            if file in hcs:
                hcs[file] += f',{scan}'
            else:
                hcs[file] = f'{scan}'
        return hcs

# miner = SpecMiner(results='Genome/genome_results_04.txt')
# miner.mine_spectra('teste_miner.txt')

# miner = MinePredicted(predicted='../for_predicting/genome_predicted.txt', results='Genome/genome_results_04.txt')
# miner.filter('../predicted/genome_04_predicted.txt')


class SpectrumMiner(object):
    def __init__(self, results, predicted):
        self.results = pd.read_csv(results, sep='\t')
        self.unfilteredScans = self.results["ScanNum"].tolist()
        self.__rename()

        self.predicted = pd.read_csv(predicted, sep='\t')


        self.highConfidenceScans = self.predicted["ScanNr"].tolist()
        self.highConfidenceFiles = self.predicted["SpecId"].tolist()
        self.__rename_predicted()
        self.fixedPredictedFiles = self.predicted["RenamedFiles"].tolist()
        print(self.fixedPredictedFiles)
        self.HCs = self.__get_hc()
        # self.results = self.results[self.results["RenamedFile"].isin(self.HCs)]

    def __get_hc(self):
        hcs = {}  # spectra predicted to have high confidence by the ML model
        for file, scan in zip(self.fixedPredictedFiles, self.highConfidenceScans):
            if file in hcs:
                hcs[file].append(scan)
            else:
                hcs[file] = [scan]
        return hcs

    def __rename_predicted(self):
        renamed = ['_'.join(file.split("_")[:-6]) for file in self.highConfidenceFiles]
        self.predicted.insert(0, 'RenamedFiles', renamed)

    def __rename(self):
        spec_files = self.results["SpecFile"].tolist()
        renamed = []
        for file, scan in zip(spec_files, self.unfilteredScans):
            renamed.append(file[:-5])
        self.results.insert(0, "RenamedFiles", renamed)
        # self.results.insert(3, "FixedScans", fixed_scan)

    def filter_results(self, output):
        df = pd.DataFrame(columns=self.results.columns)
        i = 0
        for hc in self.HCs:
            i += 1
            df = df.append(self.results[(self.results["ScanNum"].astype('int64').isin(self.HCs[hc])) & (self.results["RenamedFiles"] == hc)])
            print(i, len(self.HCs))
        df = df.drop_duplicates()
        df.to_csv(f'{output}', sep='\t', index=False)


# data = SpectrumMiner(predicted='../for_predicting/predicted_pin_files.txt', results='Transcriptome/transcriptome_results_04.txt')
# data.filter_results('Transcriptome/post_perc/transcriptome_results_05.txt')

# data = SpectrumMiner(predicted='../for_predicting/predicted_pin_files_supervised.txt', results='../post_perc/genome_results_05.txt')
# rna_data = SpectrumMiner(predicted='../for_predicting/Transcriptome/predicted_pin_files_semisupervised.txt', results='Transcriptome/transcriptome_results_04.txt')
# data ='../for_predicting/genome_predicted.txt', results='Genome/genome_results_04.txt')
# data = SpectrumMiner(results='df2', predicted='df1')
# data.filter_results()
# rna_data.filter_results('Transcriptome/transcriptome_orfs_semisupervised_24_06')

# mtb_exosomes = SpectrumMiner(predicted='predicted/exosomes_mtb_predicted.txt', results='u')