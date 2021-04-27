import pandas as pd


class StillCounting(object):
    def __init__(self, genome_df, transcriptome_df):
        """ utp_df must be a data frame containing only unique tryptic peptides. It should be generated with the
        get_utps() method from PercolatorUTP class. This class is only suitable to Percolator PSM output. """
        self.genomeDataFrame = pd.read_csv(genome_df, sep="\t")
        self.genomeDataFrame = self.genomeDataFrame[self.genomeDataFrame["proteinIds"].str.contains('ORF')]
        self.genomeDataFrame = self.genomeDataFrame[(self.genomeDataFrame["q-value"] < 0.01) & (self.genomeDataFrame["posterior_error_prob"] < 0.01)]
        self.genomeDataFrame = self.genomeDataFrame[self.genomeDataFrame["proteinIds"].str.contains('MSMEG') == False]

        self.transcriptomeDataFrame = pd.read_csv(transcriptome_df, sep="\t")
        self.transcriptomeDataFrame = self.transcriptomeDataFrame[self.transcriptomeDataFrame["proteinIds"].str.contains('ORF')]
        self.transcriptomeDataFrame = self.transcriptomeDataFrame[(self.transcriptomeDataFrame["q-value"] < 0.01) & (self.transcriptomeDataFrame["posterior_error_prob"] < 0.01)]
        self.transcriptomeDataFrame = self.transcriptomeDataFrame[self.transcriptomeDataFrame["proteinIds"].str.contains("MSMEG") == False]

        self.genomeProteins = self.genomeDataFrame["proteinIds"].tolist()
        self.transcriptomeProteins = self.transcriptomeDataFrame["proteinIds"].tolist()

    def __count_total(self, subset):
        orfs = []
        for proteins in subset:
            protein_set = proteins.split(",")
            for protein in protein_set:
                if protein not in orfs:
                    orfs.append(protein)
        return len(orfs)

    def __count_single(self, subset):
        """ Counts only a single ORF per peptide identification, excluding alternative start codons. """
        orfs = []
        for proteins in subset:
            protein = proteins.split(",")[0]
            if protein not in orfs:
                orfs.append(protein)
        return len(orfs)

    def count_orfs(self):
        genome_total = self.__count_total(self.genomeProteins)
        transcriptome_total = self.__count_total(self.transcriptomeProteins)

        genome_single = self.__count_single(self.genomeProteins)
        transcriptome_single = self.__count_single(self.transcriptomeProteins)

        print(f"__Total ORFs including alternative start codons__\nGenome: {genome_total}\nTranscriptome:"
              f" {transcriptome_total}\n__Total ORFs excluding alternative start codons__\nGenome: {genome_single}\n"
              f"Transcriptome: {transcriptome_single}")
        return self
