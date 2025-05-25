import os
import sys

import pandas as pd


class ResultsSummarizer:
    def __init__(self, transcriptome):
        self.transcriptomeResults = f'Transcriptome/Results'
        self.genomeResults = f'Genome/Results'
        self.transcriptome = transcriptome

        self.outfile = {'ORFs': [], 'Unique microprotein sequences': [], 'RBS': [], 'Subset': []}

    def get_results(self):
        subsets = {
            'Genome LC': f'{self.genomeResults}/genome_pre_validation_results.txt',
            'Genome HC': f'{self.genomeResults}/genome_post_validation_results.txt'
        }
        if self.transcriptome == 'YES':
            subsets['Transcriptome LC'] = f'{self.transcriptomeResults}/transcriptome_pre_validation_results.txt'
            subsets['Transcriptome HC'] = f'{self.transcriptomeResults}/transcriptome_post_validation_results.txt'

        for subset in subsets:
            counter = self.count_orfs(file=subsets[subset])
            self.outfile['ORFs'].append(len(counter['orfs']))
            self.outfile['Unique microprotein sequences'].append(len(counter['seqs']))
            self.outfile['RBS'].append(len(counter['shines']))
            self.outfile['Subset'].append(subset)

    @staticmethod
    def count_orfs(file):
        counter = {'orfs': [], 'seqs': [], 'shines': []}
        df = pd.read_csv(file, sep='\t')
        orfs = df["ORF name"].tolist()
        seqs = df["Protein sequence"].tolist()
        shines = df["Shine Dalgarno"].tolist()
        for i, orf in enumerate(orfs):
            if orf not in counter['orfs']:
                counter['orfs'].append(orf)
                if 'Present' in shines[i]:
                    counter['shines'].append(shines[i])
            if seqs[i] not in counter['seqs']:
                counter['seqs'].append(seqs[i])
        return counter

    def save(self):
        df = pd.DataFrame(data=self.outfile)
        df.to_csv('summarized_results.txt', sep='\t', index=False)