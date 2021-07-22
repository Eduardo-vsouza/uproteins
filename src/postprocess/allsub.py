import os
import pandas as pd
from Bio import SeqIO


class AllSub(object):
    def __init__(self, folder, filetype):
        self.folder = folder
        self.filetype = filetype
        self.fileToFix = f'{self.folder}/Percolator/{self.folder}_pin.txt'
        self.dfToFix = pd.read_csv(self.fileToFix, sep='\t', usecols=[i for i in range(35)])
        self.catProteins = []
        # self._fix_columns()
        # self._save_fixed()

        # self.pinFile = pd.read_csv(f'{self.folder}/Percolator/{self.folder}_pin.txt', sep='\t')

    def remove_annotated(self):
        removed = []
        with open(self.fileToFix, 'r') as tofix, open(f'{self.folder}/Percolator/{self.folder}_pin_removed.txt', 'w') as out:
            lines = tofix.readlines()
            for line in lines:
                if 'ANNO' not in line:
                    removed.append(line)
            out.writelines(removed)
        os.system(f'mv {self.fileToFix} {self.fileToFix}_all.txt')
        os.system(f'mv {self.folder}/Percolator/{self.folder}_pin_removed.txt {self.fileToFix}')

    def _fix_columns(self):
        with open(self.fileToFix) as unfixed:
            lines = unfixed.readlines()
            i = 0
            for line in lines:
                line = line.rstrip()
                cols = line.split("\t")
                new_proteins = ""
                # if i > 2:
                if i > 1:
                    if len(cols) >= 35:
                        for col in range(36, len(cols)):
                            if col == len(cols):
                                new_proteins += f'{cols[col]},\n'
                            else:
                                new_proteins += f'{cols[col]},'
                    else:
                        new_proteins += f'{cols}\n'
                    self.catProteins.append(new_proteins)
                elif i == 1:
                    self.catProteins.append('no')
                i += 1

    def _save_fixed(self):
        cmd_rename_old = f'mv {self.folder}/Percolator/{self.folder}_pin.txt {self.folder}/Percolator/{self.folder}_pin_unfixed.txt'
        os.system(cmd_rename_old)
        self.dfToFix.insert(35, "Proteins", self.catProteins)
        self.dfToFix.to_csv(self.fileToFix, sep='\t', index=False)

    def remove_irrelevant(self):
        """ Removes from the pin files any peptide that matches an annotated protein, leaving only those relevant
        to the hypothesis of the pipeline. """
        self.pinFile = self.pinFile[self.pinFile["Proteins"].str.contains('ANNO') == False]
        os.system(f'mv {self.folder}/Percolator/{self.folder}_pin.txt {self.folder}/Percolator/{self.folder}_all_pin.txt')

    def modify_decoy(self):
        """ Removes from the decoy the sequences of annotated proteins """
        records = SeqIO.parse(f'{self.folder}/Percolator/{self.folder}_decoy.fasta', 'fasta')
        fasta = []
        for record in records:
            if 'ANNO' not in str(record.description):
                fasta.append(f'>{str(record.description)}\n{str(record.seq)}\n')
        self._rename_decoy()
        with open(f'{self.folder}/Percolator/{self.folder}_decoy.fasta', 'w') as out:
            out.writelines(fasta)

    def _rename_decoy(self):
        cmd_mv = f'mv {self.folder}/Percolator/{self.folder}_decoy.fasta {self.folder}/Percolator/{self.folder}_decoy_all.fasta'
        os.system(cmd_mv)

    def save(self):
        self.pinFile.to_csv(f'{self.folder}/Percolator/{self.folder}_pin.txt', sep='\t', index=False)
        to_write = []
        with open(self.fileToFix, 'r') as tofix, open(f'{self.fileToFix}_refixed', 'w') as out:
            lines = tofix.readlines()
            for line in lines:
                if 'SpecId' not in line and 'Default' not in line:
                    line = line.replace(",", "\t")
                to_write.append(line)
                print(line)
            out.writelines(to_write)
        cmd_mv = f'mv {self.fileToFix} {self.fileToFix}_old.txt'
        cmd_mv2 = f'mv {self.fileToFix}_refixed {self.fileToFix}'
        os.system(cmd_mv)
        os.system(cmd_mv2)

