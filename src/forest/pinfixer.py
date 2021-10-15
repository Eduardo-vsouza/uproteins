import pandas as pd
import os


class ProteinFixer(object):
    def __init__(self, pin_folder):
        self.folder = pin_folder
        self.files = [f'{pin_folder}/{file}' for file in os.listdir(self.folder)]

        # self.done = [f'{pin_folder}/fixed{file}' for file in os.listdir(f'{self.folder}/fixed')]

    def fix_files(self, outdir):
        for file in self.files:
            filename = file.split("/")

            # if 'pin' in file and f'{outdir}/fixed_{filename[len(filename)-1]}' not in self.done:
            if 'pin' in file or 'pep_results' in file:

                print(file)
                with open(file, 'r') as unfixed, open(f'{outdir}/for_predicting/fixed_{filename[len(filename)-1]}', 'w') as out:
                    lines = unfixed.readlines()
                    i = 0
                    new_lines = []
                    for line in lines:
                        if 'Default' not in line:

                            cols = line.split("\t")
                            if i == 0:
                                proteins_col = len(cols) - 1
                            i += 1
                            if len(cols) > proteins_col + 1:
                                new_proteins = ','.join(cols[proteins_col:])
                            else:
                                new_proteins = ','.join(cols[proteins_col:])
                            new_line = '\t'.join(cols[:proteins_col])+'\t'
                            new_line += new_proteins
                            new_lines.append(new_line)
                            # new_lines.append(new_proteins)
                    out.writelines(new_lines)


if __name__ == '__main__':
    folder = 'uproteins_results_1608/tofix'
    data = ProteinFixer(pin_folder=folder)
    data.fix_files(outdir='uproteins_results_1608')