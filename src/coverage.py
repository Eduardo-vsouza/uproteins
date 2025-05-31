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
import pandas as pd
from matplotlib_venn import venn2, venn2_circles


class Results(object):
    def __init__(self, folder, filetype, orf_db):
        self.folder = folder
        self.filetype = filetype
        self.orf_db = orf_db

    def coverage(self):
        df = pd.read_csv("%s/Results/Results_with_utps.txt" % self.folder, sep='\t')

        # proteins info
        proteins = df["ORF Sequence"].tolist()
        peptides = df["pepSeq"].tolist()
        # proteins = ["KCSTMAGEVIMCSWTEDCR", "VGDGPVGHDHPQHDRGTRGRAQSR"]
        # peptides = ["AGEV", "EDCR", "KCS", "EVIM", "GHDHP"]

        coverage_pct = []

        for i in range(len(proteins)):
            prot = proteins[i]
            coverage = 0
            nums = []
            print(len(prot))
            for j in range(len(peptides)):
                pep = peptides[j]
                if pep in prot:
                    print(pep)
                    start = prot.find(pep)
                    size = len(pep)
                    end = size + start
                    for k in range(start, end):
                        if k not in nums:
                            nums.append(k)
            for m in range(len(prot)):
                if m in nums:
                    coverage += 1
            coverage = (coverage/(len(prot)))*100
            coverage_pct.append(coverage)
        df.insert(4, "Coverage", coverage_pct, allow_duplicates=True)
        df.to_csv("%s/Results/Results_with_coverage" % self.folder, sep='\t', index=False)

    def create_fasta(self):
        df = pd.read_csv("%s/Results/Results_with_coverage" % self.folder, sep='\t')
        ids = df["accession"].tolist()
        seqs = df["ORF Sequence"].tolist()
        entries = []
        with open("%s/Results/%s_ORFs.fasta" %(self.folder, self.filetype), 'w') as fa:
            for i in range(len(ids)):
                entries.append(">" + ids[i] + "\n" + seqs[i] + "\n")
            fa.writelines(entries)


def create_venn():
    # Transcriptome data
    rna_df = pd.read_csv("Transcriptome/Results/results_with_seqs", sep='\t')
    rna_seqs = rna_df["ORF Sequence"]

    # Genome data
    dna_df = pd.read_csv("Genome/Results/results_with_seqs", sep='\t')
    dna_seqs = dna_df["ORF Sequence"]

    # Building venn diagram
    venn2([set(rna_seqs), set(dna_seqs)])
    plt.title('Shared ORFs')
    matplotlib.pyplot.savefig('venn.png')
