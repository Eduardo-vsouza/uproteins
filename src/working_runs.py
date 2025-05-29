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


import pandas as pd
import matplotlib.pyplot as plt
import random
import numpy as np
import seaborn as sns
from tqdm import tqdm


class MinimalRuns(object):
    def __init__(self, file, sample_number):
        self.file = file
        self.n = sample_number

    def get_df (self):
        df = pd.read_csv(self.file, sep='\t')
        return df

    def accession(self):
        df = self.get_df()
        entries = df["accession"].tolist()
        return entries

    def run_names(self):
        df = self.get_df()
        run = df["spectrumFile"].tolist()
        return run

    def sample_runs(self, n):
        runs = self.run_names()
        indexes = []
        for i in range(len(runs)):
            indexes.append(i)
        sampled = random.sample(indexes, n)
        return sampled

    def run_numbers(self):
        entries = self.accession()
        run = self.sample_runs(self.n)
        runs = self.run_names()

        run_dic = {

        }
        for i in range(len(entries)):
            if entries[i] not in run_dic:
                run_dic[entries[i]] = ""
                for j in range(len(run)):
                    if entries[run[j]] == entries[i]:
                        if runs[run[j]] not in run_dic[entries[i]].split(","):
                            run_dic[entries[i]] += "%s," % entries[run[j]]

        check_entries = []
        new_entries = []
        number_of_runs = []
        for i in range(len(entries)):
            if entries[i] not in check_entries:
                check_entries.append(entries[i])
                listtt = run_dic[entries[i]].split(",")
                # print(listtt)
                length = len(run_dic[entries[i]].split(","))
                # print(length)
                if length >= 2:
                    number_of_runs.append(length-1)
                elif length < 2:
                    number_of_runs.append(0)
                new_entries.append(entries[i])
        new_df = pd.DataFrame()
        new_df.insert(0, "Entry", new_entries)
        new_df.insert(1, "Number of Runs", number_of_runs)

        return new_df
        # new_df.to_csv("runs_per_orf7.xls", sep='\t', index=False)

    def diminishing_returns(self, n_runs, sampling):
        df = self.run_numbers()
        # df = run_numbers(file)
        entries = df["Entry"].tolist()
        runs = df["Number of Runs"].tolist()
        entries_1_run = []

        for i in range(n_runs):
            entries_1_run.append(0)

        for i in range(len(entries)):
            for j in range(n_runs):
                if runs[i] >= j:
                    entries_1_run[j] += 1

        # print(entries_1_run)
        new_df = pd.DataFrame()
        new_df.insert(0, "Minimum Number of Runs", range(n_runs))
        new_df.insert(1, "ORFs identified", entries_1_run)
        new_df.insert(2, "Sample", sampling)
        return entries_1_run


class OrganizePlot(object):
    def __init__(self, file):
        self.file = file

    def repeat(self):
        full_orfs = []
        full_runs = []
        for i in tqdm(range(100, 700)):
            if i % 100 == 0:
                for k in range(100):
                    genome = MinimalRuns(self.file, i)
                    dna_runs = genome.diminishing_returns(20, i)
                    number_of_runs = []
                    for j in range(len(dna_runs)):
                        number_of_runs.append(i)
                    full_orfs.append(dna_runs)
                    full_runs.append(number_of_runs)
        return full_orfs, full_runs

    def add_minimal_runs(self):
        orfs, runs = self.repeat()
        min_runs = []
        final_orfs = []
        final_runs = []
        for i in range(len(orfs)):
            for j in range(len(orfs[i])):
                final_orfs.append(orfs[i][j])
                final_runs.append(runs[i][j])
                min_runs.append(j+1)
        return min_runs, final_orfs, final_runs

    def create_df(self):
        min_runs, orfs, runs = self.add_minimal_runs()
        df = pd.DataFrame()
        df.insert(0, "ORFs", orfs)
        df.insert(1, "Runs", runs)
        df.insert(2, "Minimum", min_runs)
        df.to_csv("teste_df.xls", sep="\t", index=False)
        return df

    def plot_data(self):
        df = self.create_df()
        sns.catplot(data=df, x="Minimum", y="ORFs", hue="Runs", kind="bar")
        plt.title("ORFs that appeared in a minimum number of MS Runs")
        plt.ylabel("ORFs identified")
        plt.xlabel("Minimum Runs")
        plt.ylim(0, 60)
        plt.show()


# genome = OrganizePlot("/home/eduardo/Documents/Progenitus/RunTrial/pczao/Genome/Results/genome_unique_results_summarized.xls")
# genome.plot_data()
