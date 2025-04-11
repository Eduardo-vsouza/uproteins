import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from venn import venn



class Metrics:
    def __init__(self, args):
        self.workingDir = os.getcwd()

        self.outdir = f'{self.workingDir}'
        print(self.outdir)
        self.args = args
        self.transcriptomeDir = f'{self.outdir}/Transcriptome'
        self.genomeDir = f'{self.outdir}/Genome'
        self.maxsize = args.maxSize
        self.metricsDir = f'{self.outdir}/metrics'
        self.__check_dirs([self.metricsDir])
        self.RBS = {}
        self.energies = {}
        self.data = {'microproteins': [], 'RBS': [], 'Energies': [], 'Database': [], 'Subset': []}

        self.genomeResultsDir = f'{self.genomeDir}/Results'
        self.transcriptomeResultsDir = f'{self.transcriptomeDir}/Results'
        self.genomePreValidation = f'{self.genomeDir}/Results/genome_pre_validation_results.txt'
        self.genomePostValidation = self.genomePreValidation.replace("pre_validation", "post_validation")
        self.transcriptomePreValidation = f'{self.transcriptomeDir}/Results/transcriptome_pre_validation_results.txt'
        self.transcriptomePostValidation = self.transcriptomePreValidation.replace("pre_validation", "post_validation")

        self.colors = ["#6A93E6", "lightgreen"]
        self.palette = {"Genome": "#6A93E6", "Transcriptome": "lightgreen"}

        self.fastaFiles = []


    @staticmethod
    def __check_dirs(folders):
        for folder in folders:
            if not os.path.exists(folder):
                os.mkdir(folder)

    def __check_redundant(self):
        self.df = self.df.drop_duplicates(subset=["microproteins"])

    def count_sequences(self):

        def count(df, subset, db):
            df = pd.read_csv(df, sep='\t')
            df = df[df["Protein sequence"].str.len() <= self.maxsize]
            df = df.drop_duplicates(subset=["Protein sequence"])
            # df = df.drop_duplicates(subset=["ORF name"])
            total_mips = df["Protein sequence"].tolist()
            rbs = df["Shine Dalgarno"].tolist()
            energies = df["Free Energy"].tolist()

            for i, mp in enumerate(total_mips):
                self.data['microproteins'].append(mp)
                self.data["RBS"].append(rbs[i])
                self.data["Energies"].append(energies[i])
                self.data["Subset"].append(subset)
                self.data["Database"].append(db)

        count(df=self.genomePreValidation, db="Genome", subset="Pre-validation")
        count(df=self.genomePostValidation, db="Genome", subset="Validated")
        count(df=self.transcriptomePreValidation, db="Transcriptome", subset="Pre-validation")
        count(df=self.transcriptomePostValidation, db="Transcriptome", subset="Validated")
        self.df = pd.DataFrame(data=self.data)
        # self.__check_redundant()


    def plot_aa_dist(self):
        plt.clf()
        aas = {}
        total = {}
        for i, mp in enumerate(self.data["microproteins"]):
            for aa in mp:
                db = self.data["Database"][i]
                if db not in aas:
                    aas[db] = {}
                subset = self.data["Subset"][i]
                if db not in total:
                    total[db] = 0
                if subset not in aas[db]:
                    aas[db][subset] = {}
                if aa not in aas[db][subset]:
                    aas[db][subset][aa] = 0
                aas[db][subset][aa] += 1
                total[db] += 1

        data = {'Amino acid': [], 'Subset': [], "Database": [], 'Count': []}
        plt.figure(figsize=(5, 6))

        for db in aas:
            for subset in aas[db]:
                for aa in aas[db][subset]:
                    data['Amino acid'].append(aa)
                    data['Count'].append(((int(aas[db][subset][aa]))/total[db])*100)
                    data['Subset'].append(subset)
                    data['Database'].append(db)
        ax = sns.barplot(data=data, x="Amino acid", hue="Database", ci=None, y="Count", palette=["#6A93E6", "lightgreen"],
                         edgecolor='black')
        plt.ylabel("Frequency (subset %)")
        sns.despine(top=True, right=True)

        plt.savefig(f'{self.metricsDir}/aa_distribution.pdf', dpi=600)

    def plot_orf_numbers(self):
        plt.clf()
        data = {'Microproteins': [], 'Subset': [], "Database": []}
        dbs = ["Genome", "Transcriptome"]
        subsets = ["Pre-validation", "Validated"]
        for db in dbs:
            for subset in subsets:
                df = self.df[self.df["Subset"] == subset]
                df = df[df["Database"] == db]
                data["Microproteins"].append(len(df["microproteins"].tolist()))
                data["Subset"].append(subset)
                data["Database"].append(db)

        df = pd.DataFrame(data=data)
        # plt.bar(x="Database", height="Microproteins", hue="Subset")
        plt.figure(figsize=(5, 6))

        sns.barplot(data=df, x='Subset', y='Microproteins', hue='Database', palette={"Genome": "#6A93E6", "Transcriptome": "lightgreen"},
                    edgecolor='black')

        for p in plt.gca().patches:
            if p.get_height() != 0:  # Skip annotation if height is 0

                plt.gca().annotate(str(int(p.get_height())), (p.get_x() + p.get_width() / 2., p.get_height()),
                                   ha='center', va='center', xytext=(0, 10), textcoords='offset points')

        plt.xlabel('Subset')
        plt.ylabel('Counts')
        plt.title('Barplot of Counts by Subset and Database')
        plt.legend(title='Database')
        sns.despine(top=True, right=True)
        plt.savefig(f'{self.metricsDir}/identified_microproteins.pdf', dpi=600)

    def plot_energies(self):
        plt.clf()
        # fig, ax = plt.subplots(figsize=(5, 6))
        plt.figure(figsize=(5, 6))

        # sns.set(rc={'figure.figsize': (11.7, 8.27)})
        ax = sns.boxplot(data=self.df, x="Subset", hue="Database", y="Energies", palette=self.palette)
        # plt.figure(figsize=(10,20))
        sns.despine(top=True, right=True)
        # sns.move_legend(ax, "upper right", bbox_to_anchor=(1, 1.2))
        plt.tight_layout()
        plt.legend([], [], frameon=False)

        plt.savefig(f'{self.metricsDir}/energies_distribution.pdf', dpi=600)


    def plot_rbs(self):
        plt.clf()
        # counts = self.df.groupby(["Subset", "Database", "RBS"]).size().reset_index(name="Count")
        rbs = set(self.df["RBS"].tolist())
        # Plot
        # plt.figure(figsize=(10, 6))
        # # sns.barplot(data=counts, x="Subset", y="Count", hue="Database", ci=None, edgecolor='black',
        # #             palette=self.palette)
        # counts = self.df.groupby(["Subset", "Database", "RBS"]).size().unstack(fill_value=0).reset_index()
        # counts_melted = pd.melt(counts, id_vars=["Subset", "Database"], value_vars=[rb for rb in rbs], var_name="String",
        #                         value_name="Count")
        # plt.figure(figsize=(10, 6))
        # sns.barplot(data=counts_melted, x="Subset", y="Count", hue="Database", ci=None, hue_order=["X", "Y"],
        #             palette="colorblind")
        # # plt.title("Occurrences of strings 'A', 'B', 'C' in column 'RBS' divided by 'Subset' with hue 'Database'")
        # plt.xlabel("Subset")
        # plt.ylabel("Count")
        # plt.show()
        # print(self.df)
        # counts = self.df.groupby(["Subset", "Database", "RBS"]).size().unstack(fill_value=0).reset_index()
        # df = pd.DataFrame(counts)
        # print(counts)
        rbs_done = []
        # data = {"Counts": [], "RBS": [], "Subset": []}
        data = {"RBS": [], "Database": [], "Subset": []}
        rbscount = {}
        dbs = self.df["Database"].tolist()
        subsets = self.df["Subset"].tolist()
        for i, rbs in enumerate(self.df["RBS"].tolist()):
            db = dbs[i]
            if rbs not in rbscount:
                rbscount[rbs] = {}
            if db not in rbscount[rbs]:
                rbscount[rbs][db] = {}
            subset = subsets[i]
            if subset not in rbscount[rbs][db]:
                rbscount[rbs][db][subset] = 0
            rbscount[rbs][db][subset] += 1
        for rbs in rbscount:
            for db in rbscount[rbs]:
                for sub in rbscount[rbs][db]:
                    data["RBS"].append(rbscount[rbs][db][sub])
                    data["Database"].append(db)
                    data["Subset"].append(rbs.replace("Leaderless", "Absent").replace("Present. Low to moderate binding", "Present. Low to\n moderate binding").replace("Strong binding", "\nStrong binding"))

        plt.figure(figsize=(5, 6))
        df = pd.DataFrame(data)

        sns.barplot(data=df, x="Subset", y="RBS", hue="Database", ci=None, palette=self.palette, estimator=sum,
                    edgecolor='black')
        # plt.title("Stacked Bar Plot of RBS values by Subset with Hue by Database")
        # plt.xlabel("Subset")
        # plt.ylabel("RBS Values")
        plt.legend(title="Database")
        sns.despine(top=True, right=True)
        # plt.tight_layout()
        # plt.xticks(rotation=45, ha='right')
        plt.tight_layout(h_pad=-5)

        plt.xlabel("Shine-Dalgarno")
        for p in plt.gca().patches:
            if p.get_height() != 0:  # Skip annotation if height is 0

                plt.gca().annotate(str(int(p.get_height())), (p.get_x() + p.get_width() / 2., p.get_height()),
                                   ha='center', va='center', xytext=(0, 10), textcoords='offset points')
        plt.savefig(f'{self.metricsDir}/shine-dalgarno.pdf', dpi=600)
        # total = tips.groupby('day')['total_bill'].sum().reset_index()

        # df.set_index('Subset').plot(kind='bar', stacked=True, color=['steelblue', 'red'])

        # ax = sns.histplot(data, hue='Database', weights='RBS', x="Subset", multiple='stack')

        # plt.title("Occurrences of strings 'A', 'B', 'C' in column 'RBS' divided by 'Subset' with hue 'Database'")
        # plt.xlabel("Subset")
        # plt.ylabel("Count")
        # plt.legend(title="Database")
        # Plot stacked bars
        # plt.figure(figsize=(10, 6))
        # plt.figure(figsize=(10, 6))
        # sns.barplot(data=counts, x="Subset", y=[rb for rb in rbs], hue="Database", ci=None)
        # plt.title("Occurrences of strings 'A', 'B', 'C' in column 'RBS' divided by 'Subset' with hue 'Database'")
        # plt.xlabel("Subset")
        # plt.ylabel("Count")

        # plt.title("smORF-containing transcripts with Ribosome binding sites (RBS)")
        # plt.xlabel("Subset")
        # plt.ylabel("RBS")

    def generate_fasta(self):
        files = [self.transcriptomePreValidation, self.transcriptomePostValidation, self.genomePreValidation,
                 self.genomePostValidation]
        folders = [self.transcriptomeResultsDir, self.transcriptomeResultsDir, self.genomeResultsDir,
                   self.genomeResultsDir]
        subsets = ['Pre-validation', 'Validated', 'Pre-validation', 'Validated']
        for i, file in enumerate(files):
            fasta = []
            df = pd.read_csv(file, sep='\t')
            df = df[df["Protein sequence"].str.len() <= self.maxsize]
            print("lens", len(df["Protein sequence"].tolist()))
            # df = self.df[self.df["Subset"] == subsets[i]]

            df = df.drop_duplicates(subset=["Protein sequence"])
            print("dupli", len(df["Protein sequence"].tolist()))
            entries, seqs = df["ORF name"].tolist(), df["Protein sequence"].tolist()
            fasta_file = f'{folders[i]}/{subsets[i]}.fasta'
            for entry, seq in zip(entries, seqs):
                fasta.append(f'>{entry}\n{seq}\n')
            with open(fasta_file, 'w') as handler:
                handler.writelines(fasta)
            self.fastaFiles.append(fasta_file)


    def plot_venn(self):
        validated = ''
        pre_validation = ''
        for fasta in self.fastaFiles:
            if 'Validated' in fasta:
                validated += f'{fasta},'
            else:
                pre_validation += f'{fasta},'
        val_venn = ResultsComparison(fasta_files=f'{validated}', groups='Transcriptome,Genome')
        val_venn.gather_microproteins()
        val_venn.plot_venn(palette=self.colors[::-1], output=f'{self.metricsDir}/venn_validated.png')

        pre_val_venn = ResultsComparison(fasta_files=pre_validation, groups='Transcriptome,Genome')
        pre_val_venn.gather_microproteins()
        pre_val_venn.plot_venn(palette=self.colors[::-1], output=f'{self.metricsDir}/venn_pre_validation.png')


    def plot_venn_2(self):
        print(self.df)
        # rna_pre = self.df[(self.df["Subset"== "Pre-validation"]) & (self.df["Database"] == "Transcriptome")]
        print(self.df[(self.df["Subset" == "Pre-validation"]) & (self.df["Database"] == "Transcriptome")])

        mps = {"Transcriptome": set(self.df[(self.df["Subset" == "Pre-validation"]) & (self.df["Database"] == "Transcriptome")]["Protein sequence"].tolist()),
               "Genome": set(self.df[(self.df["Subset" == "Pre-validation"]) & (self.df["Database"] == "Genome")]["Protein sequence"].tolist())}
        # venn(self.microproteins, cmap=['#1D7DBB', '#BB1D1D'], fontsize=14)
        venn(mps, fontsize=16, cmap=self.colors[::-1])
        # plt.legend([])
        plt.show()
        # plt.savefig(f'{self.metricsDir}/')

class ResultsComparison:
    def __init__(self, fasta_files, groups):
        self.folders = fasta_files.split(",")
        self.groups = groups.split(",")

        self.microproteins = {}

    def gather_microproteins(self):
        for fasta, group in zip(self.folders, self.groups):
            if group not in self.microproteins:
                self.microproteins[group] = set()
            records = SeqIO.parse(fasta, 'fasta')
            for record in records:
                # self.microproteins[group].add(str(record.seq))
                self.microproteins[group].add(str(record.seq))

    # def gather_mps_from_df(self):

