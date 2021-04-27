import pandas as pd
import glob
import numpy as np
import os
import time
import main
from matplotlib import pyplot as plt
from Bio import SeqIO
from tqdm import tqdm


class Digestion(object):
    def __init__(self, database):
        self.database = database

    def digest_orfs(self, enzymes):
        pypath = main.pypath
        enzyme_list = enzymes.split(",")
        for enzyme in enzyme_list:
            cmd_digest = 'perl %s/protein2digest.pl -in=2 -ez=%s -mc=3 %s' % (pypath, enzyme, self.database)
            os.system(cmd_digest)
            cmd_rename_dir = 'mv %s ./%s_%s' % (self.database[:-6], self.database[:-6], enzyme)
            os.system(cmd_rename_dir)


class Digested(object):
    def __init__(self, folder, filetype, orf_db, enzyme, proteome):
        self.folder = folder
        self.filetype = filetype
        self.orf_db = orf_db
        self.enzyme = enzyme
        self.proteome = proteome

    def virtual_digestion(self, mix):
        # mix must be either 'no' or the path to the other enzyme digestion folder, i.e. a second self.folder
        if not os.path.exists("Virtual_digestion/"):
            cmd_dir = 'mkdir Virtual_digestion'
            os.system(cmd_dir)

        files = glob.glob("%s/*PROT*" % self.folder)

        if mix != "no":
            second_files = glob.glob("%s/*PROT*" % mix)
            for i in range(len(second_files)):
                files.append(second_files[i])

        the_Dic = {

        }
        print("\nDigesting your ORF database.\n")
        for i in range(len(files)):
            with open(files[i], 'r') as name:
                names = []
                lines = name.readlines()
                for line in range(len(lines)):
                    lines[line] = lines[line].rstrip()
                    lines[line] = lines[line][3:]
                names.append(lines[1])
                if names[0] not in the_Dic:
                    the_Dic[names[0]] = ""
            df = pd.read_csv(files[i], sep='\t', header=3, keep_default_na=False, na_values=['_'])
            pepseq = df["peptide sequence"].tolist()
            for j in range(len(pepseq)):
                if 5 <= len(pepseq[j]) <= 20:
                    the_Dic[names[0]] += "%s," % pepseq[j]
        return the_Dic

    def virtual_coverage(self, mix):
        """ Verifies the maximum obtainable coverage for each ORF when using a specific protease for protein digestion.
        It then writes the coverage into a tab-separated file for it to be plotted in the next step. """
        ids = []
        seqs = []
        parser = SeqIO.parse(self.orf_db, 'fasta')
        for entry in parser:
            ids.append(entry.id)
            seqs.append(entry.seq)

        the_Dic = self.virtual_digestion(mix)
        all_coverage = []
        final_df = pd.DataFrame()
        print("\nPredicting virtual coverage\n")
        for entry in tqdm(the_Dic):
            peptides = the_Dic[entry].split(",")[:-1]
            coverage = 0
            nums = []
            id_index = ids.index(entry)
            seq = seqs[id_index]
            for i in range(len(peptides)):
                pep = peptides[i]
                if pep in seq:
                    start = seq.find(pep)
                    size = len(pep)
                    end = size + start
                    for k in range(start, end):
                        if k not in nums:
                            nums.append(k)
            for m in range(len(seq)):
                if m in nums:
                    coverage += 1
            coverage = (coverage / (len(seq))) * 100
            all_coverage.append(coverage)
        final_df.insert(0, "Coverages", all_coverage)
        if mix != "no":
            final_df.to_csv("Virtual_digestion/%s_coverage_mixed_tryp_lys-c.txt" % self.enzyme, sep='\t', index=False)
        elif mix == "no":
            final_df.to_csv("Virtual_digestion/%s_coverage.txt" % self.enzyme, sep='\t', index=False)

    def virtual_ups(self):
        """ Access the number of possible unique peptides from an in silico protein digestion using a specific
         enzyme. The criteria for classifying a peptide as unique follows the same rules as those for identifying unique
         peptides from MS experiments. """

        the_Dic = self.virtual_digestion()
        # print(the_Dic)
        peptides = []
        entries = []
        for entry in the_Dic:
            peps = set(the_Dic[entry].split(",")[:-1])
            for i in peps:
                peptides.append(i)
                entries.append(entry)
        df = pd.DataFrame()
        df.insert(0, "Entry", entries)
        df.insert(1, "Peptide", peptides)

        unique = []

        # database info
        orf_ids = []
        orf_seqs = []
        orf_locus_start = []
        orf_locus_end = []

        # RefSeq info
        anno_ids = []
        anno_seqs = []

        # Add ORF ids and seqs to the database info
        with open("%s_database_no_anno.fasta" % self.filetype, 'r') as db, open(self.proteome, 'r') as anno:
            records = SeqIO.parse(db, 'fasta')
            for record in records:
                orf_ids.append(str(record.id))
                orf_seqs.append(str(record.seq))
            ref_records = SeqIO.parse(anno, 'fasta')
            for record in ref_records:
                anno_ids.append(str(record.id))
                anno_seqs.append(str(record.seq))

        # Add info about the ORF loci in the genome to database info
        for i in range(len(orf_ids)):
            icon = orf_ids[i].find("_-_")
            beginning = orf_ids[i].find("[")
            ending = orf_ids[i].find("]")
            locus_start = orf_ids[i][beginning + 1:icon]
            locus_end = orf_ids[i][icon + 3:ending]
            orf_locus_start.append(locus_start)
            orf_locus_end.append(locus_end)


        checked_peptides = {

        }
        # Check if the peptide identified by MS is unique or not
        timer = 0
        print("\nCounting Virtual Unique Peptides\n")
        for pep in tqdm(range(len(peptides))):
            if peptides[pep] not in checked_peptides:
                if peptides[pep] != "":
                    # print(peptides[pep])
                    # timer += 1
                    # print((timer / len(peptides)) * 100)
                    appearances = 0
                    ref_appear = 0
                    entries = []
                    entries_start = []
                    entries_end = []
                    unique_pep = True
                    # print(orf_seqs)
                    # is_pep_in_orf = [orf for orf in orf_seqs if peptides[pep] in orf]
                    # ind = orf_seqs.index(is_pep_in_orf[0])
                    # print(is_pep_in_orf)
                    for i in range(len(orf_seqs)):
                        if peptides[pep] in orf_seqs[i]:
                            entries.append(orf_ids[i])
                            entries_start.append(orf_locus_start[i])
                            entries_end.append(orf_locus_end[i])
                            appearances += 1
                    if [orf for orf in anno_seqs if peptides[pep] in orf]:
                        # for i in range(len(anno_seqs)):
                        #     if peptides[pep] in anno_seqs[i]:
                        ref_appear += 1
                    if appearances > 1:
                        for i in range(len(entries)):
                            for j in range(len(entries_start)):
                                entries_length = len(entries)
                                if int(entries_start[entries_length - 1]) not in range(int(entries_start[j]),
                                                                                       int(entries_end[j])):
                                    unique_pep = False
                    elif appearances == 0 and ref_appear > 1:
                        unique_pep = False
                    elif appearances >= 1 and ref_appear >= 1:
                        unique_pep = False
                    unique.append(unique_pep)
                    checked_peptides[peptides[pep]] = unique_pep
            else:
                unique.append(checked_peptides[peptides[pep]])
        df.insert(2, "Unique Peptide", unique, allow_duplicates=True)
        return df
        # df.to_csv("Virtual_digestion/Unique_peptides_from_%s_%s_digestion.txt" % (self.enzyme, self.filetype), sep='\t',
        #           index=False)

    def ups_per_orf(self):
        df = self.virtual_ups()
        orf_dic = {}
        # df = pd.read_csv("Virtual_digestion/Unique_peptides_from_%s_%s_digestion.txt" % (self.enzyme, self.filetype),
        #                  sep='\t')
        df = df.loc[df['Unique Peptide'] == True]
        entries = df["Entry"].tolist()
        peptides = df["Peptide"].tolist()
        for i in range(len(entries)):
            if entries[i] not in orf_dic:
                orf_dic[entries[i]] = ""
            orf_dic[entries[i]] += "%s," % peptides[i]
        new_entries = []
        number_ups = []
        for entry in orf_dic:
            entry_peptides = orf_dic[entry].split(",")[:-1]
            number = 0
            new_entries.append(entry)
            for i in range(len(entry_peptides)):
                number += 1
            number_ups.append(number)
        df_with_ups = pd.DataFrame()
        df_with_ups.insert(0, "ORF", new_entries)
        df_with_ups.insert(1, "Number of obtainable unique peptides", number_ups)
        df_with_ups.to_csv("Virtual_digestion/%s_%s_ups" % (self.enzyme, self.filetype), sep='\t', index=False)


class PlotData(object):
    def __init__(self, enzymes):
        """List of enzymes must be separated by a comma (,) and cannot include any spaces. """
        self.enzymes = enzymes

    def plot_coverage_data(self):
        enzyme_list = self.enzymes.split(",")
        trypsin_df = pd.read_csv("Virtual_digestion/Trypsin_coverage.txt", sep="\t")
        lys_df = pd.read_csv("Virtual_digestion/Lys_C_coverage.txt", sep="\t")
        argc_df = pd.read_csv("Virtual_digestion/Arg_C_coverage.txt", sep="\t")
        gluc_df = pd.read_csv("Virtual_digestion/Glu_C_coverage.txt", sep="\t")
        mix_tl_df = pd.read_csv("Virtual_digestion/Trypsin_coverage_mixed_tryp_lys-c.txt", sep="\t")
        tryp_cov = trypsin_df["Coverages"].tolist()
        lys_cov = lys_df["Coverages"].tolist()
        argc_cov = argc_df["Coverages"].tolist()
        gluc_cov = gluc_df["Coverages"].tolist()
        tl_cov = mix_tl_df["Coverages"].tolist()
        # try_lysc = tryp_cov + lys_cov
        # try_argc = tryp_cov + argc_cov
        labels = ["Trypsin", "Lys-C", "Arg-C", "Glu-C", "Trypsin+Lysc"]
        tryp_mean = np.mean(tryp_cov)
        lys_mean = np.mean(lys_cov)
        means = [tryp_mean, lys_mean, np.mean(argc_cov), np.mean(gluc_cov), np.mean(tl_cov)]
        print(means)
        plt.title("Coverage after digestion by different enzymes")
        plt.ylabel("Coverage")
        plt.xlabel("Enzyme")
        plt.bar(labels, means, edgecolor="black", color="lightblue")
        plt.show()

    def mix_up_data(self):
        enzyme_list = self.enzymes.split(",")
        for i in range(len(enzyme_list)):
            if enzyme_list[i] == 0:
                enzyme_list[1] = "Trypsin"
            elif enzyme_list[i] == 1:
                enzyme_list[i] = "Lys_C"
            elif enzyme_list[i] == 2:
                enzyme_list[i] = "Arg_C"
            elif enzyme_list[i] == 3:
                enzyme_list[i] = "Glu_C"
        up_tables = []
        for i in range(len(enzyme_list)):
            up_tables.append("Virtual_digestion/%s_genome_ups" % enzyme_list[i])

        # trypsin column info
        tryp_df = pd.read_csv("Virtual_digestion/Trypsin_genome_ups", sep='\t')
        tryp_ups = tryp_df["Number of obtainable unique peptides"].tolist()
        tryp_orfs = tryp_df["ORF"].tolist()

        # Lys-C column info
        lysc_df = pd.read_csv("Virtual_digestion/Lys_C_genome_ups", sep='\t')
        lysc_ups = lysc_df["Number of obtainable unique peptides"].tolist()
        lysc_orfs = lysc_df["ORF"].tolist()

        # Arg-C column info
        argc_df = pd.read_csv("Virtual_digestion/Arg_C_genome_ups", sep="\t")
        argc_ups = argc_df["Number of obtainable unique peptides"].tolist()
        argc_orfs = argc_df["ORF"].tolist()

        # Glu-C column info
        gluc_df = pd.read_csv("Virtual_digestion/Glu_C_genome_ups", sep="\t")
        gluc_ups = gluc_df["Number of obtainable unique peptides"].tolist()
        gluc_orfs = gluc_df["ORF"].tolist()

        # Trypsin + Lys-C info
        tl_up_per_orf = []

        for i in range(len(tryp_orfs)):
            if tryp_orfs[i] in lysc_orfs:
                finder = lysc_orfs.index(tryp_orfs[i])
                ups = lysc_ups[finder] + tryp_ups[i]
                tl_up_per_orf.append(ups)

        # Trypsin + Arg-C info
        ta_up_per_orf = []

        for i in range(len(tryp_orfs)):
            if tryp_orfs[i] in argc_orfs:
                finder = argc_orfs.index(tryp_orfs[i])
                ups = argc_ups[finder] + tryp_ups[i]
                ta_up_per_orf.append(ups)

        # Trypsin + Glu-C info
        tg_up_per_orf = []

        for i in range(len(tryp_orfs)):
            if tryp_orfs[i] in gluc_orfs:
                finder = gluc_orfs.index(tryp_orfs[i])
                ups = gluc_ups[finder] + tryp_ups[i]
                tg_up_per_orf.append(ups)

        return tl_up_per_orf, ta_up_per_orf, tg_up_per_orf

    def organize_up_data(self):
        enzyme_list = self.enzymes.split(",")
        for i in range(len(enzyme_list)):
            if enzyme_list[i] == 0:
                enzyme_list[1] = "Trypsin"
            elif enzyme_list[i] == 1:
                enzyme_list[i] = "Lys_C"
            elif enzyme_list[i] == 2:
                enzyme_list[i] = "Arg_C"
            elif enzyme_list[i] == 3:
                enzyme_list[i] = "Glu_C"
        up_tables = []
        for i in range(len(enzyme_list)):
            up_tables.append("Virtual_digestion/%s_genome_ups" % enzyme_list[i])
        tryp_df = pd.read_csv("Virtual_digestion/Trypsin_genome_ups", sep='\t')
        tryp_ups = tryp_df["Number of obtainable unique peptides"].tolist()
        tryp_up_number = 0
        tryp_orf_number_1 = 0
        tryp_orf_number_2 = 0
        for i in range(len(tryp_ups)):
            tryp_up_number += int(tryp_ups[i])
            if int(tryp_ups[i]) >= 1:
                tryp_orf_number_1 += 1
            if int(tryp_ups[i]) >= 2:
                tryp_orf_number_2 += 1
        lysc_df = pd.read_csv("Virtual_digestion/Lys_C_genome_ups", sep='\t')
        lysc_ups = lysc_df["Number of obtainable unique peptides"].tolist()
        lysc_up_number = 0
        lysc_orf_number_1 = 0
        lysc_orf_number_2 = 0
        for i in range(len(lysc_ups)):
            lysc_up_number += int(lysc_ups[i])
            if int(lysc_ups[i]) >= 1:
                lysc_orf_number_1 += 1
            if int(lysc_ups[i]) >= 2:
                lysc_orf_number_2 += 1
        argc_df = pd.read_csv("Virtual_digestion/Arg_C_genome_ups", sep="\t")
        argc_ups = argc_df["Number of obtainable unique peptides"].tolist()
        argc_up_number = 0
        argc_orf_number_1 = 0
        argc_orf_number_2 = 0
        for i in range(len(argc_ups)):
            argc_up_number += int(argc_ups[i])
            if int(argc_ups[i]) >= 1:
                argc_orf_number_1 += 1
            if int(argc_ups[i]) >= 2:
                argc_orf_number_2 += 1
        gluc_df = pd.read_csv("Virtual_digestion/Glu_C_genome_ups", sep="\t")
        gluc_ups = gluc_df["Number of obtainable unique peptides"].tolist()
        gluc_up_number = 0
        gluc_orf_number_1 = 0
        gluc_orf_number_2 = 0
        for i in range(len(gluc_ups)):
            gluc_up_number += int(gluc_ups[i])
            if int(gluc_ups[i]) >= 1:
                gluc_orf_number_1 += 1
            if int(gluc_ups[i]) >= 2:
                gluc_orf_number_2 += 1

        # height of bars
        trypsin = [tryp_orf_number_1, tryp_orf_number_2, tryp_up_number/10]
        lysc = [lysc_orf_number_1, lysc_orf_number_2, lysc_up_number/10]
        argc = [argc_orf_number_1, argc_orf_number_2, argc_up_number/10]
        gluc = [gluc_orf_number_1, gluc_orf_number_2, gluc_up_number/10]

        return trypsin, lysc, argc, gluc

    def plot_up_orf(self):
        trypsin, lysc, argc, gluc = self.organize_up_data()

        barWidth = 0.21

        # bar position on X axis
        r1 = np.arange(len(trypsin))
        r2 = [x + barWidth for x in r1]
        r3 = [x + barWidth for x in r2]
        r4 = [x + barWidth for x in r3]

        # plot
        plt.bar(r1, trypsin, color="lightblue", width=barWidth, edgecolor="black", label="Trypsin")
        plt.bar(r2, lysc, color="darkred", width=barWidth, edgecolor="black", label="Lys_C")
        plt.bar(r3, argc, color="goldenrod", width=barWidth, edgecolor="black", label="Arg_C")
        plt.bar(r4, gluc, color="forestgreen", width=barWidth, edgecolor="black", label="Glu_C")

        # xticks in the middle of group bars
        plt.xlabel("Unique Peptides", fontweight="bold")
        plt.xticks([r+barWidth for r in range(len(trypsin))], ["ORFs with >= 1 UP", "ORFs with >= 2 UP",
                                                               "Number of UPs per ORF / 10"])

        plt.title("Variation in the number of Unique Peptides when using different proteases")
        plt.legend()
        plt.savefig("Proteases_ups.png")

    def organize_mixed_proteases(self):
        tl_up_per_orf, ta_up_per_orf, tg_up_per_orf = self.mix_up_data()

        # Trypsin + Lys-C
        tl_up_number = 0
        tl_orf_number_1 = 0
        tl_orf_number_2 = 0
        for i in range(len(tl_up_per_orf)):
            tl_up_number += int(tl_up_per_orf[i])
            if int(tl_up_per_orf[i]) >= 1:
                tl_orf_number_1 += 1
            if int(tl_up_per_orf[i]) >= 2:
                tl_orf_number_2 += 1

        # Trypsin + Arg-C
        ta_up_number = 0
        ta_orf_number_1 = 0
        ta_orf_number_2 = 0
        for i in range(len(ta_up_per_orf)):
            ta_up_number += int(ta_up_per_orf[i])
            if int(ta_up_per_orf[i]) >= 1:
                ta_orf_number_1 += 1
            if int(tg_up_per_orf[i]) >= 2:
                ta_orf_number_2 += 1

        # Trypsin + Glu-C
        tg_up_number = 0
        tg_orf_number_1 = 0
        tg_orf_number_2 = 0
        for i in range(len(tg_up_per_orf)):
            tg_up_number += int(tg_up_per_orf[i])
            if int(tg_up_per_orf[i]) >= 1:
                tg_orf_number_1 += 1
            if int(tg_up_per_orf[i]) >= 2:
                tg_orf_number_2 += 1

        tl_info = [tl_orf_number_1, tl_orf_number_2, tl_up_number/10]
        ta_info = [ta_orf_number_1, ta_orf_number_2, ta_up_number/10]
        tg_info = [tg_orf_number_1, tg_orf_number_2, tg_up_number/10]

        return tl_info, ta_info, tg_info

    def plot_mixed_data(self):
        tl, ta, tg = self.organize_mixed_proteases()

        barWidth = 0.21

        # bar position on X axis
        r1 = np.arange(len(tl))
        r2 = [x + barWidth for x in r1]
        r3 = [x + barWidth for x in r2]

        # plot
        plt.bar(r1, tl, color="lightblue", width=barWidth, edgecolor="black", label="Trypsin+Lys-C")
        plt.bar(r2, ta, color="darkred", width=barWidth, edgecolor="black", label="Trypsin+Arg-C")
        plt.bar(r3, tg, color="goldenrod", width=barWidth, edgecolor="black", label="Trypsin+Glu-C")

        # xticks in the middle of group bars
        plt.xlabel("Unique Peptides", fontweight="bold")
        plt.xticks([r + barWidth for r in range(len(tl))], ["ORFs with >= 1 UP", "ORFs with >= 2 UP",
                                                                 "Number of UPs per ORF / 10"])

        plt.title("Unique Peptides when using a combination of different proteases")
        plt.legend()
        plt.savefig("Protease_combination_ups.png")

"""
# genome_trypsin = Digested("/home/eduardo/Documents/smORFs_SMEG/MSMEG/genome_database_no_anno_trypsin/", "genome",
#                            "pipeline/genome_database_no_anno.fasta", "Trypsin", "MSMEG_proteome.fasta")
# genome_trypsin.virtual_coverage()
# genome_lysc = Digested("/home/eduardo/Documents/smORFs_SMEG/MSMEG/genome_database_no_anno_lys_c/", "genome",
#                         "pipeline/genome_database_no_anno.fasta", "Lys_C", "MSMEG_proteome.fasta")
# genome_lysc.virtual_coverage()
# 
# genome_lysc.plot_data()
# genome_trypsin.virtual_ups()
# genome_lysc.virtual_ups()
# genome_trypsin.ups_per_orf()
# genome_lysc.ups_per_orf()
# 
# up_enzymes = PlotData("Trypsin,Lys_C")
# up_enzymes.plot_up_data()
# rna_digest = Digestion("transcriptome_database_no_anno.fasta")
# enzymes = [0, 1, 2, 3]
# for i in enzymes:
#     rna_digest.digest_orfs(i)
"""
#
# genome_trypsin = Digested("/home/eduardo/Documents/smORFs_SMEG/MSMEG/genome_database_no_anno_trypsin/", "genome",
#                           "pipeline/genome_database_no_anno.fasta", "Trypsin", "MSMEG_proteome.fasta")
# genome_lysc = Digested("/home/eduardo/Documents/smORFs_SMEG/MSMEG/genome_database_no_anno_lys_c/", "genome",
#                        "pipeline/genome_database_no_anno.fasta", "Lys_C", "MSMEG_proteome.fasta")
# genome_argc = Digested("/home/eduardo/Documents/smORFs_SMEG/MSMEG/genome_database_no_anno_arg_c/", "genome",
#                        "pipeline/genome_database_no_anno.fasta", "Arg_C", "MSMEG_proteome.fasta")
# genome_gluc = Digested("/home/eduardo/Documents/smORFs_SMEG/MSMEG/genome_database_no_anno_glu_c/", "genome",
#                        "pipeline/genome_database_no_anno.fasta", "Glu_C", "MSMEG_proteome.fasta")
#
# enzyme_list = [genome_argc, genome_gluc]
# for j in range(len(enzyme_list)):
#     enz = enzyme_list[j]
#     enz.virtual_digestion()
#     enz.virtual_coverage()
#     # enz.virtual_ups()
#     enz.ups_per_orf() # do not use self.virtual_ups as ups_per_orf calls this command already
#
# print("It took ", time.time() - start, "seconds to execute this Script.")

# cov_enzymes = PlotData("Trypsin,Lys_C,Arg_C,GluC")
# # cov_enzymes.plot_coverage_data()
# cov_enzymes.plot_up_data()