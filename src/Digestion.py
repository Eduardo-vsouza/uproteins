import pandas as pd
import glob
import numpy as np
import os
from matplotlib import pyplot as plt
from Bio import SeqIO
from tqdm import tqdm


class Digestion(object):
    def __init__(self, database):
        self.database = database

    def digest_orfs(self, enzyme):
        cmd_digest = 'perl protein2digest.pl -in=2 -ez=%s -mc=3 %s' % (enzyme, self.database)
        os.system(cmd_digest)
        cmd_rename_dir = 'mv %s %s_%s' % (self.database[:-6], self.database[:-6], enzyme)
        os.system(cmd_rename_dir)


class Digested(object):
    def __init__(self, folder, filetype, orf_db, enzyme, proteome):
        self.folder = folder
        self.filetype = filetype
        self.orf_db = orf_db
        self.enzyme = enzyme
        self.proteome = proteome

    def virtual_digestion(self):
        if not os.path.exists("Virtual_digestion/"):
            cmd_dir = 'mkdir Virtual_digestion'
            os.system(cmd_dir)

        files = glob.glob("%s/*PROT*" % self.folder)

        the_Dic = {

        }
        for i in range(10):
            with open(files[i], 'r') as name:
                names = []
                lines = name.readlines()
                for line in range(len(lines)):
                    lines[line] = lines[line].rstrip()
                    lines[line] = lines[line][3:]
                names.append(lines[1])
                the_Dic[names[0]] = ""
            df = pd.read_csv(files[i], sep='\t', header=3, keep_default_na=False, na_values=['_'])
            pepseq = df["peptide sequence"].tolist()
            for j in range(len(pepseq)):
                if 5 <= len(pepseq[j]) <= 20:
                    the_Dic[names[0]] += "%s," % pepseq[j]
        return the_Dic

    def virtual_coverage(self):
        """ Verifies the maximum obtainable coverage for each ORF when using a specific protease for protein digestion.
        It then writes the coverage into a tab-separated file for it to be plotted in the next step. """
        ids = []
        seqs = []
        parser = SeqIO.parse(self.orf_db, 'fasta')
        for entry in parser:
            ids.append(entry.id)
            seqs.append(entry.seq)

        the_Dic = self.virtual_digestion()
        all_coverage = []
        final_df = pd.DataFrame()
        for entry in the_Dic:
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
                orf_ids.append(record.id)
                orf_seqs.append(record.seq)
            ref_records = SeqIO.parse(anno, 'fasta')
            for record in ref_records:
                anno_ids.append(record.id)
                anno_seqs.append(record.seq)

        # Add info about the ORF loci in the genome to database info
        for i in range(len(orf_ids)):
            icon = orf_ids[i].find("_-_")
            beginning = orf_ids[i].find("[")
            ending = orf_ids[i].find("]")
            locus_start = orf_ids[i][beginning + 1:icon]
            locus_end = orf_ids[i][icon + 3:ending]
            orf_locus_start.append(locus_start)
            orf_locus_end.append(locus_end)

        # Check if the peptide identified by MS is unique or not
        # timer = 0
        for pep in range(len(peptides)):
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
                for i in range(len(orf_seqs)):
                    if peptides[pep] in orf_seqs[i]:
                        entries.append(orf_ids[i])
                        entries_start.append(orf_locus_start[i])
                        entries_end.append(orf_locus_end[i])
                        appearances += 1
                for i in range(len(anno_seqs)):
                    if peptides[pep] in anno_seqs[i]:
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
        tryp_cov = trypsin_df["Coverages"].tolist()
        lys_cov = lys_df["Coverages"].tolist()
        labels = ["Trypsin", "Lys-C"]
        tryp_mean = np.mean(tryp_cov)
        lys_mean = np.mean(lys_cov)
        means = [tryp_mean, lys_mean]
        print(means)
        plt.title("Coverage after digestion by different enzymes")
        plt.ylabel("Coverage")
        plt.xlabel("Enzyme")
        plt.bar(labels, means)
        plt.show()

    def plot_up_data(self):
        enzyme_list = self.enzymes.split(",")
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
        labels = ["Trypsin 1 UP", "Lys C 1 UP", "Trypsin 2 UP", "Lys C 2 UP", "Trypsin Total UPs", "Lys C Total UPs"]
        values = [tryp_orf_number_1, lysc_orf_number_1, tryp_orf_number_2, lysc_orf_number_2, tryp_up_number,
                  lysc_up_number]
        plt.title("Unique Peptides per enzyme")
        plt.ylabel("Number of UPs")
        plt.xlabel("Enzymes")
        plt.bar(labels, values)
        plt.show()


# genome_trypsin = Digested("/home/eduardo/Documents/smORFs_SMEG/MSMEG/genome_database_no_anno_trypsin/", "genome",
#                            "pipeline/genome_database_no_anno.fasta", "Trypsin", "MSMEG_proteome.fasta")
# genome_trypsin.virtual_coverage()
# genome_lysc = Digested("/home/eduardo/Documents/smORFs_SMEG/MSMEG/genome_database_no_anno_lys_c/", "genome",
#                         "pipeline/genome_database_no_anno.fasta", "Lys_C", "MSMEG_proteome.fasta")
# genome_lysc.virtual_coverage()

# genome_lysc.plot_data()
# genome_trypsin.virtual_ups()
# genome_lysc.virtual_ups()
# genome_trypsin.ups_per_orf()
# genome_lysc.ups_per_orf()

# up_enzymes = PlotData("Trypsin,Lys_C")
# up_enzymes.plot_up_data()
# rna_digest = Digestion("transcriptome_database_no_anno.fasta")
# enzymes = [0, 1, 2, 3]
# for i in enzymes:
#     rna_digest.digest_orfs(i)

#
# genome_trypsin = Digested("/home/eduardo/Documents/smORFs_SMEG/MSMEG/genome_database_no_anno_trypsin/", "genome",
#                           "pipeline/genome_database_no_anno.fasta", "Trypsin", "MSMEG_proteome.fasta")
# genome_lysc = Digested("/home/eduardo/Documents/smORFs_SMEG/MSMEG/genome_database_no_anno_lys_c/", "genome",
#                           "pipeline/genome_database_no_anno.fasta", "Lys_C", "MSMEG_proteome.fasta")
# genome_argc = Digested("/home/eduardo/Documents/smORFs_SMEG/MSMEG/genome_database_no_anno_arg_c/", "genome",
#                           "pipeline/genome_database_no_anno.fasta", "Arg_C", "MSMEG_proteome.fasta")
# genome_gluc = Digested("/home/eduardo/Documents/smORFs_SMEG/MSMEG/genome_database_no_anno_glu_c/", "genome",
#                           "pipeline/genome_database_no_anno.fasta", "Glu_C", "MSMEG_proteome.fasta")
#
# enzyme_list = [genome_trypsin]
# for j in range(len(enzyme_list)):
#     enz = enzyme_list[j]
#     enz.virtual_digestion()
#     enz.virtual_coverage()
#     enz.virtual_ups()
#     enz.ups_per_orf()
