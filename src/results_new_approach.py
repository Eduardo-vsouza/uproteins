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
import os
import sys
from Bio.Blast import NCBIXML
from Bio import SeqIO
from tqdm import tqdm
from matplotlib import pyplot as plt
import matplotlib_venn as venn

path = sys.path[0]
blast_dir = f'{path}/dependencies/blast_for_uproteins'


def fix_frames():
    df = pd.read_csv("rna_results_with_rfs.ods", sep="\t")
    df_new = df[df["Reading Frame"] != "RF"]
    df_new.to_csv("rna_results_with_rfs.ods", sep="\t", index=False)


def find_nth(string, substring, n):
    """ Finds the nth occurrence of substring in string and returns its index. """
    start = string.find(substring)
    while start >= 0 and n > 1:
        start = string.find(substring, start + len(substring))
        n -= 1
    return start


class FileType(object):
    def __init__(self, dataframe):
        self.dataframe = dataframe

    def table(self):
        df = pd.read_csv(self.dataframe, sep='\t')
        return df

    def names(self):
        df = self.table()
        names = df["accession"].tolist()
        return names

    def numbers(self):
        names = self.names()
        orf_number = {}
        number = 0
        for i in range(len(names)):
            if names[i] not in orf_number:
                number += 1
                orf_number[names[i]] = number
        return orf_number

    def assign_numbers(self):
        numbers = self.numbers()
        df = self.table()
        names = self.names()
        orf_numbers = []
        for i in range(len(names)):
            orf_numbers.append(numbers[names[i]])
        df.insert(0, "ORF Number", orf_numbers)
        return df

    def entries(self):
        df = self.table()
        entries = df["ORF Number"].tolist()
        return entries

    def seqs(self):
        df = self.table()
        seqs = df["ORF Sequence"].tolist()
        return seqs


def merge_databases():
    blast_cmd = f'{blast_dir}blastp -query Genome/genome_database.fasta -subject' \
                ' Transcriptome/transcriptome_database.fasta -task blastp-short -qcov_hsp_perc 100 -evalue 0.001 -outfmt 5' \
                ' -out blasted_databases.xml'
    os.system(blast_cmd)
    handler = open("blasted_databases.xml", 'r')
    blast_parse = NCBIXML.parse(handler)
    intersec = []
    for record in blast_parse:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                if hsp.identities == hsp.align_length:
                    intersec.append(record.query)
    with open("Genome/genome_database.fasta", 'r') as dna, open("merged_databases.fasta", 'w') as merge:
        intersec_fasta = []
        print(intersec)
        parser = SeqIO.parse(dna, 'fasta')
        for record in parser:
            # print(record.id)
            if record.id in intersec:
                intersec_fasta.append(">" + str(record.id) + "\n" + str(record.seq) + "\n")
        merge.writelines(intersec_fasta)


class Results(object):
    def __init__(self, folder, filetype, orf_db):
        self.folder = folder
        self.filetype = filetype
        self.orf_db = orf_db

    def rename_files(self):
        if not os.path.exists("%s/filtered" % self.folder):
            cmd_make_dir = 'mkdir %s/filtered' % self.folder
            os.system(cmd_make_dir)
        cmd_move = 'mv %s/*FILTERED* %s/filtered/' % (self.folder, self.folder)
        os.system(cmd_move)
        mv_pepseq = 'mv %s/filtered/*pepseq.txt %s/' % (self.folder, self.folder)
        os.system(mv_pepseq)
        files = os.listdir("%s/filtered" % self.folder)
        for file in files:
            cmd_rename = 'mv %s/filtered/%s %s/filtered/%s' % (self.folder, file, self.folder, file[:-3] + "txt")
            os.system(cmd_rename)

    def write_results(self):
        files = os.listdir("%s/filtered" % self.folder)
        i = 0
        if not os.path.exists("%s/Results" % self.folder):
            cmd_mkdir = 'mkdir %s/Results' % self.folder
            os.system(cmd_mkdir)
        for file in files:
            df = pd.read_csv("%s/filtered/%s" % (self.folder, file), sep="\t")
            print(file)
            new_df = df[
                ["accession", "scan number(s)", "acquisitionNum", "calculatedMassToCharge", "experimentalMassToCharge",
                 "chargeState", "MS-GF:DeNovoScore", "MS-GF:EValue", "MS-GF:PepQValue", "MS-GF:QValue",
                 "MS-GF:RawScore", "MS-GF:SpecEValue", "length", "spectrumFile", "pepSeq"]]
            pep_list = df["pepSeq"].tolist()
            specs = []
            names = df["accession"].tolist()
            name_e_pep = [k + j for k, j in zip(names, pep_list)]
            for occurence in name_e_pep:
                spec_count = name_e_pep.count(occurence)
                specs.append(spec_count)
            new_df.insert(2, "Spec Counts (MS Run)", specs, allow_duplicates=True)
            new_df.to_csv('%s/Results/result_n%s.txt' % (self.folder, i), sep='\t', index=False)
            i += 1

    def merge(self):
        cmd_cat = 'cat %s/Results/*result* > %s/Results/cat.txt' % (self.folder, self.folder)
        os.system(cmd_cat)
        cmd_sort = 'sort %s/Results/cat.txt | uniq > %s/Results/Results_sorted.txt' % (self.folder, self.folder)
        os.system(cmd_sort)

    def add_orf(self):
        # creates a table containing the transcriptome entry for its corresponding genome entry
        # blast_handler = open("blasted_databases.xml", 'r')
        # records = NCBIXML.parse(blast_handler)
        # results = ["Genome\tTranscriptome\n"]
        # with open("corresponding_db.txt", 'w') as corre:
        #     for record in records:
        #         for alignment in record.alignments:
        #             for hsp in alignment.hsps:
        #                 if hsp.identities == hsp.align_length:
        #                     results.append(record.query + "\t" + alignment.hit_def + "\n")
        #     corre.writelines(results)

        # adds the entry to the results file
        df = pd.read_csv("%s/Results/Results_sorted.txt" % self.folder, sep="\t")
        orf_name = df["accession"].tolist()
        dbs_df = pd.read_csv("corresponding_db.txt", sep="\t")
        dna = dbs_df["Genome"].tolist()
        rna = dbs_df["Transcriptome"].tolist()

        intersec_entry = []
        for inde in range(len(orf_name)):
            intersec_entry.append("")

        if self.filetype == "genome":
            i = 0
            for ind in range(len(orf_name)):
                i += 1
                print((i / len(orf_name)) * 100)
                # for row_db in range(len(dna)):
                if orf_name[ind] in dna:
                    finder = dna.index(orf_name[ind])
                    # if orf_name[ind] == dna[row_db]:
                    intersec_entry[ind] = rna[finder]
                else:
                    intersec_entry[ind] = "Genome unique entry"
        elif self.filetype == "transcriptome":
            i = 0
            for ind in range(len(orf_name)):
                i += 1
                print((i / len(orf_name)) * 100)
                if orf_name[ind] in rna:
                    finder = rna.index(orf_name[ind])
                    # if orf_name[ind] == rna[row_db]:
                    intersec_entry[ind] = dna[finder]
                else:
                    intersec_entry[ind] = "Transcriptome unique entry"

        df.insert(1, "Merged_db_entry", intersec_entry)
        # filled_intersec = df.loc[df['Merged_db_entry'] != ""]
        df.to_csv("%s/Results/Results_with_orfs.txt" % self.folder, sep="\t", header=True, index=False)

    def add_orf_sequence(self):
        df = pd.read_table("%s/Results/Results_with_orfs.txt" % self.folder)
        orf_name = df["accession"].tolist()
        orf_db_ids = []
        orf_db_seqs = []
        seqs = []

        with open("%s" % self.orf_db, 'r') as db:
            parse_db = SeqIO.parse(db, 'fasta')
            for record in parse_db:
                orf_db_ids.append(str(record.id))
                orf_db_seqs.append(str(record.seq))
        i = 0
        for row in range(len(orf_name)):
            i += 1
            print((i / len(orf_name)) * 100)
            try:
                finder = orf_db_ids.index(orf_name[row])
            except ValueError:
                finder = "not found"
            seqs.append(orf_db_seqs[finder])

        df.insert(1, "ORF Sequence", seqs, allow_duplicates=True)
        df.to_csv('%s/Results/Results_with_seqs.txt' % self.folder, sep='\t', index=False)

    def total_spec_count(self):
        df = pd.read_csv("%s/Results/Results_with_seqs.txt" % self.folder, sep='\t')
        pep_list = df["pepSeq"].tolist()
        names = df["accession"].tolist()
        name_and_pep = [k + j for k, j in zip(names, pep_list)]

        df.insert(0, "Name+Pepseq", name_and_pep)
        df_total = df.groupby('Name+Pepseq', as_index=False)['Spec Counts (MS Run)'].sum()
        total_spec = []
        spec_run = df["Spec Counts (MS Run)"].tolist()
        total_names = df_total["Name+Pepseq"].tolist()
        total_specs_df = df["Spec Counts (MS Run)"].tolist()
        names_pepseq_original = df["Name+Pepseq"].tolist()
        runs = df["spectrumFile"].tolist()
        check_run = []
        dic_count = {

        }
        for i in range(len(spec_run)):
            if name_and_pep[i] not in dic_count:
                dic_count[name_and_pep[i]] = 0
            dic_count[name_and_pep[i]] += 1
            # if str(spec_run[i])+str(name_and_pep[i]) not in check_run:
            #
            #     if name_and_pep[i] in dic_count:
            #         dic_count[name_and_pep[i]] += 1
            #     check_run.append(str(spec_run[i])+str(name_and_pep[i]))
        for i in range(len(spec_run)):
            total_spec.append(dic_count[name_and_pep[i]])
        print(dic_count)
        # print("\nPerforming Spectral Counts\n")
        # for linha in range(len(names_pepseq_original)):
        # i += 1
        # print((i / len(names_pepseq_original)) * 100)
        # if names_pepseq_original[linha] in total_names:
        #     if runs[linha] not in check_run:
        # check_run.append(runs[linha])
        # finder = total_names.index(names_pepseq_original[linha])
        # total_spec.append(total_specs_df[finder])
        # for row in total_names:
        #     ind = total_names.index(row)
        #     if row == linha:
        #         total_spec.append(total_specs_df[ind])
        df.insert(6, "Total Spec Counts", total_spec)
        df = df.drop(columns="Name+Pepseq")
        print("\nSpectral counts added to the table\n")
        df.to_csv("%s/Results/Results_with_total_specs.txt" % self.folder, sep="\t", index=False)

    def genome_coordinates(self):
        df = pd.read_csv("%s/Results/Results_with_total_specs.txt" % self.folder, sep='\t')
        if self.filetype == "genome":
            merged_entries = df["accession"].tolist()
        elif self.filetype == "transcriptome":
            merged_entries = df["Merged_db_entry"].tolist()
        coords = []
        for i in range(len(merged_entries)):
            start = merged_entries[i].find("[")
            end = merged_entries[i].find("]")
            separator = merged_entries[i].find("_-_")
            start_c = merged_entries[i][start + 1:separator]
            end_c = merged_entries[i][separator + 3:end]
            coordinates = "%s - %s" % (start_c, end_c)
            coords.append(coordinates)
        df.insert(7, "Genome Coordinates", coords)
        df.to_csv("%s/Results/Results_with_coordinates.ods" % self.folder, sep="\t", index=False)

    def check_unique(self, args):
        # spec counts data frame
        df = pd.read_csv("%s/Results/Results_with_coordinates.ods" % self.folder, sep='\t')
        peptides = df["pepSeq"].tolist()
        entries = df["accession"].tolist()
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
        with open("%s_database_no_anno.fasta" % self.filetype, 'r') as db, open(args.proteome, 'r') as anno:
            records = SeqIO.parse(db, 'fasta')
            for record in records:
                orf_ids.append(record.id)
                orf_seqs.append(record.seq)
            ref_records = SeqIO.parse(anno, 'fasta')
            for record in ref_records:
                anno_ids.append(record.id)
                anno_seqs.append(record.seq)

        # Add info about the ORF loci in the genome to database info
        if self.filetype != "transcriptome":
            for i in range(len(orf_ids)):
                icon = orf_ids[i].find("_-_")
                beginning = orf_ids[i].find("[")
                ending = orf_ids[i].find("]")
                locus_start = orf_ids[i][beginning + 1:icon]
                locus_end = orf_ids[i][icon + 3:ending]
                orf_locus_start.append(locus_start)
                orf_locus_end.append(locus_end)
        else:
            for i in range(len(orf_ids)):
                # rna_df = df[df["accession"] == orf_ids[i]]
                # print(rna_df)
                # full_rna_coods = rna_df["Genome Coordinates"].tolist()
                # locus_start = full_rna_coods[0].split(" - ")[0]
                # locus_end = full_rna_coods[0].split(" - ")[1]
                # print(locus_start)

                icon = orf_ids[i].find("_-_")
                beginning = find_nth(orf_ids[i], "[", 2)
                ending = find_nth(orf_ids[i], "]", 2)
                locus_start = orf_ids[i][beginning + 1:icon]
                locus_end = orf_ids[i][icon + 3:ending]
                orf_locus_start.append(locus_start)
                orf_locus_end.append(locus_end)

        # create a dictionary for peptides that were already checked for their uniqueness
        # this reduced execution time from more than 5 hours to less than 30 minutes when using more than 20 MS Runs
        checked_peptides = {

        }
        # Check whether the peptide identified by MS is unique or not
        timer = 0
        for pep in tqdm(range(len(peptides))):
            if peptides[pep] not in checked_peptides:
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
                        if appearances > 2:
                            break
                for i in range(len(anno_seqs)):
                    if peptides[pep] in anno_seqs[i]:
                        ref_appear += 1
                        if ref_appear > 2:
                            break
                if appearances > 1:
                    for i in range(len(entries)):
                        for j in range(len(entries_start)):
                            entries_length = len(entries)
                            if int(entries_start[entries_length - 1]) not in range(int(entries_start[j]),
                                                                                   int(entries_end[j])):
                                unique_pep = False
                                break
                elif appearances == 0 and ref_appear > 1:
                    unique_pep = False
                elif appearances >= 1 and ref_appear >= 1:
                    unique_pep = False
                unique.append(unique_pep)
                checked_peptides[peptides[pep]] = unique_pep
            else:
                unique.append(checked_peptides[peptides[pep]])

        df.insert(2, "Unique Peptide", unique, allow_duplicates=True)
        df.to_csv("%s/Results/Results_with_utps.txt" % self.folder, sep='\t', index=False)

    # def check_unique(self, args):
    #     # spec counts data frame
    #     df = pd.read_csv("%s/Results/Results_with_coordinates.ods" % self.folder, sep='\t')
    #     peptides = df["pepSeq"].tolist()
    #     unique = []
    #
    #     # database info
    #     orf_ids = []
    #     orf_seqs = []
    #     orf_locus_start = []
    #     orf_locus_end = []
    #
    #     # RefSeq info
    #     anno_ids = []
    #     anno_seqs = []
    #
    #     # Add ORF ids and seqs to the database info
    #     with open("%s_database_no_anno.fasta" % self.filetype, 'r') as db, open(args.proteome, 'r') as anno:
    #         records = SeqIO.parse(db, 'fasta')
    #         for record in records:
    #             orf_ids.append(record.id)
    #             orf_seqs.append(record.seq)
    #         ref_records = SeqIO.parse(anno, 'fasta')
    #         for record in ref_records:
    #             anno_ids.append(record.id)
    #             anno_seqs.append(record.seq)
    #
    #     # Add info about the ORF loci in the genome to database info
    #     for i in range(len(orf_ids)):
    #         icon = orf_ids[i].find("_-_")
    #         beginning = orf_ids[i].find("[")
    #         ending = orf_ids[i].find("]")
    #         locus_start = orf_ids[i][beginning + 1:icon]
    #         locus_end = orf_ids[i][icon + 3:ending]
    #         orf_locus_start.append(locus_start)
    #         orf_locus_end.append(locus_end)
    #
    #     # Check whether the peptide identified by MS is unique or not
    #     timer = 0
    #     for pep in peptides:
    #         timer += 1
    #         print((timer / len(peptides)) * 100)
    #         appearances = 0
    #         ref_appear = 0
    #         entries = []
    #         entries_start = []
    #         entries_end = []
    #         unique_pep = True
    #         for i in range(len(orf_seqs)):
    #             if pep in orf_seqs[i]:
    #                 entries.append(orf_ids[i])
    #                 entries_start.append(orf_locus_start[i])
    #                 entries_end.append(orf_locus_end[i])
    #                 appearances += 1
    #         for i in range(len(anno_seqs)):
    #             if pep in anno_seqs[i]:
    #                 ref_appear += 1
    #         if appearances > 1:
    #             for i in range(len(entries)):
    #                 for j in range(len(entries_start)):
    #                     entries_length = len(entries)
    #                     if int(entries_start[entries_length - 1]) not in range(int(entries_start[j]),
    #                                                                            int(entries_end[j])):
    #                         unique_pep = False
    #         elif appearances == 0 and ref_appear > 1:
    #             unique_pep = False
    #         elif appearances >= 1 and ref_appear >= 1:
    #             unique_pep = False
    #         unique.append(unique_pep)
    #     df.insert(4, "Unique Peptide", unique, allow_duplicates=True)
    #     df.to_csv("%s/Results/Results_with_utps.txt" % self.folder, sep='\t', index=False)

    def coverage(self):
        df = pd.read_csv("%s/Results/Results_with_utps.txt" % self.folder, sep='\t')
        proteins = df["ORF Sequence"].tolist()

        coverage_pct = []
        coverage_dic = {

        }
        for i in range(len(proteins)):
            prot = proteins[i]
            if prot not in coverage_dic:
                coverage_dic[prot] = 0
                coverage = 0
                nums = []
                new_df = df[df["ORF Sequence"] == prot]
                peptides = new_df["pepSeq"].tolist()
                for j in range(len(peptides)):
                    pep = peptides[j]
                    if pep in prot:
                        start = prot.find(pep)
                        size = len(pep)
                        end = size + start
                        for k in range(start, end):
                            if k not in nums:
                                nums.append(k)
                for m in range(len(prot)):
                    if m in nums:
                        coverage += 1
                coverage = (coverage / (len(prot))) * 100
                coverage_pct.append(coverage)
                coverage_dic[prot] = coverage
            else:
                coverage_pct.append(coverage_dic[prot])
        df.insert(4, "Coverage", coverage_pct, allow_duplicates=True)
        df.to_csv("%s/Results/Results_with_coverage.txt" % self.folder, sep='\t', index=False)

    #  OLD version
    # def create_fasta(self):
    #     df = pd.read_csv("%s/Results/Results_with_coverage" % self.folder, sep='\t')
    #     ids = df["accession"].tolist()
    #     seqs = df["ORF Sequence"].tolist()
    #     entries = []
    #     with open("%s/Results/%s_ORFs.fasta" % (self.folder, self.filetype), 'w') as fa:
    #         for i in range(len(ids)):
    #             entries.append(">" + ids[i] + "\n" + seqs[i] + "\n")
    #         fa.writelines(entries)

    def summarize_results(self):
        df = pd.read_csv("%s/Results/Results_with_coverage.txt" % self.folder, sep='\t')
        summarized = df.loc[df['Unique Peptide'] == True]
        if self.filetype == "transcriptome":
            summarized = summarized.loc[summarized['Merged_db_entry'] != "Transcriptome unique entry"]
        if self.filetype == "genome":
            summarized = summarized.loc[summarized['Merged_db_entry'] != "Genome unique entry"]
        summarized.to_csv("%s/Results/%s_summarized_final_results.xlsl" % (self.folder, self.filetype), sep='\t',
                          index=False)

    def unique_orfs(self):
        df = pd.read_csv("%s/Results/%s_summarized_final_results.xlsl" % (self.folder, self.filetype), sep="\t")
        orf_names = df["accession"].tolist()
        seqs = df["ORF Sequence"].tolist()

        check_entry = []
        entries = []

        for i in range(len(orf_names)):
            if orf_names[i] not in check_entry:
                check_entry.append(orf_names[i])
                entries.append(">" + orf_names[i] + "\n" + seqs[i] + "\n")

        with open("%s/Results/%s_unique_ORFs_both.fasta" % (self.folder, self.filetype), 'w') as out:
            out.writelines(entries)
        if self.filetype == "genome":
            df = pd.read_csv("Genome/Results/both_summarized_final_results.xlsl", sep="\t")

            orf_names = df["accession"].tolist()
            seqs = df["ORF Sequence"].tolist()

            check_entry = []
            entries = []

            for i in range(len(orf_names)):
                if orf_names[i] not in check_entry:
                    check_entry.append(orf_names[i])
                    entries.append(">" + orf_names[i] + "\n" + seqs[i] + "\n")

            with open("%s/Results/both_ORFs.fasta" % (self.folder), 'w') as out:
                out.writelines(entries)

    def summarize_dna_rna(self):
        """ Creates summarized results for ORFs uniquely found on Genome or Transcriptome databases. These results
        contain only ORFs with Unique Peptides. """
        df = pd.read_csv("%s/Results/Results_with_coverage.txt" % self.folder, sep="\t")
        summarized = df.loc[df['Unique Peptide'] == True]
        if self.filetype == "transcriptome":
            summarized = summarized.loc[summarized['Merged_db_entry'] == "Transcriptome unique entry"]
        elif self.filetype == "genome":
            summarized = summarized.loc[summarized['Merged_db_entry'] == "Genome unique entry"]

        if self.filetype == "transcriptome":
            summarized = summarized[summarized["accession"].str.match("TRINITY")]
        # make sure the entries in the reference genome start with NC_
        elif self.filetype == "genome":
            summarized = summarized[summarized["accession"].str.match("NC_")]

        summarized.to_csv("%s/Results/%s_unique_results_summarized.xlsl" % (self.folder, self.filetype),
                          sep='\t', index=False)

    # @staticmethod
    # def add_names():
    #     dna = FileType("Genome/Results/genome_unique_results_summarized.xlsl")
    #     both = FileType("Genome/Results/both_summarized_final_results.ods")
    #     rna = FileType("Transcriptome/Results/transcriptome_unique_results_summarized.xlsl")
    #     dna_df = dna.assign_numbers()
    #     both_df = both.assign_numbers()
    #     rna_df = rna.assign_numbers()
    #     rna_df_new = rna_df[rna_df["Genome Coordinates"] != "Not found - Not found"]
    #     dna_df.to_csv("Genome/Results/Genome_unique_results_numbers.xlsl", sep='\t', index=False)
    #     both_df.to_csv("Genome/Results/Both_results_numbers.xlsl", sep='\t', index=False)
    #     rna_df_new.to_csv("Transcriptome/Results/Rna_unique_results_numbers.xlsl", sep='\t', index=False)

    # def create_fast(self, subset):
    #     """ Subset must be either Genome, Both or Rna. """
    #     if subset == "Both":
    #         orfs = FileType("%s/Results/%s_results_numbers.xlsl" % (self.folder, subset))
    #     else:
    #         orfs = FileType("%s/Results/%s_unique_results_numbers.xlsl" % (self.folder, subset))
    #
    #     entries = orfs.entries()
    #     seqs = orfs.seqs()
    #     check = []
    #     fasta = []
    #
    #     for i in range(len(entries)):
    #         if entries[i] not in check:
    #             check.append(entries[i])
    #             fasta.append(">ORF " + str(entries[i]) + "\n" + seqs[i] + "\n")
    #
    #     with open("%s/Results/%s_ORFs.fasta" % (self.folder, subset), 'w') as outfile:
    #         outfile.writelines(fasta)

    # @staticmethod
    def merge_both_results(self):
        """ Uses the strategy employed to reduce ORF loss due to decoy .
        """
        pei = self.filetype
        genome_df = pd.read_csv("Genome/Results/genome_summarized_final_results.xlsl", sep="\t")
        rna_df = pd.read_csv("Transcriptome/Results/transcriptome_summarized_final_results.xlsl", sep="\t")

        dna_names = genome_df["accession"].tolist()
        rna_names = rna_df["accession"].tolist()
        rna_merged = rna_df["Merged_db_entry"].tolist()

        columns = []
        indexes = []
        rna_entries = []
        new_df = pd.DataFrame()
        new_df = new_df.append(genome_df, ignore_index=True)
        new_df = new_df.append(rna_df, ignore_index=True)
        # genome_df.append(rna_df)
        new_df.to_csv("Genome/Results/both_summarized_final_results.xlsl", sep="\t", index=False)

        # for i in tqdm(range(len(rna_names))):
        #     if rna_merged[i] in dna_names:
        #         indexes.append(i)

        # for i in range(len(rna_df.head(0))):
        #     genome_df.append

        # print(rna_df.loc[rna_df['Merged_db_entry'] in genome_df["accession"]])

        # print(new_df)

    def minimal_results(self):
        if self.filetype == "Both":
            df = pd.read_csv("%s/Results/%s_results_final.xlsl" % (self.folder, self.filetype), sep="\t")
            accession = df["accession"].tolist()
            seqs = df["ORF Sequence"].tolist()
            coords = df["Genome Coordinates"].tolist()
            coverage = df["Coverage"].tolist()
            specs = df["Total Spec Counts"].tolist()
            new_df = pd.DataFrame()
            new_df.insert(0, "accession", accession)
            new_df.insert(1, "ORF Sequence", seqs)
            new_df.insert(2, "Genome Coordinates", coords)
            new_df.insert(3, "Coverage", coverage)
            new_df.insert(4, "Total Spec Counts", specs)
            new_df = new_df.drop_duplicates()
            new_df.to_csv("%s/Results/%s_results_compact.xlsl" % (self.folder, self.filetype), sep="\t")
        else:
            df = pd.read_csv("%s/Results/%s_unique_results_summarized.xlsl" % (self.folder, self.filetype), sep="\t")
            accession = df["accession"].tolist()
            seqs = df["ORF Sequence"].tolist()
            coords = df["Genome Coordinates"].tolist()
            coverage = df["Coverage"].tolist()
            specs = df["Total Spec Counts"].tolist()
            new_df = pd.DataFrame()
            new_df.insert(0, "accession", accession)
            new_df.insert(1, "ORF Sequence", seqs)
            new_df.insert(2, "Genome Coordinates", coords)
            new_df.insert(3, "Coverage", coverage)
            new_df.insert(4, "Total Spec Counts", specs)
            new_df = new_df.drop_duplicates()
            new_df.to_csv("%s/Results/%s_results_compact.xlsl" % (self.folder, self.filetype), sep="\t", index=False)

    def minimum_runs_compact(self, runs, min_type):
        df = pd.read_csv("%s_minimum_%s_runs.xlsl" % (min_type, runs), sep="\t")
        accession = df["accession"].tolist()
        seqs = df["ORF Sequence"].tolist()
        coords = df["Genome Coordinates"].tolist()
        coverage = df["Coverage"].tolist()
        specs = df["Total Spec Counts"].tolist()
        new_df = pd.DataFrame()
        new_df.insert(0, "accession", accession)
        new_df.insert(1, "ORF Sequence", seqs)
        new_df.insert(2, "Genome Coordinates", coords)
        new_df.insert(3, "Coverage", coverage)
        new_df.insert(4, "Total Spec Counts", specs)
        new_df = new_df.drop_duplicates()
        if not os.path.exists("Results"):
            cmd_dir = 'mkdir Results'
            os.system(cmd_dir)
        new_df.to_csv("Results/%s_minimum_%s_runs_compact.xlsl" % (min_type, runs), sep="\t", index=False)


class GenomicContext(object):
    def __init__(self, args, db, folder):
        """ db parameter corresponds to the set of ORFs being worked on, i.e "both", "genome" or "transcriptome" ORFs.
        """
        self.args = args
        self.df = "{}/Results/{}_summarized_final_results.xlsl".format(folder, db)
        self.gb = args.genbank
        self.products = []
        self.genes = []
        self.coords = []
        self.hits = {}
        self.db = db
        pass

    def get_info(self):
        df = pd.read_csv(self.df, sep="\t")
        names = df["accession"].tolist()
        seqs = df["ORF Sequence"].tolist()
        coords = df["Genome Coordinates"].tolist()
        return names, seqs, coords

    def parse_genbank(self):
        gb_file = self.gb
        for seq_record in SeqIO.parse(gb_file, "genbank"):
            for feature in seq_record.features:
                if feature.type != "CDS" and feature.type != "gene":
                    self.genes.append(str(feature.type))
                    if feature.qualifiers.get("note") is not None:
                        self.products.append(str(feature.qualifiers.get("note")[0]))
                    else:
                        self.products.append(str("Product not found"))
                    if feature.location is not None:
                        self.coords.append(str(feature.location)[1:-4])
                    else:
                        self.coords.append("Location not found")

    def match_features(self):
        """ See if the identified ORFs match to the location of a previously annotated, non-coding feature.
        """
        names, seqs, coords = self.get_info()
        self.parse_genbank()
        for i in range(len(coords)):
            if coords[i] != "Not found - Not found":
                data_locus = coords[i].split(" - ")
                my_start = int(data_locus[0])
                my_end = int(data_locus[1])
                if int(my_start) > int(my_end):
                    my_start = int(data_locus[1])
                    my_end = int(data_locus[0])
                # print(my_start, my_end)
                check = 0
                for j in range(len(self.coords)):
                    locus = self.coords[j].split(":")
                    if self.products[j] != "Product not found":
                        if len(locus) < 3:
                            # print(locus)
                            start = int(locus[0])
                            end = int(locus[1])
                            # print(start, end)
                            if start in range(my_start, my_end) or end in range(my_start, my_end):
                                check += 1
                            if my_start in range(start, end) or my_end in range(start, end):
                                check += 1
                            if check > 0:
                                if names[i] not in self.hits:
                                    self.hits[names[i]] = "{3}, {2}, {0} - {1}".format(start, end, self.products[j],
                                                                                       self.genes[j])
                                else:
                                    self.hits[names[i]] += "|{3}, {2}, {0} - {1}".format(start, end, self.products[j],
                                                                                         self.genes[j])
                            check = 0
                if check < 1:
                    if names[i] not in self.hits:
                        self.hits[names[i]] = "No feature"
                    # else:
                    #     self.hits[names[i]] += "No feature"
            else:
                self.hits[names[i]] = "No feature"

    def add_to_df(self):
        names, seqs, coords = self.get_info()
        feature_product = []

        for i in range(len(names)):
            feature_product.append(self.hits[names[i]])
        new_df = pd.DataFrame()
        new_df.insert(0, "accession", names)
        new_df.insert(1, "ORF Sequence", seqs)
        new_df.insert(2, "Genome Coordinates", coords)
        new_df.insert(3, "Genomic Feature", feature_product)
        new_df.to_csv("{}_orfs_and_features.xlsl".format(self.db), sep="\t", index=False)


class TranscriptLocalizer(object):
    def __init__(self, args):
        self.args = args

    def find_transcript(self):
        genome = self.args.genome
        transcriptome = "trinity/Trinity.fasta"
        cmd = 'nucmer %s %s' % (genome, transcriptome)
        os.system(cmd)
        cmd_coords = 'show-coords -r -c -l -T out.delta > Transcriptome/Results/nucmer.coords'
        os.system(cmd_coords)

    def fix_nucmer(self):
        with open("Transcriptome/Results/nucmer.coords", 'r') as nucmer, open("Transcriptome/Results/nucmer_fixed.txt",
                                                                              'w') as fixed:
            lines = nucmer.readlines()
            new_lines = []
            for line in lines:
                if line.startswith("[S1]"):
                    line = line.rstrip()
                    line += "\t[ENTRY]\n"
                    new_lines.append(line)
                else:
                    new_lines.append(line)
            fixed.writelines(new_lines)

    def find_coords(self):
        df = pd.read_csv("Transcriptome/Results/nucmer_fixed.txt", sep="\t", header=2)
        start_rna = df["[S1]"].tolist()
        end_rna = df["[E1]"].tolist()
        tags = df["[ENTRY]"].tolist()
        start_in_rna = df["[S2]"].tolist()
        end_in_rna = df["[E2]"].tolist()
        coord_dic = {

        }
        coord_rna_dic = {}

        for i in range(len(tags)):
            if tags[i] not in coord_dic:
                coord_dic[tags[i]] = "%s - %s" % (start_rna[i], end_rna[i])
                coord_rna_dic[tags[i]] = "%s - %s" % (start_in_rna[i], end_in_rna[i])

        return coord_dic, coord_rna_dic

    def add_coords(self):
        coord_dic, coord_rna = self.find_coords()
        df = pd.read_csv("Transcriptome/Results/Results_with_coordinates.ods", sep="\t")
        df = df[df["accession"].str.match("TRINITY")]
        entries = df["accession"].tolist()
        coords_in_rna = df["Genome Coordinates"].tolist()

        # info about genome coordinates of the RNAs
        rna_starts = []
        rna_ends = []
        in_rna_start = []
        in_rna_end = []

        for i in range(len(entries)):
            finder = find_nth(entries[i], "_", 5)
            entry = entries[i][:finder]
            if entry in coord_dic:
                locus = coord_dic[entry].split(" - ")
                start = locus[0]
                end = locus[1]
                rna_starts.append(start)
                rna_ends.append(end)
            if entry in coord_rna:
                locus = coord_rna[entry].split(" - ")
                start = locus[0]
                end = locus[1]
                in_rna_start.append(start)
                in_rna_end.append(end)
            else:
                rna_starts.append("Not found")
                rna_ends.append("Not found")
                in_rna_start.append("Not found")
                in_rna_end.append("Not found")
        # info about orfs in the rnas
        orf_rna_starts = []
        orf_rna_ends = []

        for i in range(len(coords_in_rna)):
            # locus = coords_in_rna[i].split(" - ")
            finder = find_nth(entries[i], "[", 2)
            finder_end = find_nth(entries[i], "]", 2)
            coordinates = entries[i][finder + 1:finder_end]

            locus = coordinates.split("_-_")

            rna_orf_start = locus[0]
            rna_orf_end = locus[1]
            orf_rna_starts.append(rna_orf_start)
            orf_rna_ends.append(rna_orf_end)

        # info about orfs in the genome
        orf_gen_starts = []
        orf_gen_ends = []
        genome_coords = []
        for i in range(len(entries)):
            if rna_starts[i] != "Not found":

                if in_rna_start[i] > in_rna_end[i]:
                    orf_genome_start = int(rna_ends[i]) - int(orf_rna_ends[i])
                    orf_genome_end = orf_genome_start + (int(orf_rna_ends[i]) - int(orf_rna_starts[i]))
                    orf_gen_ends.append(orf_genome_start + 1)
                    orf_gen_starts.append(orf_genome_end + 1)
                elif in_rna_start[i] < in_rna_end[i]:
                    orf_genome_start = int(rna_starts[i]) + int(orf_rna_starts[i])
                    orf_genome_end = orf_genome_start + (int(orf_rna_ends[i]) - int(orf_rna_starts[i]))
                    orf_gen_ends.append(orf_genome_end - 1)
                    orf_gen_starts.append(orf_genome_start - 1)
            else:
                orf_gen_starts.append("Not found")
                orf_gen_ends.append("Not found")

        for i in range(len(orf_gen_starts)):
            genome_coords.append("%s - %s" % (orf_gen_starts[i], orf_gen_ends[i]))
        # create table
        df = df.drop(columns="Genome Coordinates")
        df.insert(9, "Genome Coordinates", genome_coords)
        df.to_csv("Transcriptome/Results/Results_with_coordinates.ods", sep="\t", index=False)


def create_venn():
    # # Transcriptome data
    # rna_df = pd.read_csv("Transcriptome/Results/Results_with_utps.txt", sep='\t')
    # rna_access = rna_df["accession"].tolist()
    # rna_seqs = rna_df["ORF Sequence"].tolist()
    # rna_unique = rna_df["Unique Peptide"].tolist()
    # rna_novel = []
    #
    # for i in range(len(rna_access)):
    #     if rna_access[i].startswith("TRINITY"):
    #         if rna_unique[i] == True:
    #             rna_novel.append(rna_seqs[i])
    #
    # rna_novel_df = pd.DataFrame()
    # rna_novel_df.insert(0, "ORF", rna_novel)
    # rna_novel_df.to_csv("Transcriptome/Results/ORFs.txt", sep='\t', index=False)
    # cmd_rna_unique = 'sort Transcriptome/Results/ORFs.txt | uniq > Transcriptome/Results/Uniq_orfs.txt'
    # os.system(cmd_rna_unique)
    # unique_rna_seqs = []
    # with open("Transcriptome/Results/Uniq_orfs.txt", 'r') as rna:
    #     for line in rna:
    #         unique_rna_seqs.append(line)
    #
    # # Genome data
    # dna_df = pd.read_csv("Genome/Results/Results_with_utps.txt", sep='\t')
    # dna_access = dna_df["accession"].tolist()
    # dna_seqs = dna_df["ORF Sequence"].tolist()
    # dna_uniq = dna_df["Unique Peptide"].tolist()
    # dna_novel = []
    #
    # for i in range(len(dna_access)):
    #     if dna_access[i].startswith("NC_"):
    #         if dna_uniq[i] == True:
    #             dna_novel.append(dna_seqs[i])
    #
    # dna_novel_df = pd.DataFrame()
    # dna_novel_df.insert(0, "ORF", dna_novel)
    # dna_novel_df.to_csv("Genome/Results/ORFs.txt", sep='\t', index=False)
    #
    # cmd_dna_unique = 'sort Genome/Results/ORFs.txt | uniq > Genome/Results/Uniq_orfs.txt'
    # os.system(cmd_dna_unique)
    # unique_dna_seqs = []
    # with open("Genome/Results/Uniq_orfs.txt", 'r') as dna:
    #     for line in dna:
    #         unique_dna_seqs.append(line)

    # transcriptome data
    rna_df = pd.read_csv("Transcriptome/Results/transcriptome_unique_results_summarized.xlsl", sep="\t")
    rna_orfs = rna_df["ORF Sequence"].tolist()
    rna_unique_orfs = []
    for i in range(len(rna_orfs)):
        if rna_orfs[i] not in rna_unique_orfs:
            rna_unique_orfs.append(rna_orfs[i])

    # genome data
    dna_df = pd.read_csv("Genome/Results/genome_unique_results_summarized.xlsl", sep="\t")
    dna_orfs = dna_df["ORF Sequence"].tolist()
    dna_unique_orfs = []
    for i in range(len(dna_orfs)):
        if dna_orfs[i] not in dna_unique_orfs:
            dna_unique_orfs.append(dna_orfs[i])

    # both data
    both_df = pd.read_csv("Genome/Results/Both_results_final.xlsl", sep="\t")
    both_orfs = both_df["ORF Sequence"].tolist()
    both_unique = []
    for i in range(len(both_orfs)):
        if both_orfs[i] not in both_unique:
            both_unique.append(both_orfs[i])
    for i in range(len(both_unique)):
        rna_unique_orfs.append(both_unique[i])
        dna_unique_orfs.append(both_unique[i])

    # Building venn diagram
    venn.venn2(subsets=[set(rna_unique_orfs), set(dna_unique_orfs)], set_labels=('Transcriptome ORFs', 'Genome ORFs'))
    plt.title('Shared ORFs')
    plt.savefig("venn.pdf")

    # def create_fasta(self):
    #     with open("%s" % self.orf_db, 'r') as db, \
    #             open("%s/Results/results_with_seqs.txt" % self.folder, 'r') as results, \
    #             open("%s/Results/%s_novel_orfs.fasta" % (self.folder, self.filetype), 'w') as fna, \
    #             open("%s/Results/%s_orf_subset" % (self.folder, self.filetype), 'w') as subset:
    #         parse_db = SeqIO.parse(db, 'fasta')
    #         df = pd.read_table(results)
    #         orf_name = df["accession"].tolist()
    #         sequences = []
    #         seq_for_venn = []
    #         for record in parse_db:
    #             if record.id in orf_name:
    #                 sequences.append(">" + str(record.id) + "\n")
    #                 sequences.append(str(record.seq) + "\n")
    #                 seq_for_venn.append(str(record.seq) + "\n")
    #         fna.writelines(sequences)
    #         subset.writelines(seq_for_venn)
    #
    # # def blast_orfs(self, genome):
    # #     """ Aligns the fasta file containing the ORFs that matched to peptides identified in the MS run
    # #     to the genome. Only for genome results. """
    # #     Cmd_matched_pep = '/home/farminf/Downloads/ncbi-blast-2.9.0+/bin/tblastn -query ' \
    # #                       '%s/Results/%s_novel_orfs.fasta -subject %s -db_gencode 11 -evalue 0.01 -qcov_hsp_perc 100' \
    # #                       ' -outfmt 5 -out %s/Results/%s_orfs_blasted_to_genome.xml' \
    # #                       % (self.folder, self.filetype, genome, self.folder, self.filetype)
    # #     os.system(Cmd_matched_pep)
    #
    # def add_genome_coordinates(self):
    #     """ Do not use on Transcriptome results. These need to be adressed differently, as isoforms may be present in
    #     the results. """
    #     df = pd.read_csv("%s/Results/Results_with_total_specs.txt" % self.folder, sep='\t')
    #     name = df["accession"].tolist()
    #     name_subs = df["accession"].tolist()
    #     start = []
    #     end = []
    #     i = 0
    #     for row in name:
    #         start.append("Not found")
    #         end.append("Not found")
    #     for row in name:
    #         i += 1
    #         print((i / len(name)) * 100)
    #         ind = name.index(row)
    #         result_handle = open("%s/Results/%s_orfs_blasted_to_genome.xml" % (self.folder, self.filetype), 'r')
    #         blast_records = NCBIXML.parse(result_handle)
    #         for record in blast_records:
    #             if row == record.query:
    #                 if record.alignments:
    #                     for alignment in record.alignments:
    #                         for hsp in alignment.hsps:
    #                             if hsp.identities == hsp.align_length:
    #                                 start[ind] = str(hsp.sbjct_start)
    #                                 end[ind] = str(hsp.sbjct_end)
    #     print(len(start))
    #     print(start)
    #     df.insert(6, "Start", start)
    #     df.insert(7, "End", end)
    #     df.to_csv("%s/Results/Results_with_coordinates_tblastn" % self.folder, sep='\t', index=False)
    #
    # def rna_coordinates(self):
    #     df = pd.read_csv("/home/farminf/Eduardo/smORFs/protopype/rna_mapping", delim_whitespace=True)
    #     coord = df.iloc[:, 8].tolist()
    #     old_name = df.iloc[:, 0].tolist()
    #     name = []
    #     for item in old_name:
    #         item = item.replace(">", "")
    #         name.append(item)
    #     results = pd.read_csv("Transcriptome/Results/Results_with_total_specs.txt", sep='\t')
    #     orf_name = results["accession"].tolist()
    #     new_orf_name = []
    #
    #     # Get genome coordinates for RNA #
    #     for row in orf_name:
    #         _ind = find_nth(row, "_", 5)
    #         new_orf_name.append(row[:_ind])
    #     genome_coord = []
    #     for row in new_orf_name:
    #         if row in name:
    #             local = name.index(row)
    #             if row == name[local]:
    #                 genome_coord.append(coord[local])
    #         elif row not in name:
    #             genome_coord.append("Not found")
    #
    #     # RNA Coordinates in the genome separated by start and end #
    #     start = []
    #     end = []
    #     for item in genome_coord:
    #         item = item.split("..")
    #         start.append(item[0])
    #         end.append(item[1])
    #
    #     # ORF coordinates in the transcript #
    #     orf_start = []
    #     orf_end = []
    #     for row in orf_name:
    #         first_ind = find_nth(row, "_", 6)
    #         second_ind = find_nth(row, "_", 9)
    #         coordinates = row[first_ind + 2:second_ind - 1]
    #         subordinates = coordinates.split("_-_")
    #         orf_start.append(subordinates[0])
    #         orf_end.append(subordinates[1])
    #
    #     # ORF Coordinates in the genome #
    #     orf_genome_start = []
    #     orf_genome_end = []
    #     for row in orf_start:
    #         ind = orf_start.index(row)
    #         genstart = int(start[ind]) + int(row)
    #         genend = int(start[ind]) + int(orf_end[ind])
    #         orf_genome_start.append(genstart)
    #         orf_genome_end.append(genend)
    #     print(orf_genome_start)
    #     print(orf_genome_end)
    #
    #     # Add orf coordinates in the genome to the results table #
    #     results.insert(6, "Start", orf_genome_start)
    #     results.insert(7, "End", orf_genome_end)
    #     results.to_csv("%s/Results/Results_with_coordinates.txt" % self.folder, sep="\t", index=False)
    #
    # def summarized_results(self):
    #     df = pd.read_csv("%s/Results/Results_with_coordinates_tblastn" % self.folder, sep="\t")
    #     df = df.drop(columns="Spec Counts (MS Run)")
    #     df = df.drop(columns="spectrumFile")
    #     df.to_csv("%s/Results/Summarized_results.txt" % self.folder, sep="\t", index=False)
    #     cmd_uniq = 'sort %s/Results/Summarized_results.txt | uniq > %s/Results/Summarized_results_unique.txt' \
    #                % (self.folder, self.folder)
    #     os.system(cmd_uniq)
    #     new_df = pd.DataFrame({})
    #     accession = df["accession"].tolist()
    #     orf_name = df["ORF Sequence"].tolist()
    #     pepseq = df["pepSeq"].tolist()
    #     spec = df["Total Spec Counts"].tolist()
    #     starter = df["Start"].tolist()
    #     ender = df["End"].tolist()
    #
    #     # new set #
    #     name = []
    #     accs = []
    #     seq = []
    #     specs = []
    #     start = []
    #     end = []
    #
    #     for row in starter:
    #         ind = starter.index(row)
    #         if row != "Not found":
    #             name.append(accession[ind])
    #             accs.append(orf_name[ind])
    #             seq.append(pepseq[ind])
    #             specs.append(spec[ind])
    #             start.append(starter[ind])
    #             end.append(ender[ind])
    #
    #     new_df.insert(0, "accession", name)
    #     new_df.insert(1, "ORF Sequence", accs)
    #     new_df.insert(2, "pepSeq", seq)
    #     new_df.insert(3, "Total Spec Counts", specs)
    #     new_df.insert(4, "Start", start)
    #     new_df.insert(5, "End", end)
    #
    #     new_df.to_csv("%s/Results/Summarized_found_results.txt" % self.folder, sep="\t", index=False)
    #     cmd_uniq = 'sort %s/Results/Summarized_found_results.txt | uniq > %s/Results/Unique_found_results.txt' \
    #                % (self.folder, self.folder)
    #     os.system(cmd_uniq)
    #
    # def create_venn_subsets(self):
    #     df = pd.read_csv("%s/Results/Summarized_found_results.txt" % self.folder, sep="\t")
    #     start = df["Start"].tolist()
    #     end = df["End"].tolist()
    #     start_end = ["%s" % k + "%s" % j for k, j in zip(start, end)]
    #     dataframe = pd.DataFrame({})
    #     dataframe.insert(0, "coordinates", start_end)
    #     dataframe.to_csv("%s/Results/%s_orf_subset" % (self.folder, self.filetype), sep="\t", index=False)

# def teste():
#     with open("Genome/Results/genome_orf_subset", 'r') as genome, \
#             open("Transcriptome/Results/transcriptome_orf_subset", 'r') as rna:
#         genome_subset = genome.readlines()
#         rna_st = rna.readlines()
#         common = []
#         unique = []
#         for line in rna_st:
#             if line in genome_subset:
#                 common.append(line)
#             elif line not in genome_subset:
#                 unique.append(line)
#         print(common)
#         print(len(common))
#         print(len(unique))
#         print(unique)

# ISTSATPAVHPRYPSRSPAAPARSSSRYWTSRRPPVSETRPPAPHRASSTPAPAARLDSETSAPTIRASGTPLSPA
# teste()
# os.chdir("/home/farminf/Eduardo/smORFs_SMEG/pipeline/")

# merge_databases()
# genome_results = Results("Genome", "genome", "Genome/genome_database.fasta")
# genome_results.rename_files()
# genome_results.write_results()
# genome_results.merge()
# genome_results.add_orf()
# genome_results.add_orf_sequence()
# genome_results.total_spec_count()
# genome_results.genome_coordinates()
# genome_results.check_unique()
# transcriptome_results = Results("Transcriptome", "transcriptome",
# "Transcriptome/transcriptome_database.fasta")
# transcriptome_results.rename_files()
# transcriptome_results.write_results()
# transcriptome_results.merge()
# transcriptome_results.add_orf()
# transcriptome_results.add_orf_sequence()
# transcriptome_results.total_spec_count()
# transcriptome_results.genome_coordinates()
# transcriptome_results.check_unique()
