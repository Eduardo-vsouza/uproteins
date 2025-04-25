import os

from tqdm import tqdm
from Bio import SeqIO

from .f_translation_biopy import Translator


# class Translator(object):
#     def __init__(self, genome):
#         self.genome = genome
#
#     def translate(self):
#         table = {
#             'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
#             'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
#             'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
#             'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
#             'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
#             'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
#             'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
#             'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
#             'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
#             'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
#             'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
#             'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
#             'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
#             'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
#             'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
#             'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
#         }
#         protein = ""
#         for i in range(0, len(self.genome), 3):
#             codon = self.genome[i:i + 3]
#             if len(codon) == 3:
#                 protein += table[codon]
#         return protein
#
#     def complement(self):
#         nuc_table = {
#             'A': 'T',  'T': 'A', 'G': 'C', 'C': 'G'
#         }
#         nuc_seq = ""
#         for i in self.genome:
#             nuc_seq += nuc_table[i]
#         return nuc_seq


class GenomeReader(object):
    def __init__(self, fasta, starts, ends, filetype, args):
        """ start and end coordinates must be separated by a comma, nothing else. """
        self.fasta = fasta
        self.orfs = []
        self.starts = starts
        self.stops = ends
        self.sequence = []
        self.entries = []
        self.c_entries = []
        self.c_seqs = []
        self.filetype = filetype
        self.args = args

    def get_sequence(self):
        sequence = []
        entries = []
        comp_seqs = []
        comp_entries = []
        records = SeqIO.parse(self.fasta, 'fasta')
        for record in records:
            entries.append(str(record.description))
            sequence.append(str(record.seq))
            complement = Translator(str(record.seq))
            complement_seq = complement.complement()
            comp_seqs.append(complement_seq[::-1])
            comp_entries.append(str(record.description))

        self.sequence = sequence
        self.entries = entries
        self.c_entries = comp_entries
        self.c_seqs = comp_seqs

    def three_frame_tr(self, entries, sequences, strand):
        self.orfs = []
        for seq in tqdm(range(len(sequences))):
            orf_number = 0
            for i in range(len(sequences[seq])):
                if sequences[seq][i:i + 3] in self.starts:
                    orf = ""
                    orf += sequences[seq][i:i + 3]
                    j = 3
                    cod = sequences[seq][i + j:i + 3 + j]
                    orf += cod
                    # if len(self.orfs) == 2000 or i == len(sequences[seq]):
                    while cod not in self.stops and len(cod) == 3:
                        # print(cod)
                        if j > len(sequences[seq]):
                            break
                        j += 3
                        cod = sequences[seq][i + j:i + 3 + j]
                        orf += cod
                        # if j > len(sequences[seq]) and cod not in stops:
                        #     orf = ""
                        # elif j > len(sequences[seq]) and cod in stops:
                        #     break
                    # if orf[len(orf) - 3:] not in stops:
                    #     orf = ""
                    #
                    if orf[-3:] not in self.stops:
                        orf = ""
                    for k in range(len(self.stops)):
                        if orf.endswith(self.stops[k]):
                            orf = orf[:-3]

                    if len(orf) >= self.args.minsize and len(orf) <= self.args.maxsize:
                        orf_number += 1
                        if strand == "+":
                            # start = i+1
                            start = sequences[seq].find(orf) + 1
                            end = start + len(orf) - 1
                            if (">%s_ORF_%s_[%s_-_%s]\n%s\n" % (entries[seq], orf_number, start, end, orf)) not in self.orfs:
                                self.orfs.append(">%s_ORF_%s_[%s_-_%s]\n%s\n" % (entries[seq], orf_number, start, end, orf))
                        elif strand == "-":
                            # print(orf)
                            # print(sequences[seq].find(orf))
                            # end = i-4
                            start = len(sequences[seq]) - i
                            end = start - len(orf) + 1
                            # end = sequences[seq].find(orf) + 1
                            # start = end + len(orf)
                            if (">%s_ORF_%s_[%s_-_%s]\n%s\n" % (entries[seq], orf_number, start, end, orf)) not in self.orfs:
                                self.orfs.append(">%s_ORF_%s_REVERSE_[%s_-_%s]\n%s\n" % (entries[seq], orf_number, start, end, orf))
                        else:
                            print("Incorrect strand")
                            break
        if not os.path.exists("%s_CDS.fasta" % self.filetype):
            with open("%s_CDS.fasta" % self.filetype, 'w') as file:
                file.writelines(self.orfs)
        elif os.path.exists("%s_CDS.fasta" % self.filetype):
            with open("%s_CDS.fasta" % self.filetype, 'a') as file:
                file.writelines(self.orfs)
        # if os.path.exists("testes_tr/CDS.fasta"):
        #     with open("testes_tr/CDS.fasta", 'a') as file:
        #         file.writelines(self.orfs)

    def translate(self):
        entries = []
        replaced_entries = []
        nuc_seqs = []
        aa_seqs = []
        records = SeqIO.parse("%s_CDS.fasta" % self.filetype, 'fasta')

        for record in records:
            entries.append(str(record.description))
            nuc_seqs.append(str(record.seq))
        for i in range(len(entries)):
            entry = entries[i].replace(" ", "_")
            replaced_entries.append(entry)
        for i in range(len(nuc_seqs)):
            orf = Translator(nuc_seqs[i])
            translated = orf.translate()
            m_fixed = "M"+translated[1:]
            aa_seqs.append(">" + replaced_entries[i] +"\n" + m_fixed + "\n")

        with open("%s_ORFs.fasta" % self.filetype, 'w') as out:
            out.writelines(aa_seqs)




