from Bio import SeqIO


class Translator(object):
    def __init__(self, genome):
        self.genome = genome

    def translate(self):
        table = {
            'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
            'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
        }
        protein = ""
        for i in range(0, len(self.genome), 3):
            codon = self.genome[i:i + 3]
            if len(codon) == 3:
                protein += table[codon]
        return protein

    def complement(self):
        nuc_table = {
            'A': 'T',  'T': 'A', 'G': 'C', 'C': 'G'
        }
        nuc_seq = ""
        for i in self.genome:
            nuc_seq += nuc_table[i]
        return nuc_seq


class ReadingFrame(object):
    def __init__(self, genome):
        self.genome = genome
        # self.orf = orf
        # self.coordinates = coordinates

    def find_frame(self):
        records = SeqIO.read(self.genome, 'fasta')

        proteins = []

        seq_len = len(records.seq)
        answer = []
        for strand, nuc in [(+1, records.seq), (-1, records.seq.reverse_complement())]:
            for frame in range(3):
                protein = Translator(str(nuc[frame:]))
                trans = str(protein.translate())
                print(trans)
                trans_len = len(trans)
                aa_start = 0
                aa_end = 0
                while aa_start < trans_len:
                    aa_end = trans.find("*", aa_start)
                    if aa_end == -1:
                        aa_end = trans_len
                    # 10 refers to the min protein length
                    if aa_end - aa_start >= 10:
                        if strand == 1:
                            start = frame + aa_start * 3
                            end = min(seq_len, frame + aa_end * 3)
                        else:
                            start = seq_len - frame - aa_end * 3
                            end = seq_len - frame - aa_start * 3
                        answer.append("%i\t%i\t%i\t%s\t%i\n" % (start, end, strand, trans[aa_start:aa_end], frame))
                    aa_start = aa_end + 1

        output = []


                    # length = 3 * ((len(records)-frame) // 3)
                    # for pro in nuc[frame:frame+length].translate(11).split("*"):
                    #     proteins.append("%s\tstrand %i\t frame %i\n" % (pro, strand, frame))

        with open("teste_Reading_frames.ods", 'w') as out:
            out.writelines(answer)
        #
        #
        # orf_seq = self.orf
        #
        # coordenadas = self.coordinates
        # coord = coordenadas.split(" - ")
        # start = coord[0]
        # end = coord[1]
        # startc = ""
        # endc = ""
        # if start > end:
        #     startc = end
        #     endc = start
        # elif start < end:
        #     startc = start
        #     endc = end
        #
        # frame1 = []
        # for line in lines:
        #     line = line.rstrip()
        #     frame1.append(line)



            """ new try, now translating after finding the frame """
            # first_frame = Translator(locus1)
            # first_tr = first_frame.translate()
            # second_frame = Translator(locus1[1:])
            # second_tr = second_frame.translate()
            # third_frame = Translator(locus1[2:])
            # third_tr = third_frame.translate()
            # first_translated = first_tr[int(int(startc)/3-15):int(int(endc)/3+15)]
            # second_translated = second_tr[int(int(startc)/3-15):int(int(endc)/3+15)]
            # third_translated = third_tr[int(int(startc)/3-15):int(int(endc)/3+15)]
            # # print(first_translated)
            # # print(second_translated)
            # # print(third_translated)
            #
            # fff_locus = str(locus1)[::-1]
            # fff_third_reading_frame = fff_locus[6:]
            # fff_second_reading_frame = fff_locus[4:]
            # fff_first_reading_frame = fff_locus[2:]
            # fff_1 = Translator(fff_first_reading_frame)
            # fff_2 = Translator(fff_second_reading_frame)
            # fff_3 = Translator(fff_third_reading_frame)
            # first_rv_genome = fff_1.complement()
            # second_rv_genome = fff_2.complement()
            # third_rv_genome = fff_3.complement()
            # rv1 = Translator(first_rv_genome)
            # rv2 = Translator(second_rv_genome)
            # rv3 = Translator(third_rv_genome)
            # first_rv_tr = rv1.translate()
            # second_rv_tr = rv2.translate()
            # third_rv_tr = rv3.translate()
            # first_rv_tr_rv = first_rv_tr[::-1]
            # first_rv_locus = first_rv_tr_rv[int(int(startc)/3-15):int(int(endc)/3+15)][::-1]
            # second_rv_tr_rv = second_rv_tr[::-1]
            # second_rv_locus = second_rv_tr_rv[int(int(startc)/3-15):int(int(endc)/3+15)][::-1]
            # third_rv_tr_rv = third_rv_tr[::-1]
            # third_rv_locus = third_rv_tr_rv[int(int(startc)/3-15):int(int(endc)/3+15)][::-1]
            # # print(first_rv_locus)
            # # print(second_rv_locus)
            # # print(third_rv_locus)
            #
            # reading_frame = ""
            # # print("\nORF", orf_seq, "\n")
            # if orf_seq in first_translated:
            #     reading_frame = "RF +1"
            # elif orf_seq in second_translated:
            #     reading_frame = "RF +2"
            # elif orf_seq in third_translated:
            #     reading_frame = "RF +3"
            # elif orf_seq in first_rv_locus:
            #     reading_frame = "RF -1"
            # elif orf_seq in second_rv_locus:
            #     reading_frame = "RF -2"
            # elif orf_seq in third_rv_locus:
            #     reading_frame = "RF -3"
            # return reading_frame

# # +3
# sequence = ReadingFrame("MSMEG_genome.fasta", "MPRTCWRSIRTMPRCAGSPTVIWSPWPAVSGRPRCTRRPPIGCRPASCTPRSTIP",
#                         "187458 - 187622")
#
# # +1
# sequence2 = ReadingFrame("MSMEG_genome.fasta", "VRRCAKAPTREMNRSKNTVGSTSMPISALGRRPQCPRQSATVIPRVRVPVSSTSTVVWRHRTNSGSNATAAASGSCRCAR", "1711525 - 1711764")
#
# # -2
# sequence3 = ReadingFrame("MSMEG_genome.fasta", "VAPPPGWIPSFWKRGSNSACGPAMRTSAASARLSPAPTAAPFTAAMVGSVQWATARKPS", "4799831 - 4799655")
#
# # -1
# sequence4 = ReadingFrame("MSMEG_genome.fasta", "VRRPGARSPRSTTCGIRRSRDASACCPTSRTAWA", "3361927 - 3361826")
#
# # -3
# sequence5 = ReadingFrame("MSMEG_genome.fasta", "LRARIPGFREGGRLDRTADSWWRATPLRLERMLQRALP", "2170962 - 2170849")
# rf = sequence.find_frame()
# rf2 = sequence2.find_frame()
# rf3 = sequence3.find_frame()
# rf4 = sequence4.find_frame()
# rf5 = sequence5.find_frame()
# print(rf)
# print(rf2)
# print(rf3)
# print(rf4)
# print(rf5)
#
# sequence = ReadingFrame("MSMEG_genome.fasta")
# sequence.find_frame()