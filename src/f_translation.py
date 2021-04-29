import pandas as pd
from Bio import SeqIO
from difflib import SequenceMatcher


# def fix_weird_characters():


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
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
            'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
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
    def __init__(self, genome, orf, coordinates):
        self.genome = genome
        self.orf = orf
        self.coordinates = coordinates

    def find_frame(self, locus1, first_tr, second_tr, third_tr, first_rv_tr_rv, second_rv_tr_rv, third_rv_tr_rv):

        orf_seq = self.orf

        coordenadas = self.coordinates
        coord = coordenadas.split(" - ")
        print(coord)
        start = coord[0]
        end = coord[1]
        startc = ""
        endc = ""
        if start > end:
            startc = end
            endc = start
        elif start < end:
            startc = start
            endc = end

        if len(start) < len(end):
            startc = start
            endc = end



        """ new try, now translating after finding the frame """
        # print(first_tr)
        first_translated = first_tr[int(int(startc)/3-15):int(int(endc)/3+15)]

        second_translated = second_tr[int(int(startc)/3-15):int(int(endc)/3+15)]
        third_translated = third_tr[int(int(startc)/3-15):int(int(endc)/3+15)]
        # print(first_translated)
        # print(second_translated)
        # print(third_translated)


        first_rv_locus = first_rv_tr_rv[int(int(startc)/3-15):int(int(endc)/3+15)][::-1]
        second_rv_locus = second_rv_tr_rv[int(int(startc)/3-15):int(int(endc)/3+15)][::-1]
        third_rv_locus = third_rv_tr_rv[int(int(startc)/3-15):int(int(endc)/3+15)][::-1]
        reading_frame = ""
        # print("\nORF", orf_seq, "\n")
        translations_dict = {first_translated: "RF +1", second_translated: "RF +2", third_translated: "RF +3",
                             first_rv_locus: "RF -1", second_rv_locus: "RF -2", third_rv_locus: "RF -3"}
        translations = [first_translated, second_translated, third_translated, first_rv_locus, second_rv_locus,
                        third_rv_locus]
        # for trans in translations:
        #     diff_seq = SequenceMatcher(None, orf_seq[1:], trans).ratio()
        #     # print(diff_seq*100)
        #     print(diff_seq)
        #     print(orf_seq[1:])
        #     print(trans)
        #     break
        #     if (diff_seq*100) >= 90:
        #         reading_frame = translations_dict[trans]
        #     else:
        #         reading_frame = "RF"
        if orf_seq[1:] in first_translated:
            reading_frame = "RF +1"
        elif orf_seq[1:] in second_translated:
            reading_frame = "RF +2"
        elif orf_seq[1:] in third_translated:
            reading_frame = "RF +3"
        elif orf_seq[1:] in first_rv_locus:
            reading_frame = "RF -1"
        elif orf_seq[1:] in second_rv_locus:
            reading_frame = "RF -2"
        elif orf_seq[1:] in third_rv_locus:
            reading_frame = "RF -3"
        else:
            reading_frame = "RF"
        return reading_frame



def find_orfs_rfs(file, genome, filetype):
    with open(genome, 'r') as genome_file:

        lines = genome_file.readlines()[1:]
        frame1 = []
        for line in lines:
            line = line.rstrip()
            frame1.append(line)

        locus1 = str(''.join(map(str, frame1)))
        first_frame = Translator(locus1)
        first_tr = first_frame.translate()
        second_frame = Translator(locus1[1:])
        second_tr = second_frame.translate()
        third_frame = Translator(locus1[2:])
        third_tr = third_frame.translate()

        fff_locus = str(locus1)[::-1]
        fff_third_reading_frame = fff_locus[6:]
        fff_second_reading_frame = fff_locus[4:]
        fff_first_reading_frame = fff_locus[2:]
        fff_1 = Translator(fff_first_reading_frame)
        fff_2 = Translator(fff_second_reading_frame)
        fff_3 = Translator(fff_third_reading_frame)
        first_rv_genome = fff_1.complement()
        second_rv_genome = fff_2.complement()
        third_rv_genome = fff_3.complement()
        rv1 = Translator(first_rv_genome)
        rv2 = Translator(second_rv_genome)
        rv3 = Translator(third_rv_genome)
        first_rv_tr = rv1.translate()
        second_rv_tr = rv2.translate()
        third_rv_tr = rv3.translate()
        first_rv_tr_rv = first_rv_tr[::-1]
        second_rv_tr_rv = second_rv_tr[::-1]
        third_rv_tr_rv = third_rv_tr[::-1]
        df = pd.read_csv(file, sep="\t")
        seqs = df['ORF Sequence'].tolist()
        coordenadas = df['Genome Coordinates'].tolist()
        rfs = []
        for i in range(len(seqs)):
            sequence = ReadingFrame(genome, seqs[i], coordenadas[i])
            try:
                rf = sequence.find_frame(locus1, first_tr, second_tr, third_tr, first_rv_tr_rv, second_rv_tr_rv,
                                         third_rv_tr_rv)

                rfs.append(rf)
            except:
                rfs.append("RF")

        df.insert(4, "Reading Frame", rfs)
        # print(df)
        df.to_csv("%s_results_with_rfs.ods" % filetype, sep="\t", header=True, index=False)


def find_proteome_rfs(genome, gb_file, folder):
    with open(genome, 'r') as genome_file:

        lines = genome_file.readlines()[1:]
        frame1 = []
        for line in lines:
            line = line.rstrip()
            frame1.append(line)

        locus1 = str(''.join(map(str, frame1)))
        first_frame = Translator(locus1)
        first_tr = first_frame.translate()
        second_frame = Translator(locus1[1:])
        second_tr = second_frame.translate()
        third_frame = Translator(locus1[2:])
        third_tr = third_frame.translate()

        fff_locus = str(locus1)[::-1]
        fff_third_reading_frame = fff_locus[6:]
        fff_second_reading_frame = fff_locus[4:]
        fff_first_reading_frame = fff_locus[2:]
        fff_1 = Translator(fff_first_reading_frame)
        fff_2 = Translator(fff_second_reading_frame)
        fff_3 = Translator(fff_third_reading_frame)
        first_rv_genome = fff_1.complement()
        second_rv_genome = fff_2.complement()
        third_rv_genome = fff_3.complement()
        rv1 = Translator(first_rv_genome)
        rv2 = Translator(second_rv_genome)
        rv3 = Translator(third_rv_genome)
        first_rv_tr = rv1.translate()
        second_rv_tr = rv2.translate()
        third_rv_tr = rv3.translate()
        first_rv_tr_rv = first_rv_tr[::-1]
        second_rv_tr_rv = second_rv_tr[::-1]
        third_rv_tr_rv = third_rv_tr[::-1]


        seqs = []
        loci = []
        names = []
        for seq_record in SeqIO.parse(gb_file, "genbank"):
            for feature in seq_record.features:
                if feature.qualifiers.get("translation") is not None:
                    seqs.append(str(feature.qualifiers.get("translation"))[2:-2])
                    loci.append(str(feature.location))
                    names.append(str(feature.qualifiers.get("locus_tag")))

        starts = []
        ends = []
        for i in range(len(loci)):
            strand = loci[i].split("(")
            if strand[1] == "+)":
                st = loci[i].split(":")
                starts.append(int(st[0][1:]))
                ends.append(int((st[1][:-4])))
            elif strand[1] == "-)":
                st = loci[i].split(":")
                starts.append(int(st[1][:-4]))
                ends.append(int((st[0][1:])))
        coords_final = []
        for i in range(len(starts)):
            coords_final.append("%i - %i" %(starts[i], ends[i]))

        rfs = []
        for i in range(len(seqs)):
            sequence = ReadingFrame(genome, seqs[i], coords_final[i])
            rf = sequence.find_frame(locus1, first_tr, second_tr, third_tr, first_rv_tr_rv, second_rv_tr_rv,
                                     third_rv_tr_rv)
            rfs.append(rf)
        df = pd.DataFrame()
        df.insert(0, "Gene", names)
        df.insert(1, "Sequence", seqs)
        df.insert(2, "Start", starts)
        df.insert(3, "End", ends)
        df.insert(4, "Reading Frame", rfs)
        df = df[df["Reading Frame"] != "RF"]
        df.to_csv("RefSeq_Reading_Frames.ods", sep="\t", header=True, index=False)
