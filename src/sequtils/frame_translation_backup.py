import sys

from Bio import SeqIO

from conversion.translate import Translator
from conversion import StrandConverter
from orflib import ORF, ORFs


class FrameTranslator(object):
    def __init__(self, sequence, form, frames, minsize=0, maxsize=0):
        self.sequence = sequence
        self.form = form
        self.frames = frames
        self.minSize = minsize
        self.maxSize = maxsize
        self.orfs = []

    def read_genome(self):
        if self.format == 'fasta':
            records = SeqIO.parse(self.sequence, 'fasta')
            return [record.seq for record in records]
        elif self.format == 'list':
            return self.sequence
        else:
            print("Invalid format. Please specify either 'fasta' or 'list'. ")

    def parse_frames(self, starts=['ATG'], stops=['TGA', 'TAA', 'TAG'], seqtype='both'):
        """ Specify a list of start and stop codons. By default, only 'ATG' is included as a start. TGA, TAA and TAG
         are the default stop codons.'seqtype' accepts 'aa', 'cds', or 'both'. It changes which type of sequence the
         function adds to the 'ORF' object. To save memory, you may exclude one of them. """
        sequences = self.seqsToTranslate
        forward_orfs = self.__get_cds(sequences, starts, stops, 'forward', seqtype)
        self.__add_orfs(forward_orfs)
        if self.frames == 6:
            reverse_sequences = StrandConverter(sequences)
            rev = reverse_sequences.complement().reverse()
            reverse_orfs = self.__get_cds(rev, starts, stops, 'reverse', seqtype)
            self.__add_orfs(reverse_orfs)
        return ORFs().add_orfs(self.orfs)

    def __add_orfs(self, orf_list):
        """ Appends all ORFs in a list to the 'orfs' attribute of this class' instance. """
        for orf in orf_list:
            self.orfs.append(orf)
        return self

    def __get_cds(self, sequences, starts, stops, strand, seqtype):
        """ Returns all possible CDS for the three frames of a nucleotide sequence. """
        orfs = []
        orf_number = 0
        for seq in range(len(sequences)):
            for i in range(len(sequences[seq])):
                if sequences[seq][i:i + 3] in starts:
                    orf = ""
                    orf += sequences[seq][i:i + 3]
                    j = 3
                    cod = sequences[seq][i + j:i + 3 + j]
                    orf += cod
                    while cod not in stops and len(cod) == 3:
                        if j > len(sequences[seq]):
                            break
                        j += 3
                        cod = sequences[seq][i + j:i + 3 + j]
                        orf += cod
                    if orf[-3:] not in stops:
                        orf = ""
                    for k in range(len(stops)):
                        if orf.endswith(stops[k]):
                            orf = orf[:-3]

                    if self.__within_limits(orf):
                        orf_number += 1
                        start = sequences[seq].find(orf) + 1
                        end = start + len(orf) - 1
                        orf_i = self.__check_seqtype(seqtype=seqtype, start=start, end=end, orf_number=orf_number,
                                                     strand=strand, cds=orf, chromosome=seq+1)
                        orfs.append(orf_i)

        return orfs

    def __check_seqtype(self, seqtype, **kwargs):
        """ Check which type of sequence should be added to the 'ORF' object. """
        start = kwargs.get("start")
        end = kwargs.get("end")
        orf_number = kwargs.get("orf_number")
        strand = f'_{kwargs.get("strand")}'
        cds = None
        orf = kwargs.get("cds")
        protein = None
        chromosome = kwargs.get('chromosome')
        if self.frames == 6:
            na_type = 'chromosome'
        elif self.frames == 3:
            na_type = 'transcript'
            strand = ""
        if seqtype == "cds" or seqtype == "both":
            cds = kwargs.get("cds")
        if seqtype == "aa" or seqtype == "both":
            aa = Translator(orf)
            protein = aa.translate()
        orf_i = ORF(name=f'ORF{orf_number}_{na_type}{chromosome}{strand}', start=start, end=end, cds=cds, seq=protein)
        return orf_i

    def __adapt_coordinates(self, start, end, strand):
        """ Fixes start and end region of the ORF in the sequence, related to its current strand. """
        if strand == "forward":
            return start, end
        elif strand == "reverse":
            pass


    def __within_limits(self, orf):
        """ Checks if the predicted ORF is within the size constraints specified by 'minsize' and 'maxsize'. """
        if self.minSize <= len(orf) <= self.maxSize:
            return True
        else:
            return False


class GenomeTranslator(FrameTranslator):
    def __init__(self, sequence=None, form='fasta', minsize=300, maxsize=30000):
        """ Sequence must be either a list or the path to a fasta file. The format must be specified by 'form', which is
        set to 'fasta' by default. If it is a list, each element must be a chromosome and thus contain all its
        nucleotide sequences. minsize and maxsize refer to the ORF length in nucleotides, and are set to 30 and 300
         by default, respectively, in order to detect only small ORFs. """
        self.frames = 6
        self.format = form
        self.sequence = sequence
        super().__init__(self.sequence, form, self.frames, minsize=minsize, maxsize=maxsize)
        self.seqsToTranslate = self.read_genome()


class TranscriptomeTranslator(FrameTranslator):
    def __init__(self, sequence=None, form='fasta', minsize=300, maxsize=30000):
        """ Sequence must be either a list or the path to a fasta file. The format must be specified by 'form', which is
        set to 'fasta' by default. If it is a list, each element must be a chromosome and thus contain all its
        nucleotide sequences. minsize and maxsize refer to the ORF length in nucleotides, and are set to 30 and 300
         by default, respectively, in order to detect only small ORFs. """
        self.frames = 3
        self.format = form
        self.sequence = sequence
        super().__init__(self.sequence, form, self.frames, minsize=minsize, maxsize=maxsize)
        self.seqsToTranslate = self.read_genome()



seq = ["ATGTGCTGATGCATGATGCTATAAGGCCAT", 'ATGTTTTAGTTT']

rna = TranscriptomeTranslator(sequence='mtb_genome.fasta', form='fasta', minsize=300, maxsize=4000000)
orfs = rna.parse_frames(starts=['ATG', 'TGT'], stops=['TGA', 'TAA', 'TAG'], seqtype='aa')

dna = GenomeTranslator(sequence=seq, form='list', minsize=300)
d_orfs = dna.parse_frames(starts=['ATG', 'TGT'], stops=['TGA', 'TAA', 'TAG'])


