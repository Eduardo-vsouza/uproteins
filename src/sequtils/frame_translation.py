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


import sys

from Bio import SeqIO
import regex as re

from .conversion.translate import Translator
from .conversion import StrandConverter
from .orflib import ORF, ORFCollection
from .__helpers import FormatError


class FrameTranslator(object):
    def __init__(self, sequence, form, frames, minsize=0, maxsize=0):
        self.sequence = sequence
        self.form = form
        self.frames = frames
        self.minSize = minsize
        self.maxsize = maxsize
        self.orfs = []

    def read_genome(self):
        if self.format == 'fasta':
            records = SeqIO.parse(self.sequence, 'fasta')
            entries = []
            seqs = []
            for record in records:
                entries.append(record.id)
                seqs.append(record.seq)
            return seqs, entries
        elif self.format == 'list':
            return self.sequence
        else:
            raise FormatError

    def parse_frames(self, starts=['ATG'], stops=['TGA', 'TAA', 'TAG'], seqtype='both', **kwargs):
        """ Specify a list of start and stop codons. By default, only 'ATG' is included as a start. TGA, TAA and TAG
         are the default stop codons.'seqtype' accepts 'aa', 'cds', or 'both'. It changes which type of sequence the
         function adds to the 'ORF' object. To save memory, you may exclude one of them. """
        sequences = self.seqsToTranslate
        # defines pattern for regex
        start_pattern = self.__get_pattern(starts, codon_type="start")
        stop_pattern = self.__get_pattern(stops, codon_type="stop")

        if kwargs.get("entry") == "full" and self.form == "fasta":
            entries = self.entries
        else:
            entries = ["" for i in range(len(sequences))]


        # retrieves a list of instances of the ORF class
        forward_orfs = self.__get_cds(sequences, 'forward', start_pattern, stop_pattern, seqtype, entries)
        # adds the orfs to the self.orfs attribute
        self.__add_orfs(forward_orfs)

        # does the same for the complementar strand
        if self.frames == 6:
            # generates a complement to the sequence and reverses it
            reverse_sequences = StrandConverter(sequences)
            rev = reverse_sequences.complement().reverse()

            reverse_orfs = self.__get_cds(rev, 'reverse', start_pattern, stop_pattern, seqtype, entries)
            self.__add_orfs(reverse_orfs)
        # returns a instance of the iterator class 'ORFCollection'
        return ORFCollection().add_orfs(self.orfs)

    def __add_orfs(self, orf_list):
        """ Appends all ORFs in a list to the 'orfs' attribute of this class' instance. """
        for orf in orf_list:
            self.orfs.append(orf)
        return self

    @staticmethod
    def __get_pattern(codons, codon_type=None):
        """ codon_type is either 'start' or 'stop'. 'codons' must refer to a list of start or stop codons."""
        pattern = ""
        for i in range(len(codons)):
            if i != 0:
                if codon_type == "stop":
                    pattern += f"|(?={codons[i]})"
                elif codon_type == "start":
                    pattern += f"|({codons[i]})"
            else:
                if codon_type == "stop":
                    pattern += f"(?={codons[i]})"
                elif codon_type == "start":
                    pattern += f"({codons[i]})"
        return pattern

    def __get_cds(self, sequences, strand, start_pattern, stop_pattern, seqtype, entries):
        """ Returns all possible CDS for the three frames of a nucleotide sequence. """
        orfs = []
        orf_number = 0
        for seq in range(len(sequences)):
            aa = re.finditer('(%s)(...)+?(%s)' % (start_pattern, stop_pattern), "%s" % sequences[seq], overlapped=True)
            for a in aa:
                orf = a.group()
                start = a.start()
                end = a.end()
                if strand == "reverse":
                    start = len(sequences[seq]) - start - 1
                    end = len(sequences[seq]) - end + 1
                if self.__within_limits(orf) and "*" not in orf:
                    orf_number += 1
                    # orf_i = self.__check_seqtype(seqtype=seqtype, start=start+1, end=end, orf_number=orf_number,
                    #                              strand=strand, cds=orf, chromosome=seq+1)
                    orf_i = self.__check_seqtype(seqtype=seqtype, start=start+1, end=end, orf_number=f'{entries[seq]}_{orf_number}',
                                                 strand=strand, cds=orf, chromosome=seq+1)
                    if orf_i is not None:
                        orfs.append(orf_i)

        return orfs

    def __check_seqtype(self, seqtype, **kwargs):
        """ Check which type of sequence should be added to the 'ORF' object. """
        start = kwargs.get("start")
        end = kwargs.get("end")
        orf_number = kwargs.get("orf_number")
        strand = kwargs.get("strand")
        cds = None
        orf = kwargs.get("cds")
        protein = None
        chromosome = None
        transcript = None
        origin = None
        if self.frames == 6:
            na_type = 'chromosome'
            chromosome = kwargs.get('chromosome')
            origin = "Genome"
        elif self.frames == 3:
            na_type = 'transcript'
            strand = ""
            transcript = kwargs.get('chromosome')
            origin = "Transcriptome"
        if seqtype == "cds" or seqtype == "both":
            cds = kwargs.get("cds")
        if seqtype == "aa" or seqtype == "both":
            aa = Translator(orf)
            protein = aa.translate()
        if "*" not in protein:
            orf_i = ORF(name=f'ORF_{orf_number}', start=start, end=end, cds=cds, seq=protein,
                        strand=strand, chromosome=chromosome, transcript=transcript, origin=origin)
            return orf_i

    def __adapt_coordinates(self, start, end, strand):
        """ Fixes start and end region of the ORF in the sequence, related to its current strand. """
        if strand == "forward":
            return start, end
        elif strand == "reverse":
            pass


    def __within_limits(self, orf):
        """ Checks if the predicted ORF is within the size constraints specified by 'minsize' and 'maxsize'. """
        if self.minSize <= len(orf) <= self.maxsize:
            return True
        else:
            return False


class GenomeTranslator(FrameTranslator):
    def __init__(self, sequence=None, form='fasta', minsize=30, maxsize=300):
        """ Sequence must be either a list or the path to a fasta file. The format must be specified by 'form', which is
        set to 'fasta' by default. If it is a list, each element must be a chromosome and thus contain all its
        nucleotide sequences. minsize and maxsize refer to the ORF length in nucleotides, and are set to 30 and 300
         by default, respectively, in order to detect only small ORFs. """
        self.frames = 6
        self.format = form
        self.sequence = sequence
        super().__init__(self.sequence, form, self.frames, minsize=minsize, maxsize=maxsize)
        self.seqsToTranslate, self.entries = self.read_genome()


class TranscriptomeTranslator(FrameTranslator):
    def __init__(self, sequence=None, form='fasta', minsize=30, maxsize=300):
        """ Sequence must be either a list or the path to a fasta file. The format must be specified by 'form', which is
        set to 'fasta' by default. If it is a list, each element must be a chromosome and thus contain all its
        nucleotide sequences. minsize and maxsize refer to the ORF length in nucleotides, and are set to 30 and 300
         by default, respectively, in order to detect only small ORFs. """
        self.frames = 3
        self.format = form
        self.sequence = sequence
        super().__init__(self.sequence, form, self.frames, minsize=minsize, maxsize=maxsize)
        self.seqsToTranslate, self.entries = self.read_genome()



# seq = ["ATGTGCTGATGCATGATGCTATAAGGCCAT", 'ATGTTTTAGTTT']
#
# rna = TranscriptomeTranslator(sequence='mtb_genome.fasta', form='fasta', minsize=30, maxsize=300)
# orfs = rna.parse_frames(starts=['ATG', 'TGT'], stops=['TGA', 'TAA', 'TAG'], seqtype='aa')
# #
# # dna = GenomeTranslator(sequence='mtb_genome.fasta', form='fasta', minsize=30, maxsize=300)
# # d_orfs = dna.parse_frames(starts=['ATG', 'TGT'], stops=['TGA', 'TAA', 'TAG'], seqtype='aa')
#
# for orf in orfs:
#     print(orf.name, orf.seq, orf.origin, orf.transcript)