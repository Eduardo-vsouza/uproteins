import pandas as pd
from Bio import SeqIO

from . import ORF, ORFCollection


class AltCodons(object):
    def __init__(self, file, genome):
        self.df = pd.read_csv(file, sep='\t')
        self.coordinates = self.df["Genome Coordinates"].tolist()
        self.names = self.df["db entry"].tolist()
        self.alternatives = self.__fetch_orfs()
        self.__genome_records = SeqIO.parse(genome, 'fasta')
        self.genome_seq = [str(record.seq) for record in self.__genome_records]

    def __split_coords(self, i):
        splat = self.coordinates.split("-")
        if 'reverse' in self.names[i]:
            start = splat[1]
            end = splat[0]
            strand = 'reverse'
        else:
            start = splat[0]
            end = splat[1]
            strand = 'forward'
        return start, end, strand

    def __fetch_orfs(self):
        """
        :returns a dictionary containing all ORFs with alternative START codons for a given STOP codon.
        """
        alternatives = {}
        for i in range(len(self.names)):
            start, end, strand = self.__split_coords(i)
            orf = ORF(name=self.names[i], start=start, end=end, strand=strand)
            orf = self.__fetch_codons(orf)
            if end not in alternatives:
                alternatives[end] = ORFCollection()
            alternatives[end].add_orf(orf)
        return alternatives

    def sort_starts_by_coordinates(self):
        """ Sorts the alternative ORFs inside self.alternatives by their start codons. """
        for alts in self.alternatives:
            starts = []
            for alt in self.alternatives[alts]:
                if len(starts) == 0:
                    starts.append(alt)
                else:
                    if alt.strand == 'forward':
                        if starts[0].start > alt.start:
                            starts.insert(0, alt)
                    else:
                        if starts[0].start < alt.start:
                            starts.insert(0, alt)
            self.alternatives[alts] = ORFCollection().add_orf(starts)

    def __fetch_codons(self, orf):
        """ :returns the nucleotide sequence of the start codon for a given ORF. """
        s_codon = self.genome_seq[0][orf.start: orf.start+3]
        orf.start_codon = s_codon
        return orf
