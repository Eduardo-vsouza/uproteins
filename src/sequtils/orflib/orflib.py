from Bio import SeqIO


class ORF(object):
    def __init__(self, name=None, seq=None, start=None, end=None, cds=None, strand=None, chromosome=None,
                 transcript=None, origin=None, appearances=1, experiment=None, protein_sequence=None):
        self.name = name
        self.seq = seq
        self.start = start
        self.end = end
        self.cds = cds
        self.strand = strand
        self.chromosome = chromosome
        self.transcript = transcript
        self.origin = origin
        self.start_codon = None
        self.shineDalgarno = None
        self.freeEnergy = None
        self.proteinSequence = protein_sequence

        self.MSPeptides = []

        # for spectral counting
        self.appearances = appearances
        self.nsaf = None
        self.experiment = experiment
        self.normalizedSpec = None

        # after using msprocess
        self.peptides = {}
        self.closestToStart = self.find_ms_peptides()

    def set_coordinates(self, genome):
        return self

    def __len__(self):
        length = len(self.seq)
        return length

    def find_ms_peptides(self):

        # if len(self.MSPeptides) >= 1:
        if self.strand == 'reverse':
            closest_to_start = int(self.end)
        else:
            closest_to_start = int(self.end)
        for peptide in self.MSPeptides:
            seq_start = self.seq.find(peptide) * 3
            if self.strand == 'reverse':
                genome_start = int(self.start) - seq_start

                if genome_start > closest_to_start:
                    closest_to_start = genome_start
            else:
                genome_start = int(self.start) + seq_start
                if genome_start <= closest_to_start:
                    closest_to_start = genome_start
        self.closestToStart = closest_to_start
        # if len(self.MSPeptides) > 1:
        #     print(self.end, self.start, closest_to_start, self.MSPeptides, self.seq, self.strand)
        return closest_to_start



class ORFCollection(object):
    def __init__(self):
        self.entries = None
        self.seqs = None
        self.starts = None
        self.ends = None
        self.coordinates = None
        self.orfs = []
        self.priority = None

    def __iter__(self):
        i = 0
        while True:
            yield self.orfs[i]
            i += 1
            if i >= self.__len__():
                break

    def set_priority(self, orf):
        self.priority = orf

    def read_fasta(self, file):
        """ Adds all entries and sequences in a fasta file to this class instance. """
        records = SeqIO.parse(file, 'fasta')
        for record in records:
            self.orfs.append(ORF(name=record.name, seq=record.seq))
        return self

    def add_orfs(self, orfs):
        """ Add ORFs translated with GenomeTranslator or TranscriptomeTranslator method. """
        for orf in orfs:
            self.orfs.append(orf)
        return self

    def add_orf(self, orf):
        """ Add a single ORF to collection."""
        self.orfs.append(orf)

    def __len__(self):
        return len(self.orfs)

    def to_fasta(self, filename="orfs.fasta"):
        """ Generates a fasta file containing all ORFs in this instance. """
        to_write = [f">{orf.name}_{orf.start}-{orf.end}_{orf.strand}\n{orf.seq}\n" for orf in self.orfs]
        with open(filename, 'w') as fa:
            fa.writelines(to_write)


# orfs = ORFs()
# orfs.read_fasta('test.fasta.txt')
# for orf in orfs:
#     print(orf.name)
