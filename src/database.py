from pyteogenomics import TranscriptomeTranslator, GenomeTranslator, ORFCollection, DatabaseGenerator
from Bio import SeqIO

class Database(object):
    def __init__(self, args):
        self.args = args
        self.genome = args.genome

    def translate(self):
        if self.args.Transcriptome is not None:
            rna = TranscriptomeTranslator(sequence="HISAT/transcripts.fasta", form='fasta',
                                                 minsize=int(self.args.minsize), maxsize=int(self.args.maxsize))
            # records = SeqIO.parse("HISAT/transcripts.fasta", 'fasta')
            # for record in records:
            #     print(record.seq)
            rna_orfs = rna.parse_frames(starts=self.args.starts.split(","), stops=self.args.stops.split(","), seqtype='aa', entry='full')

            db = DatabaseGenerator(name="transcriptome", db_type="sql")
            db.add_orfs(rna_orfs)
            db.to_fasta(filename="transcriptome_ORFs.fasta", identifier='t')
        dna = GenomeTranslator(sequence=self.args.genome, form='fasta', minsize=int(self.args.minsize),
                                      maxsize=int(self.args.maxsize))
        dna_orfs = dna.parse_frames(starts=self.args.starts.split(","), stops=self.args.stops.split(","),
                                    seqtype='aa')
        genome_db = DatabaseGenerator(name="genome", db_type="sql")
        genome_db.add_orfs(dna_orfs)
        genome_db.to_fasta(filename="genome_ORFs.fasta", identifier='g')
