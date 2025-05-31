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


from ..sequtils import TranscriptomeTranslator, GenomeTranslator, ORFCollection, DatabaseGenerator
from Bio import SeqIO


class Database(object):
    def __init__(self, args):
        self.args = args
        self.genome = args.genome

    def translate(self):
        if self.args.transcriptome:
            rna = TranscriptomeTranslator(sequence="HISAT/transcripts.fasta", form='fasta',
                                                 minsize=int(self.args.minsize), maxsize=int(self.args.maxsize))
            # records = SeqIO.parse("HISAT/transcripts.fasta", 'fasta')
            # for record in records:
            #     print(record.seq)
            rna_orfs = rna.parse_frames(starts=self.args.starts, stops=self.args.stops, seqtype='aa', entry='full')

            db = DatabaseGenerator(name="transcriptome", db_type="sql")
            db.add_orfs(rna_orfs)
            db.to_fasta(filename="transcriptome_ORFs.fasta", identifier='t')
        dna = GenomeTranslator(sequence=self.args.genome, form='fasta', minsize=int(self.args.minsize),
                                      maxsize=int(self.args.maxsize))
        dna_orfs = dna.parse_frames(starts=self.args.starts, stops=self.args.stops,
                                    seqtype='aa')
        genome_db = DatabaseGenerator(name="genome", db_type="sql")
        genome_db.add_orfs(dna_orfs)
        genome_db.to_fasta(filename="genome_ORFs.fasta", identifier='g')
