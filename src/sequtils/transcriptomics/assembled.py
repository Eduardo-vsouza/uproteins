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


from Bio import SeqIO


class TranscriptExtractor(object):
    def __init__(self, assembly):
        self.assembly = SeqIO.parse(assembly, 'fasta')

    def get_transcripts(self):
        rnas = {}
        for record in self.assembly:
            gene = str(record.id)
            # if '0001' in gene:
            #     print('WORKING', gene)
            if 'gene' in gene:
                gene = gene[5:]
            elif 'uproteins' in gene:
                splat = gene.split(".")
                gene = f'{splat[0]}.{splat[1]}'
            elif 'rna' in gene:
                gene = gene[4:]
            if gene not in rnas:
                rnas[gene] = str(record.seq)
        return rnas
