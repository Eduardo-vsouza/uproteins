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


from pyteogenomics import StringTieGFF, GenomeCoordinates, RefSeqGFF, GenomeCoordinatesRNA


class PostMS(object):
    def __init__(self, args):
        self.outdir = args.outdir
        self.transcriptomeDB = f'{self.outdir}/transcriptome_database.fasta'
        self.genomeDB = f'{self.outdir}/genome_database.fasta'
        # self.

# str_gff = StringTieGFF(gff='assembled.gtf')
# orf_dict = str_gff.get_dict()
# ref_gff = RefSeqGFF(gff='msmeg.gff3')
# ref_dict = ref_gff.get_dict()
# coordinates = GenomeCoordinatesRNA('2611_ncbi/transcriptome_converted_psm.txt', orf_dict, ref_dict)
# coordinates.get_proteins().save_table(output='2611_ncbi/transcriptome_psm_coords')