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



class StrandConverter(object):
    def __init__(self, sequence):
        self.sequence = sequence
        self.complementary = []

    def complement(self):
        nuc_table = {
            'A': 'T',  'T': 'A', 'G': 'C', 'C': 'G'
        }
        for chromosome in self.sequence:
            nuc_seq = ""
            for nuc in chromosome:
                nuc_seq += nuc_table[nuc]
            if "*" not in nuc_seq:
                self.complementary.append(nuc_seq)
        return self

    def reverse(self):
        rev = [i[::-1] for i in self.complementary]
        return rev
