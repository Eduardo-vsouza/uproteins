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


class ProgressPanel(object):
    def __init__(self, mode):
        self.mode = mode

    def __create_panels(self, step):
        up = '\x1B[1A'
        CLR = "\x1B[0K"
        # print("\n\n")  # set up blank lines so cursor moves work
        assembly = f'{up}__Assembly__\n| Read alignment {CLR}\n| Transcriptome assembly {CLR}\n| Merging assemblies {CLR}\n| Writing fasta ' \
                   f'file for spliced exons {CLR}\n'
        steps = {1: assembly.find("Read alignment"), 2: assembly.find("Transcriptome assembly"),
                 3: assembly.find("Merging"), 4: assembly.find("Writing")}
        arrow_pos = steps[step]
        assembly = self.__add_arrow(arrow_pos, assembly)
        database = f''
        mode_panel = {'assembly': assembly}
        return mode_panel

    def show(self, step):
        panels = self.__create_panels(step)
        print(panels[self.mode])

    @staticmethod
    def __add_arrow(arrow_pos, mode):
        return mode[:arrow_pos-2] + '----> | ' + mode[arrow_pos:]


if __name__ == '__main__':
    panel = ProgressPanel('assembly')
    panel.show(step=1)
    panel.show(step=2)
