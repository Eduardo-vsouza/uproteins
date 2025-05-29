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


import os
import sys


class TSVConverter(object):
    def __init__(self, folder):
        self.path = sys.path[0]
        self.folder = folder

    def convert_files(self):
        if not os.path.exists(f'{self.folder}/tsv_msgf'):
            cmd_dir = f'mkdir {self.folder}/tsv_msgf'
            os.system(cmd_dir)
        files = os.listdir(f'{self.folder}')
        for file in files:
            if '.mzid' in file:
                cmd = f'java -cp {self.path}/dependencies/MSGF/MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i ' \
                      f'{self.folder}/{file} -o {self.folder}/tsv_msgf/{file}.tsv -showQValue 1'
                os.system(cmd)
        return self
