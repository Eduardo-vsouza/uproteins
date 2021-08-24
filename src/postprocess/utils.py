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
