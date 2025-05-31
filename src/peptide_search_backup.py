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


class PeptideSearch(object):
    def __init__(self, database_type, ms_files_folder, orf_file, args):
        self.database_type = database_type
        self.ms_files_folder = ms_files_folder
        self.orf_file = orf_file
        self.args = args

    def peptide_identification(self):
        print("\nPerforming peptide search using %s database\n" % self.database_type)
        cmd_dir_ms = 'mkdir %s' % self.database_type
        os.system(cmd_dir_ms)
        cmd_copy_orf = 'cp %s %s/' % (self.orf_file, self.database_type)
        os.system(cmd_copy_orf)
        cmd_copy_db = f'cp {self.orf_file} {os.path.abspath(self.ms_files_folder)}/.'
        os.system(cmd_copy_db)

        # list of arguments passed by argparser
        arg_list = []
        for item in vars(self.args).items():
            item_list = [None, "Mass_spec", "outdir", "Transcriptome", "mode"]
            if item[0] not in item_list:
                arg_list.append("--%s %s" % (item[0], item[1]))
        arg_string = ""
        for i in range(len(arg_list)):
            arg_string += "%s " % arg_list[i]

        cmd_pep_search = 'Rscript %s/mzid_workflow.R %s--database %s --folder %s'\
                         % (sys.path[0], arg_string, os.path.abspath(self.orf_file), self.ms_files_folder)
        os.system(cmd_pep_search)
        cmd_move = 'mv %s/*.mzid %s/' % (self.ms_files_folder, self.database_type)
        os.system(cmd_move)

    def peptide_filtering(self):
        print("Filtering peptide identification")

        # list of arguments passed by argparser
        arg_list = []
        for item in vars(self.args).items():
            item_list = [None, "Mass_spec", "outdir", "Transcriptome", "mode"]
            if item[0] not in item_list:
                arg_list.append("--%s %s" % (item[0], item[1]))
        arg_string = ""
        for i in range(len(arg_list)):
            arg_string += "%s " % arg_list[i]

        cmd_ms_filter = 'Rscript %s/workflow_filter.R %s--folder %s/' % (sys.path[0], arg_string, self.database_type)
        os.system(cmd_ms_filter)
        cmd_cat = 'cat %s/*pepseq.txt > %s/clustered_peptides.txt' % (self.database_type, self.database_type)
        os.system(cmd_cat)

