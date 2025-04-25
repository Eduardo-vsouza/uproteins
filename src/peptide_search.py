import os
import sys


class PeptideSearch(object):
    def __init__(self, database_type, ms_files_folder, orf_file, args, decoy=False):
        self.database_type = database_type
        self.ms_files_folder = ms_files_folder
        self.orf_file = orf_file
        self.args = args
        self.path = sys.path[0]
        self.decoy = decoy

    def peptide_identification(self):
        print("\nPerforming peptide search using %s database\n" % self.database_type)
        cmd_dir_ms = 'mkdir %s' % self.database_type
        os.system(cmd_dir_ms)
        cmd_copy_orf = 'cp %s %s/' % (self.orf_file, self.database_type)
        os.system(cmd_copy_orf)
        cmd_copy_db = f'cp {self.orf_file} {os.path.abspath(self.ms_files_folder)}/.'
        os.system(cmd_copy_db)

        # list of arguments passed by argparser
        # arg_list = []
        # for item in vars(self.args).items():
        #     item_list = [None, "Mass_spec", "outdir", "Transcriptome", "mode"]
        #     if item[0] not in item_list:
        #         arg_list.append("--%s %s" % (item[0], item[1]))
        # arg_string = ""
        # for i in range(len(arg_list)):
        #     arg_string += "%s " % arg_list[i]

        # cmd_pep_search = 'Rscript %s/mzid_workflow.R %s--database %s --folder %s'\
        #                  % (sys.path[0], arg_string, os.path.abspath(self.orf_file), self.ms_files_folder)
        # os.system(cmd_pep_search)
        self.loop_search()
        cmd_move = 'mv %s/*.mzid %s/' % (self.ms_files_folder, self.database_type)
        os.system(cmd_move)

    def loop_search(self):
        files = [i for i in os.listdir(os.path.abspath(self.ms_files_folder)) if i.endswith('mzML')]
        for file in files:
            self.run_msgf(file)
        return self

    def run_msgf(self, file):
        self._check_folder()
        output = ""
        if self.decoy:
            output = f" -o {self.args.mass_spec}/{file}_decoy.mzid"
        ms_args = ""
        item_list = [None, "Mass_spec", "outdir", "Transcriptome", "mode", 'skip_assembly', 'skip_db', 'skip_ms',
                     'skip_postms', 'skip_validation', 'gtf', 'single', 'reads1', 'reads2', 'strandness', 'gff_compare_path',
                     'gffread_path', 'genome', 'proteome', 'minsize', 'maxsize', 'starts', 'stops', 'threads']
        for arg in vars(self.args).items():
            if arg[0] not in item_list and arg[1] is not None:
                ms_args += f" -{arg[0]} {arg[1]}"
        db = os.path.abspath(self.orf_file)
        cmd = f'java -Xmx48G -jar {self.path}/dependencies/MSGF/MSGFPlus.jar -d {db}{output} -tda 0 -s {os.path.abspath(self.args.mass_spec)}/{file} -addFeatures 1{ms_args}'
        os.system(cmd)
        return self

    def _check_folder(self):
        if not os.path.exists('mzid'):
            os.system('mkdir mzid')
        if self.args.transcriptome and not os.path.exists('mzid/Transcriptome'):
            os.system('mkdir mzid/Transcriptome')
        if not os.path.exists('mzid/Genome'):
            os.system('mkdir mzid/Genome')


    def peptide_filtering(self):
        """ DEPRECATED FOR NOW """
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

