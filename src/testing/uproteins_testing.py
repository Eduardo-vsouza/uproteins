import os
import sys

from src.database import database_generator as dg
from src import peptide_search as ps
from src import results_new_approach as pms
from . import TestArgs
from ..assembly import TranscriptAssembly, CompareTranscripts
from ..master import Archives


class TestingOutput(object):
    def __init__(self, args):

        self.args = args
        self.boolDict = {"TRUE": "skipped", "FALSE": "not tested"}
        self.assembly = self.boolDict[self.args.skip_assembly]
        self.db = self.boolDict[self.args.skip_db]
        self.ms = self.boolDict[self.args.skip_ms]
        self.postms = self.boolDict[self.args.skip_postms]


    def db_checkpoint(self):
        """ Check whether the databases are empty or not. """
        out = os.path.abspath(self.args.outdir)
        content = None
        if os.path.exists(f'{out}/transcriptome_database.fasta'):
            if os.stat(f"{out}/transcriptome_database.fasta").st_size != 0:
                content = "Filled"
            else:
                content = "Empty"
        return content

    def ms_checkpoint(self):
        """ Checks whether the clustered_peptides file is at its proper place and is not empty. This is the final result
         of the MS step. """
        out = os.path.abspath(self.args.outdir)
        content = None
        file = f'{out}/Transcriptome/clustered_peptides.txt'
        if os.path.exists(file):
            if os.stat(file).st_size != 0:
                content = "Filled"
            else:
                content = "Empty"
        return content

    def postms_checkpoint(self):
        out = os.path.abspath(self.args.outdir)
        content = None
        file = f'{out}/Results/both_ORFs.fasta'
        if os.path.exists(file):
            if os.stat(file).st_size != 0:
                content = "Filled"
            else:
                content = "Empty"
        return content


    def testing_stdout(self, mode):

        stdout = f"======================== Testing uProteInS {mode}=======================\n" \
                                "|\n" \
                                f"|- Assembly: {self.assembly}\n" \
                                "|\n" \
                                f"|- Database: {self.db}\n" \
                                "|\n" \
                                f"|- MS: {self.ms}\n" \
                                "|\n" \
                                f"|- PostMS: {self.postms}\n" \
                                "================================================================================\n"
        return stdout

class PipelineTesting(object):
    def __init__(self, args):
        self.args = args

    def test_assembly(self):
        ...

    def check_skipping(self, mode):
        """
        :param mode: skip_assembly, skip_db
        :return:
        """
        if not self.args.skip_assembly:
            arxivs = Archives()
            assembly = TranscriptAssembly(self.args)
            assembly.assemble()
            assembly.create_gtf_list()
            gtf = assembly.merge_transcripts()
            arxivs.assembledGTF = gtf
            comparisons = CompareTranscripts(self.args, arxivs.assembledGTF)
            compare_dir = comparisons.run_gffcompare()
            arxivs.comparisonsDirectory = compare_dir
            sequences = comparisons.extract_sequences()
            arxivs.RNASequences = sequences



class PipelineTestingOld(object):
    def __init__(self, outdir, skip_assembly, skip_db, skip_ms, skip_postms):
        self.outdir = outdir
        self.path = sys.path
        self.testFolder = f"{self.path[0][:-3]}testing_kit"
        self.args = TestArgs(self.testFolder, self.outdir, skip_assembly, skip_db, skip_ms, skip_postms)

    def test_assembly(self):
        args = self.args
        rna_seq = dg.Assembly(args)
        rna_seq.denovo_assembly()

    def summarize_assembly(self):
        cmd_rename = f'mv {self.args.outdir}/trinity/Trinity.fasta {self.args.outdir}/trinity/Trinity_cut.fasta'
        os.system(cmd_rename)
        cmd_new_rna = f'cp {self.testFolder}/test_transcriptome.fasta {self.args.outdir}/trinity/Trinity.fasta'
        os.system(cmd_new_rna)

    def test_database(self):
        args = self.args
        genome = dg.OrfPrediction(args)
        genome.identify_orfs()
        print("\nORFs identified \nNow performing steps to generate the GENOME database\n")
        genome_db = dg.Database("genome_ORFs.fasta", args.proteome, "genome")
        genome_db.blast_to_Proteome()
        genome_db.extract_entries()
        cmd_cat_entries = 'cat genome_entries.txt genome_entries_short.txt > genome_entries_cat.txt'
        os.system(cmd_cat_entries)
        genome_db.remove_entries()
        genome_db.create_custom()
        if args.Transcriptome is not None:
            print("\nGenome database generated. Now performing steps to generate the TRANSCRIPTOME database.\n")
            transcriptome_db = dg.Database("transcriptome_ORFs.fasta", args.proteome, "transcriptome")
            transcriptome_db.blast_to_Proteome()
            transcriptome_db.extract_entries()
            cmd_cat_entries_rna = 'cat transcriptome_entries.txt transcriptome_entries_short.txt > ' \
                                  'transcriptome_entries_cat.txt'
            os.system(cmd_cat_entries_rna)
            transcriptome_db.remove_entries()
            transcriptome_db.create_custom()
        print("Database generation step complete. Look for databases in %s" % args.outdir)

    def test_ms(self):
        """ Tests the peptide search mode of uProteInS. """
        args = self.args.msArgs
        ms_folder = args.get_ms_folder(self.args.path)
        genome = ps.PeptideSearch("Genome", ms_folder, "genome_database.fasta", args)
        # genome.peptide_identification()
        genome.peptide_filtering()
        if args.Transcriptome is not None:
            transcriptome = ps.PeptideSearch("Transcriptome", ms_folder, "transcriptome_database.fasta", args)
            # transcriptome.peptide_identification()
            transcriptome.peptide_filtering()

    def test_postms(self):
        args = self.args
        pms.merge_databases()
        genome_results = pms.Results("Genome", "genome", "genome_database.fasta")
        genome_results.rename_files()
        genome_results.write_results()
        genome_results.merge()
        genome_results.add_orf()
        genome_results.add_orf_sequence()
        genome_results.total_spec_count()
        genome_results.genome_coordinates()
        genome_results.check_unique(args)
        genome_results.coverage()
        genome_results.summarize_results()
        genome_results.summarize_dna_rna()
        # genome_results.create_fast("Genome")
        # genome_results.unique_orfs()

        # genome_results.minimal_results()
        # genome_results.minimum_runs_compact("15", "genome")
        if args.genbank is not None:
            genome_ncrnas = pms.GenomicContext(args, "genome", "Genome")
            genome_ncrnas.match_features()
            genome_ncrnas.add_to_df()
        if args.Transcriptome == "YES":
            transcriptome_results = pms.Results("Transcriptome", "transcriptome",
                                                "transcriptome_database.fasta")
            transcriptome_results.rename_files()
            transcriptome_results.write_results()
            transcriptome_results.merge()
            transcriptome_results.add_orf()
            transcriptome_results.add_orf_sequence()
            transcriptome_results.total_spec_count()
            transcriptome_results.genome_coordinates()
            transcripts = pms.TranscriptLocalizer(args)
            transcripts.find_transcript()
            transcripts.fix_nucmer()
            transcripts.add_coords()
            transcriptome_results.check_unique(args)
            transcriptome_results.coverage()
            transcriptome_results.summarize_results()
            transcriptome_results.summarize_dna_rna()
            # pms.create_venn()
            # transcriptome_results.add_names()
            # genome_results.create_fast("Both")
            # transcriptome_results.create_fast("Rna")
            # transcriptome_results.merge_both_results()
            # transcriptome_results.unique_orfs()

            both_results = pms.Results("Genome", "Both", "_")
            """
            both_results.minimal_results()
            both_results.minimum_runs_compact("15", "both")
            transcriptome_results.minimal_results()
            transcriptome_results.minimum_runs_compact("15", "transcriptome")
            """
            # if args.genbank is not None:
            #     rna_features = pms.GenomicContext(args, "transcriptome", "Transcriptome")
            #     rna_features.match_features()
            #     rna_features.add_to_df()
            #     both_features = pms.GenomicContext(args, "both", "Genome")
            #     both_features.match_features()
            #     both_features.add_to_df()
