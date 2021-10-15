import os
import sys
from src.database import database_generator as dg
from . import peptide_search as ps
from . import results_new_approach as pms
from . import the_visualizer as vis
# from . import orthologs as phylo
from . import f_translation as trans
from .working_runs import OrganizePlot
from .Digestion_sets import Digestion, Digested, PlotData
# import prowser.gene_organizer as gorg
# import prowser.browser_gui_v16 as prsr
# from .testing import uproteins_testing as test
from .percolator import Decoy
from .assembly import TranscriptAssembly, CompareTranscripts
from .master import Archives
from .database import Database
from .postprocess import PostPercolator, ExtendedInformation
from .sequtils.orflib import AltCodons
from .upstream import SDInspection
from .testing import PipelineTesting
from .postms import TSVConverter
from .pipelines import PostMSPipeline, ValidatePipeline
from .forest import ProteinFixer, PreFiltering, FeatureFishing, SpectralForest


pypath = sys.path[0]


def run_workflow(args):
    run = True
    # argfix = ArgumentsFixer(parser)
    # args = argfix.get_abspath()
    mode = args.mode
    if run:
        if not os.path.exists(args.outdir):
            cmd_dir = "mkdir %s" % args.outdir
            os.system(cmd_dir)
        os.chdir(args.outdir)
        # if args.smorfs is not None:
        #     if args.smorfs != "YES" or args.smorfs != "NO":
        #         print("\nInvalid argument. smorfs argument must be YES or NO.\n")
        #         return os.EX_USAGE

        if mode == "assembly":
            arxivs = Archives()
            # mapping = ReadMapper(args)
            # mapping.create_indexes()
            # mapping.align_reads()
            # mapping.sort_aln(i)
            assembly = TranscriptAssembly(args)
            assembly.assemble()
            assembly.create_gtf_list()
            gtf = assembly.merge_transcripts()
            arxivs.assembledGTF = gtf
            comparisons = CompareTranscripts(args, arxivs.assembledGTF)
            compare_dir = comparisons.run_gffcompare()
            arxivs.comparisonsDirectory = compare_dir
            sequences = comparisons.extract_sequences()
            arxivs.RNASequences = sequences

        elif mode == "testing":
            testing = PipelineTesting(args)
            testing.run()


        elif mode == "database":
            #genome = dg.OrfPrediction(args)
            #genome.identify_orfs()
            db = Database(args)
            db.translate()
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

        elif mode == "ms":
            genome = ps.PeptideSearch("Genome", args.Mass_spec, "genome_database.fasta", args)
            genome.peptide_identification()
            genome_decoy = Decoy(db="genome_database.fasta", db_type="Genome")
            genome_decoy.reverse_sequences().to_fasta()
            genome_decoy_search = ps.PeptideSearch("Genome", args.Mass_spec, "Genome/Percolator/Genome_decoy.fasta", args, decoy=True)
            genome_decoy_search.peptide_identification()
            #genome.peptide_filtering()
            if args.Transcriptome is not None:
                transcriptome = ps.PeptideSearch("Transcriptome", args.Mass_spec, "transcriptome_database.fasta", args)
                transcriptome.peptide_identification()
                transcriptome_decoy = Decoy(db="transcriptome_database.fasta", db_type="Transcriptome")
                transcriptome_decoy.reverse_sequences().to_fasta()
                transcriptome_decoy_search = ps.PeptideSearch("Transcriptome", args.Mass_spec, "Transcriptome/Percolator/Transcriptome_decoy.fasta", args, decoy=True)
                transcriptome_decoy_search.peptide_identification()
                #transcriptome.peptide_filtering()

        elif mode == "postms":
            """ newest method """
            genome = PostMSPipeline(args=args, filetype='genome', folder='Genome')
            genome.run()
            if args.Transcriptome == 'YES':
                transcriptome = PostMSPipeline(args=args, filetype='transcriptome', folder='Transcriptome')
                transcriptome.run()
            """ new method """

            print("DONE. Results are inside Genome or Transcriptome folders.")



        elif mode == "validate":
            validation = ValidatePipeline(args=args)
            validation.validate_genome()
            validation.validate_transcriptome()

        elif mode == "visualization":
            if not os.path.exists("RefSeq_Reading_Frames.ods"):
                print(
                    "\nWe are preparing the proteogenomics visualizer for its first time use. This will happen only once.\n")
                trans.find_proteome_rfs(args.genome, args.genbank, args.outdir)
            if args.type == "dna":
                if not os.path.exists("dna_results_with_rfs.ods"):
                    print("\nPreparing the genome unique ORFs for proper visualization\n")
                    trans.find_orfs_rfs("Genome/Results/Genome_unique_results_numbers.xls", args.genome, "dna")
                vis.prog_browser(args.type)
            elif args.type == "both":
                if not os.path.exists("both_results_with_rfs.ods"):
                    trans.find_orfs_rfs("Genome/Results/Both_results_numbers.xls", args.genome, "both")
                vis.prog_browser(args.type)

            elif args.type == "rna":
                if not os.path.exists("rna_results_with_rfs.ods"):
                    trans.find_orfs_rfs("Transcriptome/Results/Rna_unique_results_numbers.xls", args.genome, "rna")
                    pms.fix_frames()
                vis.prog_browser(args.type)

        elif mode == "prowser":
            if not os.path.exists("RefSeq_Reading_Frames.ods"):
                print(
                    "\nWe are preparing the proteogenomics visualizer for its first time use. This will happen only once.\n")
                trans.find_proteome_rfs(args.genome, args.genbank, args.outdir)
            if args.type == "genome":
                if not os.path.exists("dna_results_with_rfs.ods"):
                    print("\nPreparing the genome unique ORFs for proper visualization\n")
                    trans.find_orfs_rfs("Genome/Results/genome_unique_results_summarized.xls", args.genome, "dna")
                    orfome = gorg.ORFome("dna_results_with_rfs.ods", "dna")
                    orfome.define_gene_names()
                    orfome.get_gene_info()
                    orfome.get_specs_counts()
                    orfome.organize_df()
                prsr.run_prowser("dna")
            elif args.type == "transcriptome":
                if not os.path.exists("rna_results_with_rfs"):
                    print("\nPreparing the transcriptome unique ORFs for browsing.\n")
                    trans.find_orfs_rfs("Transcriptome/Results/transcriptome_unique_results_summarized.xls", args.genome,
                                        "rna")
                    pms.fix_frames()
                    orfome = gorg.ORFome("rna_results_with_rfs.ods", "rna")
                    orfome.define_gene_names()
                    orfome.get_gene_info()
                    orfome.get_specs_counts()
                    orfome.organize_df()
                prsr.run_prowser("rna")
            elif args.type == "both":
                if not os.path.exists("both_results_with_rfs.ods"):
                    trans.find_orfs_rfs("Genome/Results/both_summarized_final_results.xls", args.genome, "both")
                if not os.path.exists("both_df_to_visualize.ods"):
                    orfome = gorg.ORFome("both_results_with_rfs.ods", "both")
                    orfome.define_gene_names()
                    orfome.get_gene_info()
                    orfome.get_specs_counts()
                    orfome.organize_df()
                prsr.run_prowser("both")


            # elif args.type == "rna":
            #     if not os.path.exists("rna_results_with_rfs.ods"):
            #         trans.find_orfs_rfs("Transcriptome/Results/transcriptome_unique_results_summarized.xls", args.genome, "rna")
            #         pms.fix_frames()
            else:
                print("\nPlease inform a valid data subset.")

        # elif mode == "orthologs":
        #     database = phylo.BlastDB("Orthologs")
        #     database.download_assemblies()
        #     database.create_database()
        #     rna = phylo.Orthologs(args.organism, args, "transcriptome")
        #     rna.blast_sequences()
        #     rna.identify_hits()
        #     both = phylo.Orthologs(args.organism, args, "both")
        #     both.blast_sequences()
        #     both.identify_hits()
        #     dna = phylo.Orthologs(args.organism, args, "genome")
        #     dna.blast_sequences()
        #     dna.identify_hits()
        #     results = phylo.Motifs(args)
        #     results.find_motifs()

        #
        # elif mode == "digestion":
        #     if not args.coverage and not args.ups:
        #         print("\nPlease choose either --coverage and/or --ups in order to perform these analyzes. You may choose "
        #               "both.\n")
        #     else:
        #         peptides = Digestion("genome_database_no_anno.fasta")
        #         peptides.digest_orfs(args.enzymes)
        #         genome_trypsin = Digested("genome_database_no_anno_0/", "genome",
        #                                   "genome_database_no_anno.fasta", "Trypsin", args.proteome)
        #         genome_lysc = Digested("genome_database_no_anno_1/", "genome",
        #                                "genome_database_no_anno.fasta", "Lys_C", args.proteome)
        #         genome_argc = Digested("genome_database_no_anno_2/", "genome",
        #                                "genome_database_no_anno.fasta", "Arg_C", args.proteome)
        #         genome_gluc = Digested("genome_database_no_anno_3/", "genome",
        #                                "genome_database_no_anno.fasta", "Glu_C", args.proteome)
        #         enzymes = args.enzymes.split(",")
        #         if args.coverage is True:
        #             if "0" in enzymes:
        #                 genome_trypsin.virtual_coverage("genome_database_no_anno_3")
        #             if "1" in enzymes:
        #                 genome_lysc.virtual_coverage()
        #             if "2" in enzymes:
        #                 genome_argc.virtual_coverage()
        #             if "3" in enzymes:
        #                 genome_gluc.virtual_coverage()
        #         if args.ups is True:
        #             if "0" in enzymes:
        #                 genome_trypsin.ups_per_orf()
        #             if "1" in enzymes:
        #                 genome_lysc.ups_per_orf()
        #             if "2" in enzymes:
        #                 genome_argc.ups_per_orf()
        #             if "3" in enzymes:
        #                 genome_gluc.ups_per_orf()
        #         enzymes = PlotData(args.enzymes)
        #         enzymes.plot_up_orf()
        #         enzymes.plot_mixed_data()
        #         enzymes.plot_coverage_data()
        #
        # elif mode == "runs":
        #     genome = OrganizePlot(
        #         "Genome/Results/genome_unique_results_summarized.xls")
        #     genome.plot_data()
        #     if args.Transcriptome is not None:
        #         rna = OrganizePlot("Transcriptome/Results/transcriptome_unique_results_summarized.xls")
        #         rna.plot_data()
        #         both = OrganizePlot("Genome/Results/Summarized_final_results.ods")
        #         both.plot_data()

        else:
            print("Invalid mode. Please choose a valid mode.")
