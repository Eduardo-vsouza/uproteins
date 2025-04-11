import os
import sys

from src.database import database_generator as dg
from . import peptide_search as ps
from . import results_new_approach as pms
# from . import the_visualizer as vis
# from . import orthologs as phylo
from . import f_translation as trans
# from .working_runs import OrganizePlot
# from .Digestion_sets import Digestion, Digested, PlotData
# import prowser.gene_organizer as gorg
# import prowser.browser_gui_v16 as prsr
# from .testing import uproteins_testing as test
from .percolator import Decoy
from .assembly import TranscriptAssembly, CompareTranscripts, ReadMapper
from .master import Archives
from .database import Database
from .postprocess import PostPercolator, ExtendedInformation, ResultsWrapper, ResultsSummarizer
from .sequtils.orflib import AltCodons
from .upstream import SDInspection
from .testing import PipelineTesting
from .postms import TSVConverter
from .pipelines import PostMSPipeline, ValidatePipeline
from .forest import ProteinFixer, PreFiltering, FeatureFishing, SpectralForest
from src.sequtils.__helpers import ExternalAssemblyError
from src.metrics import Metrics


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
            if args.step <= 1:
                mapping = ReadMapper(args)
                mapping.create_indexes()
                mapping.align_reads()
            #mapping.sort_aln(i)
            assembly = TranscriptAssembly(args)
            if args.step <= 2:
                assembly.assemble()
            if args.step <= 3:
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
            if args.external_gtf is not None and args.external_transcriptome is not None:
                os.system(f'cp {args.external_gtf} {args.outdir}/assembled.gtf')
                if not os.path.exists(f'{args.outdir}/HISAT'):
                    os.system(f'mkdir {args.outdir}/HISAT')
                os.system(f'cp {args.external_transcriptome} {args.outdir}/HISAT/transcripts.fasta')
            elif args.external_gtf is not None and args.external_transcriptome is None:
                raise ExternalAssemblyError
            elif args.external_transcriptome is not None and args.external_gtf is None:
                raise ExternalAssemblyError
            # genome = dg.OrfPrediction(args)
            # genome.identify_orfs()
            print("Generating the genome database.")
            db = Database(args)
            db.translate()
            # print("\nORFs identified \nNow performing steps to generate the GENOME database\n")
            genome_db = dg.Database("genome_ORFs.fasta", args.proteome, "genome")
            genome_db.mark_annotated()
            genome_db.unify()
            print("Genome database generated.")
            if args.Transcriptome is not None:
                print("Generating the transcriptome database.")
                transcriptome_db = dg.Database("transcriptome_ORFs.fasta", args.proteome, "transcriptome")
                transcriptome_db.mark_annotated()
                transcriptome_db.unify()
                print("Transcriptome database generated.")
            print("Database generation step complete. Look for databases in %s" % args.outdir)

        elif mode == "ms":
            genome = ps.PeptideSearch("Genome", args.Mass_spec, "genome_database.fasta", args)
            genome.peptide_identification()
            genome_decoy = Decoy(db="genome_database.fasta", db_type="Genome")
            genome_decoy.reverse_sequences().to_fasta()
            genome_decoy_search = ps.PeptideSearch("Genome", args.Mass_spec, "Genome/Percolator/Genome_decoy.fasta", args, decoy=True)
            genome_decoy_search.peptide_identification()
            if args.Transcriptome is not None:
                transcriptome = ps.PeptideSearch("Transcriptome", args.Mass_spec, "transcriptome_database.fasta", args)
                transcriptome.peptide_identification()
                transcriptome_decoy = Decoy(db="transcriptome_database.fasta", db_type="Transcriptome")
                transcriptome_decoy.reverse_sequences().to_fasta()
                transcriptome_decoy_search = ps.PeptideSearch("Transcriptome", args.Mass_spec, "Transcriptome/Percolator/Transcriptome_decoy.fasta", args, decoy=True)
                transcriptome_decoy_search.peptide_identification()

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
            summ = ResultsSummarizer()
            summ.get_results()
            summ.save()

        elif mode == 'metrics':
            metrics = Metrics(args=args)
            metrics.count_sequences()
            metrics.plot_aa_dist()
            metrics.plot_orf_numbers()
            metrics.plot_energies()
            metrics.plot_rbs()
            metrics.generate_fasta()
            # metrics.plot_venn_2()

        else:
            print("Invalid mode. Please choose a valid mode.")
