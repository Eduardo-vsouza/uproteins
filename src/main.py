import os
import sys
from src.database import database_generator as dg
from . import peptide_search as ps
from . import results_new_approach as pms
from . import the_visualizer as vis
from . import orthologs as phylo
from . import f_translation as trans
from .postms import TSVConverter
from .working_runs import OrganizePlot
from .Digestion_sets import Digestion, Digested, PlotData
# import prowser.gene_organizer as gorg
# import prowser.browser_gui_v16 as prsr
from . import uProteins_testing as test
from .percolator import Decoy
from .assembly import TranscriptAssembly, CompareTranscripts
from .master import Archives
from .database import Database
from .postprocess import PostPercolator, ExtendedInformation
from .sequtils.postsearch import DecoyVoid
from .sequtils.orflib import AltCodons
from .upstream import SDInspection


pypath = sys.path[0]


def run_workflow(args):
    mode = args.mode
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
        test_args = test.PipelineTesting(args.outdir, args.skip_assembly, args.skip_db, args.skip_ms,
                                         args.skip_postms)
        state = test.TestingOutput(args)
        if args.skip_assembly != "TRUE":
            test_args.test_assembly()
            state.assembly = "OK"
        ass_out = state.testing_stdout("assembly step ")
        print(ass_out)
        if args.skip_db != "TRUE":
            test_args.summarize_assembly()
            test_args.test_database()
            db_content = state.db_checkpoint()
            if db_content == "Filled":
                state.db = "OK"
            else:
                state.db = "Failed"
        db_out = state.testing_stdout("database step ")
        print(db_out)

        if args.skip_ms != "TRUE":
            test_args.test_ms()
            state.ms = "OK"
            ms_content = state.ms_checkpoint()
            if ms_content == "Filled":
                state.ms = "OK"
            else:
                state.ms = "Failed"
        ms_out = state.testing_stdout("MS step ======")
        print(ms_out)

        if args.skip_postms != "TRUE":
            test_args.test_postms()
            postms_content = state.postms_checkpoint()
            if postms_content == "Filled":
                state.postms = "OK"
            else:
                state.postms = "Failed"
        postms_out = state.testing_stdout("PostMS step ==")
        print(postms_out)



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
        """ new method """
        # genome_perc = PercolatorProcessing("Genome", filetype="genome")
        # genome_perc.create_metafiles().convert_to_pin()
        # genome_perc.percolate()
        # if args.Transcriptome is not None:
        #     rna_perc = PercolatorProcessing("Transcriptome", filetype="transcriptome")
        #     rna_perc.create_metafiles().convert_to_pin()
        #     rna_perc.percolate()
        # genome_tsv = TSVConverter("Genome")
        # genome_tsv.convert_files()
        # if args.Transcriptome is not None:
        #     transcriptome_tsv = TSVConverter("Transcriptome")
        #     transcriptome_tsv.convert_files()
        # genome_filter = PostPercolator(args, "Genome", filetype='genome')
        # genome_filter.convert_output()
        # genome_filter.get_coordinates_genome()
        # genome_filter.filter_novel()
        # genome_filter.unique_peptides()
        # genome_filter.msgf_info()
        # genome_filter.protein_seqs()
        # genome_filter.protein_threshold()
        # genome_alts_pre_rf = AltCodons(file='Genome/post_perc/genome_results_02.txt', genome=args.genome,
        #                                maxsize=args.maxsize)
        # genome_alts_pre_rf.extend_orfs(args=args)
        # if args.rrna is not None:
        #     genome_rbs = SDInspection(args, filetype="genome", folder="Genome", alternatives=genome_alts_pre_rf.alternatives)
        #     alts = genome_rbs.get_free_energy()
        #     genome_alts_pre_rf.alternatives = alts
        # genome_alts_pre_rf.sort_by_coordinates()
        # genome_alts_pre_rf.sort_by_atg()
        # genome_alts_pre_rf.sort_by_shine()
        # genome_alts_pre_rf.sort_by_peptides()
        # priorities = genome_alts_pre_rf.get_priorities()
        # ext = ExtendedInformation(folder='Genome', filetype='genome', alternatives=genome_alts_pre_rf.alternatives)
        # ext.filter_alternatives(priorities)
        # ext.extract_spectra()
        if args.Transcriptome is not None:
            transcriptome_filter = PostPercolator(args, 'Transcriptome', filetype='transcriptome')
        #     transcriptome_filter.convert_output()
        #     transcriptome_filter.get_coordinates_rna()
        #     transcriptome_filter.filter_novel()
        #     transcriptome_filter.unique_peptides()
        #     transcriptome_filter.msgf_info()
        #     transcriptome_filter.protein_seqs()
        #     transcriptome_filter.protein_threshold()
            transcriptome_alts_pre_rf = AltCodons(file='Transcriptome/post_perc/transcriptome_results_02.txt',
                                                  genome=args.genome, maxsize=args.maxsize)
            transcriptome_alts_pre_rf.extend_orfs(args=args)
            if args.rrna is not None:
                transcriptome_rbs = SDInspection(args, filetype='transcriptome', folder='Transcriptome',
                                                 alternatives=transcriptome_alts_pre_rf.alternatives)
                rna_alts = transcriptome_rbs.get_free_energy()
                transcriptome_alts_pre_rf.alternatives = rna_alts
            transcriptome_alts_pre_rf.sort_by_coordinates()
            transcriptome_alts_pre_rf.sort_by_atg()
            transcriptome_alts_pre_rf.sort_by_shine()
            transcriptome_alts_pre_rf.sort_by_peptides()
            rna_priorities = transcriptome_alts_pre_rf.get_priorities()
            rna_ext = ExtendedInformation(folder='Transcriptome', filetype='transcriptome',
                                          alternatives=transcriptome_alts_pre_rf.alternatives)
            rna_ext.filter_alternatives(rna_priorities)
            rna_ext.extract_spectra()
    # decoy = DecoyVoid('Genome/genome_proteined.tsv', 'Transcriptome/transcriptome_proteined.tsv',
            #                   'genome_database.fasta', 'transcriptome_database.fasta')
            # decoy.check_decoy()
            # decoy.subset_single_df('all_subsets')
            # decoy.create_venn()
            # decoy.separate_subsets('.')
            # decoy.count_peptides()

        """ old"""
        # pms.merge_databases()
        # genome_results = pms.Results("Genome", "genome", "genome_database.fasta")
        # genome_results.rename_files()
        # genome_results.write_results()
        # genome_results.merge()
        # genome_results.add_orf()
        # genome_results.add_orf_sequence()
        # genome_results.total_spec_count()
        # genome_results.genome_coordinates()
        # genome_results.check_unique(args)
        # genome_results.coverage()
        # genome_results.summarize_results()
        # genome_results.summarize_dna_rna()
        # genome_results.create_fast("Genome")
        # genome_results.unique_orfs()
        #
        # # genome_results.minimal_results()
        # # genome_results.minimum_runs_compact("15", "genome")
        # if args.genbank is not None:
        #     genome_ncrnas = pms.GenomicContext(args, "genome", "Genome")
        #     genome_ncrnas.match_features()
        #     genome_ncrnas.add_to_df()
        # if args.Transcriptome == "YES":
        #     transcriptome_results = pms.Results("Transcriptome", "transcriptome",
        #                                         "transcriptome_database.fasta")
        #     transcriptome_results.rename_files()
        #     transcriptome_results.write_results()
        #     transcriptome_results.merge()
        #     transcriptome_results.add_orf()
        #     transcriptome_results.add_orf_sequence()
        #     transcriptome_results.total_spec_count()
        #     transcriptome_results.genome_coordinates()
        #     transcripts = pms.TranscriptLocalizer(args)
        #     transcripts.find_transcript()
        #     transcripts.fix_nucmer()
        #     transcripts.add_coords()
        #     transcriptome_results.check_unique(args)
        #     transcriptome_results.coverage()
        #     transcriptome_results.summarize_results()
        #     transcriptome_results.summarize_dna_rna()
        #     pms.create_venn()
        #     transcriptome_results.add_names()
        #     genome_results.create_fast("Both")
        #     transcriptome_results.create_fast("Rna")
        #     transcriptome_results.merge_both_results()
        #     transcriptome_results.unique_orfs()
        #
        #     both_results = pms.Results("Genome", "Both", "_")
        #     """
        #     both_results.minimal_results()
        #     both_results.minimum_runs_compact("15", "both")
        #     transcriptome_results.minimal_results()
        #     transcriptome_results.minimum_runs_compact("15", "transcriptome")
            
            # if args.genbank is not None:
            #     rna_features = pms.GenomicContext(args, "transcriptome", "Transcriptome")
            #     rna_features.match_features()
            #     rna_features.add_to_df()
            #     both_features = pms.GenomicContext(args, "both", "Genome")
            #     both_features.match_features()
            #     both_features.add_to_df()"""
        print("DONE. Results are inside Genome or Transcriptome folders.")

    # elif mode == "onestep":
    #     rna_seq = dg.Assembly(args)
    #     rna_seq.denovo_assembly()
    #     genome = dg.OrfPrediction(args)
    #     genome.Identify_ORFs()
    #     Cmd2 = 'Rscript replace_fr2.R'
    #     os.system(Cmd2)
    #     if args.Transcriptome is not None:
    #         transcriptome = dg.OrfPrediction(args)
    #         Cmd_Replace_RNA = 'Rscript replace_RNA.R'
    #         os.system(Cmd_Replace_RNA)
    #         transcriptome.Identify_ORFs()
    #     print("\nORFs identified \nNow performing steps to generate the GENOME database\n")
    #     genome_db = dg.Database("ORFs_replaced.fasta", args.proteome, "genome")
    #     genome_db.blast_to_Proteome()
    #     genome_db.extract_entries()
    #     genome_db.remove_entries()
    #     if args.Transcriptome is not None:
    #         print("\nGenome database generated. Now performing steps to generate the TRANSCRIPTOME database.\n")
    #         transcriptome_db = dg.Database("RNA_ORFs_replaced.fasta", args.proteome, "transcriptome")
    #         transcriptome_db.blast_to_Proteome()
    #         transcriptome_db.extract_entries()
    #         transcriptome_db.remove_entries()
    #     print("Database generation step complete. Look for databases in %s" % args.outdir)
    #     genome = ps.PeptideSearch("Genome", args.Mass_spec, "genome_database.fasta")
    #     genome.peptide_identification()
    #     genome.peptide_filtering()
    #     if args.Transcriptome is not None:
    #         transcriptome = ps.PeptideSearch("Transcriptome", args.Mass_spec, "transcriptome_database.fasta")
    #         transcriptome.peptide_identification()
    #         transcriptome.peptide_filtering()
    #     genome = pms.PostMS("Genome", args.genome, "genome")
    #     genome.cluster_peptides()
    #     genome.pep_to_fasta()
    #     genome.blast_to_orfs()
    #     genome.match_orfs()
    #     genome.remove_unmatched_orfs()
    #     genome.orf_pep_to_genome()
    #     genome.write_results()
    #     genome.add_pep_spec()
    #     if args.Transcriptome is not None:
    #         transcriptome = pms.PostMS("Transcriptome", args.genome, "transcriptome")
    #         transcriptome.cluster_peptides()
    #         transcriptome.pep_to_fasta()
    #         transcriptome.blast_to_orfs()
    #         transcriptome.match_orfs()
    #         transcriptome.remove_unmatched_orfs()
    #         transcriptome.orf_pep_to_genome()
    #         transcriptome.write_results()
    #         transcriptome.add_pep_spec()
    #     print("DONE. Results written to results_with_specs.txt in Genome or Transcriptome folders.")

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

    elif mode == "orthologs":
        database = phylo.BlastDB("Orthologs")
        database.download_assemblies()
        database.create_database()
        rna = phylo.Orthologs(args.organism, args, "transcriptome")
        rna.blast_sequences()
        rna.identify_hits()
        both = phylo.Orthologs(args.organism, args, "both")
        both.blast_sequences()
        both.identify_hits()
        dna = phylo.Orthologs(args.organism, args, "genome")
        dna.blast_sequences()
        dna.identify_hits()
        results = phylo.Motifs(args)
        results.find_motifs()


    elif mode == "digestion":
        if not args.coverage and not args.ups:
            print("\nPlease choose either --coverage and/or --ups in order to perform these analyzes. You may choose "
                  "both.\n")
        else:
            peptides = Digestion("genome_database_no_anno.fasta")
            peptides.digest_orfs(args.enzymes)
            genome_trypsin = Digested("genome_database_no_anno_0/", "genome",
                                      "genome_database_no_anno.fasta", "Trypsin", args.proteome)
            genome_lysc = Digested("genome_database_no_anno_1/", "genome",
                                   "genome_database_no_anno.fasta", "Lys_C", args.proteome)
            genome_argc = Digested("genome_database_no_anno_2/", "genome",
                                   "genome_database_no_anno.fasta", "Arg_C", args.proteome)
            genome_gluc = Digested("genome_database_no_anno_3/", "genome",
                                   "genome_database_no_anno.fasta", "Glu_C", args.proteome)
            enzymes = args.enzymes.split(",")
            if args.coverage is True:
                if "0" in enzymes:
                    genome_trypsin.virtual_coverage("genome_database_no_anno_3")
                if "1" in enzymes:
                    genome_lysc.virtual_coverage()
                if "2" in enzymes:
                    genome_argc.virtual_coverage()
                if "3" in enzymes:
                    genome_gluc.virtual_coverage()
            if args.ups is True:
                if "0" in enzymes:
                    genome_trypsin.ups_per_orf()
                if "1" in enzymes:
                    genome_lysc.ups_per_orf()
                if "2" in enzymes:
                    genome_argc.ups_per_orf()
                if "3" in enzymes:
                    genome_gluc.ups_per_orf()
            enzymes = PlotData(args.enzymes)
            enzymes.plot_up_orf()
            enzymes.plot_mixed_data()
            enzymes.plot_coverage_data()

    elif mode == "runs":
        genome = OrganizePlot(
            "Genome/Results/genome_unique_results_summarized.xls")
        genome.plot_data()
        if args.Transcriptome is not None:
            rna = OrganizePlot("Transcriptome/Results/transcriptome_unique_results_summarized.xls")
            rna.plot_data()
            both = OrganizePlot("Genome/Results/Summarized_final_results.ods")
            both.plot_data()

    else:
        print("Invalid mode. Please choose a valid mode.")
