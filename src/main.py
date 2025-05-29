# Copyright © 2021-2025 Eduardo Vieira de Souza
# Copyright © 2021-2025 Adriana Canedo
# Copyright © 2021-2025 Cristiano Valim Bizarro
# Copyright © 2025 Bruno Maestri A Becker
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


import argparse
import os
import sys

from src.cli import validators
from src.database import database_generator as dg
from src import peptide_search as ps
from src.percolator import Decoy
from src.assembly import TranscriptAssembly, CompareTranscripts, ReadMapper
from src.master import Archives
from src.database import Database
from src.postprocess import ResultsSummarizer
from src.testing import PipelineTesting
from src.pipelines import PostMSPipeline, ValidatePipeline
from src.sequtils.__helpers import ExternalAssemblyError
from src.metrics import Metrics


pypath = sys.path[0]


def run_workflow(
    args,
    subparser: argparse.ArgumentParser
):
    # argfix = ArgumentsFixer(parser)
    # args = argfix.get_abspath()
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
        validators.validate_assembly(args, subparser)

        arxivs = Archives()
        if args.step <= 1:
            mapping = ReadMapper(args)
            mapping.create_indexes()
            mapping.align_reads()
        # mapping.sort_aln(i)
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
        validators.validate_database(args, subparser)

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
        if args.transcriptome:
            print("Generating the transcriptome database.")
            transcriptome_db = dg.Database("transcriptome_ORFs.fasta", args.proteome, "transcriptome")
            transcriptome_db.mark_annotated()
            transcriptome_db.unify()
            print("Transcriptome database generated.")

    elif mode == "ms":
        genome = ps.PeptideSearch("Genome", args.mass_spec, "genome_database.fasta", args)
        genome.peptide_identification()
        genome_decoy = Decoy(db="genome_database.fasta", db_type="Genome")
        genome_decoy.reverse_sequences().to_fasta()
        genome_decoy_search = ps.PeptideSearch("Genome", args.mass_spec, "Genome/Percolator/Genome_decoy.fasta", args, decoy=True)
        genome_decoy_search.peptide_identification()
        if args.transcriptome:
            transcriptome = ps.PeptideSearch("Transcriptome", args.mass_spec, "transcriptome_database.fasta", args)
            transcriptome.peptide_identification()
            transcriptome_decoy = Decoy(db="transcriptome_database.fasta", db_type="Transcriptome")
            transcriptome_decoy.reverse_sequences().to_fasta()
            transcriptome_decoy_search = ps.PeptideSearch("Transcriptome", args.mass_spec, "Transcriptome/Percolator/Transcriptome_decoy.fasta", args, decoy=True)
            transcriptome_decoy_search.peptide_identification()

    elif mode == "postms":
        """ newest method """
        validators.validate_postms(args, subparser)

        genome = PostMSPipeline(args=args, filetype='genome', folder='Genome')
        genome.run()
        if args.transcriptome:
            transcriptome = PostMSPipeline(args=args, filetype='transcriptome', folder='Transcriptome')
            transcriptome.run()
        """ new method """

    elif mode == "validate":
        validation = ValidatePipeline(args=args)
        validation.validate_genome()
        validation.validate_transcriptome()
        summ = ResultsSummarizer(args.transcriptome)
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

    print(f'uProteInS {mode} step completed.')
