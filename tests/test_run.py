# Copyright © 2025 Eduardo Vieira de Souza
# Copyright © 2025 Adriana Canedo
# Copyright © 2025 Cristiano Valim Bizarro
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


import pathlib
from importlib import resources as rsrc

import pandas as pd
import pytest

from src import uproteins
from tests import resources


# Since the validate method requires a venv, it is not being run at the moment
@pytest.mark.fullrun
def test_full_run(tmp_path):
    genome = rsrc.files(resources).joinpath("genome.fasta")
    proteome = rsrc.files(resources).joinpath("proteome.fasta")
    gtf = rsrc.files(resources).joinpath("mtb.gtf")
    mzml = rsrc.files(resources).joinpath("mzml")
    read1 = rsrc.files(resources).joinpath("ERR262980.fastq")
    read2 = rsrc.files(resources).joinpath("ERR262982.fastq")
    read3 = rsrc.files(resources).joinpath("ERR262983.fastq")
    rrna = rsrc.files(resources).joinpath("rrna.fna")
    results = rsrc.files(resources).joinpath("results")

    assembly_args = [
        'assembly',
        '--outdir', str(tmp_path),
        '--single',
            str(read1) + ',' +  # noqa: E131
            str(read2) + ',' +  # noqa: E131
            str(read3),         # noqa: E131
        '--genome', str(genome),
        '--gtf', str(gtf),
        '--strandness', 'F',
    ]
    uproteins(assembly_args)

    database_args = [
        'database',
        '--outdir', str(tmp_path),
        '--genome', str(genome),
        '--proteome', str(proteome),
        '--starts', 'ATG,GTG,TTG,CTG',
        '--minsize', '30',
        '--maxsize', '300',
        '--Transcriptome', 'YES'
    ]
    uproteins(database_args)

    ms_args = [
        'ms',
        '--outdir', str(tmp_path),
        '--Mass_spec', str(mzml),
        '--inst', '2',
        '--t', '800ppm',
        '--Transcriptome', 'YES'
    ]
    uproteins(ms_args)

    postms_args = [
        'postms',
        '--outdir', str(tmp_path),
        '--genome', str(genome),
        '--proteome', str(proteome),
        '--rrna', str(rrna),
        '--gff', str(gtf),
        '--Mass_spec', str(mzml),
        '--Transcriptome', 'YES'
    ]
    uproteins(postms_args)

    validate_args = [
        'validate',
        '--outdir', str(tmp_path),
        '--transcriptome'
    ]
    # uproteins(validate_args)

    # Paths to the output
    genome_pre_path: pathlib.Path = (
        tmp_path
        / 'Genome'
        / 'Results'
        / 'genome_pre_validation_results.txt'
    )
    genome_post_path: pathlib.Path = (
        tmp_path
        / 'Genome'
        / 'Results'
        / 'genome_post_validation_results.txt'
    )

    transcriptome_pre_path: pathlib.Path = (
        tmp_path
        / 'Transcriptome'
        / 'Results'
        / 'transcriptome_pre_validation_results.txt'
    )
    transcriptome_post_path: pathlib.Path = (
        tmp_path
        / 'Transcriptome'
        / 'Results'
        / 'transcriptome_post_validation_results.txt'
    )

    # Make sure the result files were created
    assert genome_pre_path.is_file()
    # assert genome_post.is_file()
    assert transcriptome_pre_path.is_file()
    # assert transcriptome_post.is_file()

    # And that they have the expect results
    ok_genome_pre_path = (
        results
        / 'Genome'
        / 'Results'
        / 'genome_pre_validation_results.txt'
    )
    ok_genome_post_path = (
        results
        / 'Genome'
        / 'Results'
        / 'genome_post_validation_results.txt'
    )

    ok_transcriptome_pre_path = (
        results
        / 'Transcriptome'
        / 'Results'
        / 'transcriptome_pre_validation_results.txt'
    )
    ok_transcriptome_post_path = (
        results
        / 'Transcriptome'
        / 'Results'
        / 'transcriptome_post_validation_results.txt'
    )

    genome_pre = pd.read_csv(
        genome_pre_path,
        comment='#',
        sep='\t'
    )
    ok_genome_pre = pd.read_csv(
        ok_genome_pre_path,  # pyright: ignore
        comment='#',
        sep='\t'
    )

    transcriptome_pre = pd.read_csv(
        transcriptome_pre_path,
        comment='#',
        sep='\t'
    )
    ok_transcriptome_pre = pd.read_csv(
        ok_transcriptome_pre_path,  # pyright: ignore
        comment='#',
        sep='\t'
    )

    pd.testing.assert_frame_equal(
        genome_pre,
        ok_genome_pre,
        check_like=True
    )
    pd.testing.assert_frame_equal(
        transcriptome_pre,
        ok_transcriptome_pre,
        check_like=True
    )


@pytest.mark.fullrun
def test_external_assembly(tmp_path):
    assembled = rsrc.files(resources).joinpath("assembled.gtf")
    transcripts = rsrc.files(resources).joinpath("transcripts.fasta")
    genome = rsrc.files(resources).joinpath("genome.fasta")
    proteome = rsrc.files(resources).joinpath("proteome.fasta")
    gtf = rsrc.files(resources).joinpath("mtb.gtf")
    mzml = rsrc.files(resources).joinpath("mzml")
    rrna = rsrc.files(resources).joinpath("rrna.fna")
    results = rsrc.files(resources).joinpath("results")

    database_args = [
        'database',
        '--outdir', str(tmp_path),
        '--genome', str(genome),
        '--proteome', str(proteome),
        '--starts', 'ATG,GTG,TTG,CTG',
        '--minsize', '30',
        '--maxsize', '300',
        '--Transcriptome', 'YES',
        '--external_transcriptome', str(transcripts),
        '--external_gtf', str(assembled),
    ]
    uproteins(database_args)

    ms_args = [
        'ms',
        '--outdir', str(tmp_path),
        '--Mass_spec', str(mzml),
        '--inst', '2',
        '--t', '800ppm',
        '--Transcriptome', 'YES'
    ]
    uproteins(ms_args)

    postms_args = [
        'postms',
        '--outdir', str(tmp_path),
        '--genome', str(genome),
        '--proteome', str(proteome),
        '--rrna', str(rrna),
        '--gff', str(gtf),
        '--Mass_spec', str(mzml),
        '--Transcriptome', 'YES'
    ]
    uproteins(postms_args)

    validate_args = [
        'validate',
        '--outdir', str(tmp_path),
        '--transcriptome'
    ]
    # uproteins(validate_args)

    # Paths to the output
    genome_pre_path: pathlib.Path = (
        tmp_path
        / 'Genome'
        / 'Results'
        / 'genome_pre_validation_results.txt'
    )
    genome_post: pathlib.Path = (
        tmp_path
        / 'Genome'
        / 'Results'
        / 'genome_post_validation_results.txt'
    )

    transcriptome_pre_path: pathlib.Path = (
        tmp_path
        / 'Transcriptome'
        / 'Results'
        / 'transcriptome_pre_validation_results.txt'
    )
    transcriptome_post: pathlib.Path = (
        tmp_path
        / 'Transcriptome'
        / 'Results'
        / 'transcriptome_post_validation_results.txt'
    )

    # Make sure the result files were created
    assert genome_pre_path.is_file()
    # assert genome_post.is_file()
    assert transcriptome_pre_path.is_file()
    # assert transcriptome_post.is_file()

    # And that they have the expect results
    ok_genome_pre_path = (
        results
        / 'Genome'
        / 'Results'
        / 'genome_pre_validation_results.txt'
    )
    ok_genome_post = (
        results
        / 'Genome'
        / 'Results'
        / 'genome_post_validation_results.txt'
    )

    ok_transcriptome_pre_path = (
        results
        / 'Transcriptome'
        / 'Results'
        / 'transcriptome_pre_validation_results.txt'
    )
    ok_transcriptome_post = (
        results
        / 'Transcriptome'
        / 'Results'
        / 'transcriptome_post_validation_results.txt'
    )

    genome_pre = pd.read_csv(
        genome_pre_path,
        comment='#',
        sep='\t'
    )
    ok_genome_pre = pd.read_csv(
        ok_genome_pre_path,  # pyright: ignore
        comment='#',
        sep='\t'
    )

    transcriptome_pre = pd.read_csv(
        transcriptome_pre_path,
        comment='#',
        sep='\t'
    )
    ok_transcriptome_pre = pd.read_csv(
        ok_transcriptome_pre_path,  # pyright: ignore
        comment='#',
        sep='\t'
    )

    pd.testing.assert_frame_equal(
        genome_pre,
        ok_genome_pre,
        check_like=True
    )
    pd.testing.assert_frame_equal(
        transcriptome_pre,
        ok_transcriptome_pre,
        check_like=True
    )
