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

from src import uproteins, cli, assembly  # noqa: F401
from tests import resources


def test_assembly_parser(assembly_args, tmp_file):
    assembly_args += ['--single', str(tmp_file)]
    parser, subparsers = cli.get_parsers()
    args = parser.parse_args(assembly_args)
    cli.validate_assembly(args, subparsers['assembly'])


def test_read_type(assembly_args, tmp_file):
    single = [*assembly_args, '--single', str(tmp_file)]
    paired = [
        *assembly_args, '--reads1', str(tmp_file), '--reads2', str(tmp_file)
    ]

    parser, subparsers = cli.get_parsers()
    single_args = parser.parse_args(single)
    cli.validate_assembly(single_args, subparsers['assembly'])
    single_mapper = assembly.ReadMapper(single_args)
    assert not single_mapper.is_paired

    parser, subparsers = cli.get_parsers()
    paired_args = parser.parse_args(paired)
    cli.validate_assembly(paired_args, subparsers['assembly'])
    paired_mapper = assembly.ReadMapper(paired_args)
    assert paired_mapper.is_paired


def test_assembly_mode(tmp_path):
    genome = rsrc.files(resources).joinpath("genome.fasta")
    gtf = rsrc.files(resources).joinpath("mtb.gtf")
    read1 = rsrc.files(resources).joinpath("ERR262980.fastq")
    read2 = rsrc.files(resources).joinpath("ERR262982.fastq")
    read3 = rsrc.files(resources).joinpath("ERR262983.fastq")

    args = [
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
    uproteins(args)

    assembled_path: pathlib.Path = tmp_path / 'assembled.gtf'
    transcripts_path: pathlib.Path = tmp_path / 'HISAT' / 'transcripts.fasta'

    # Make sure the expected files exist
    assert assembled_path.is_file()
    assert transcripts_path.is_file()

    # Make sure they have the expected results
    with rsrc.path(resources, 'results') as results:
        ok_assembled_path = results / 'assembled.gtf'
        ok_transcripts_path = results / 'HISAT' / 'transcripts.fasta'

        # Down the pipeline, the assembled file has the column labels inserted
        # into it. For the equals method to work, we need those labels.
        assembled = pd.read_csv(
            assembled_path,
            comment='#',
            names=[
                'seqname',
                'source',
                'feature',
                'start',
                'end',
                'score',
                'strand',
                'frame',
                'attributes'
            ],
            sep='\t'
        )
        ok_assembled = pd.read_csv(ok_assembled_path, comment='#', sep='\t')
        assert assembled.equals(ok_assembled)
        assert transcripts_path.read_text() == ok_transcripts_path.read_text()
