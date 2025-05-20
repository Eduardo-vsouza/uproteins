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
    args = [
        '--outdir', str(tmp_path),
        'assembly',
        '--single', 'tests/resources/testing_reads.fastq',
        '--genome', 'tests/resources/testing_genome.fasta',
        '--gtf', 'tests/resources/test_gtf.gtf',
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

        assembled = pd.read_csv(assembled_path, comment='#', header=1)
        ok_assembled = pd.read_csv(ok_assembled_path, comment='#', header=1)
        assert assembled.equals(ok_assembled)
        assert transcripts_path.read_text() == ok_transcripts_path.read_text()
