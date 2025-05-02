from src import uproteins, cli, assembly  # noqa: F401


# TODO
# def test_assembly(assembly_args, tmp_file):
#     assembly_args += ['--single', str(tmp_file)]
#     uproteins(assembly_args)


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


# TODO
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
    assert (tmp_path / 'assembled.gtf').is_file()
    assert (tmp_path / 'transcripts.fasta').is_file()
