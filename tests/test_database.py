from src import uproteins, cli, assembly  # noqa: F401


def test_database_parser(database_args, tmp_file):
    database_args += ['--genome', str(tmp_file), '--proteome', str(tmp_file)]
    parser, subparsers = cli.get_parsers()
    args = parser.parse_args(database_args)
    cli.validate_database(args, subparsers['database'])
