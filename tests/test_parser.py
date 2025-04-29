import itertools

import pytest

from src import cli
from src.cli import _types


# We define some fixtures for validation parameters so that we can rerun the
# test function individualy for each potential bad and good parameter
# combinations
#
# Most parameter validation is made by argparse itself or the _types submodule,
# so `validate_mode` functions have fairly few cases they test for


@pytest.fixture(params=[
    [],
    ['--reads1'],
    ['--reads2'],
    ['--single', '--reads1'],
    ['--single', '--reads2'],
    ['--single', '--reads1', '--reads2'],
])
def bad_validate_assembly(request, assembly_args, tmp_file):
    args = assembly_args
    tmp_file = str(tmp_file)

    for arg in request.param:
        args.append(arg)
        args.append(tmp_file)
    return args


@pytest.fixture(params=[['--reads1', '--reads2'], ['--single']])
def good_validate_assembly(request, assembly_args, tmp_file):
    args = assembly_args
    tmp_file = str(tmp_file)

    for arg in request.param:
        args.append(arg)
        args.append(tmp_file)
    return args


@pytest.fixture(params=[['--external-transcriptome'], ['--external-gtf']])
def bad_validate_database(request, database_args, tmp_file):
    args = database_args
    tmp_file = str(tmp_file)

    for arg in request.param:
        args.append(arg)
        args.append(tmp_file)
    return args


@pytest.fixture(params=[[], ['--external-gtf', '--external-transcriptome']])
def good_validate_database(request, database_args, tmp_file):
    args = database_args
    tmp_file = str(tmp_file)

    for arg in request.param:
        args.append(arg)
        args.append(tmp_file)
    return args


class TestTypes:
    # For those functions, argparse deals with raising a SystemExit with
    # exitcode 2 (what we want) if they simply raise TypeError, so we just make
    # sure to test every possible scenario specific to _types themselves

    # This snippet generates a list with all possible 64 codons
    codons = [''.join(codon) for codon in itertools.product('ATCG', repeat=3)]
    bad_codons = [
        'ATD',
        'ADG',
        'DTG',
        'ATGGTA',
        'A*G',
        'aTG',
        'atg',
        'CA',
        'A',
        'AA',
        'ATGA'
    ]

    def test_Executable(self, tmp_path):
        _types.Executable('py')
        with pytest.raises(TypeError):
            _types.Executable(str(tmp_path))

    def test_DirectoryName(self, tmp_path, tmp_file, inexistent_path):
        _types.DirectoryName(str(tmp_path))
        _types.DirectoryName(str(inexistent_path))
        with pytest.raises(TypeError):
            _types.DirectoryName(str(tmp_file))

    def test_DirectoryPath(self, tmp_path, inexistent_path, tmp_file):
        _types.DirectoryPath(str(tmp_path))
        with pytest.raises(TypeError):
            _types.DirectoryPath(str(inexistent_path))
        with pytest.raises(TypeError):
            _types.DirectoryPath(str(tmp_file))

    def test_FileName(self, tmp_path, inexistent_path, tmp_file):
        _types.FileName(str(inexistent_path))
        _types.FileName(str(tmp_file))
        with pytest.raises(TypeError):
            _types.FileName(str(tmp_path))

    def test_FilePath(self, tmp_path, inexistent_path, tmp_file):
        _types.FilePath(str(tmp_file))
        with pytest.raises(TypeError):
            _types.FilePath(str(inexistent_path))
        with pytest.raises(TypeError):
            _types.FilePath(str(tmp_path))

    @pytest.mark.parametrize('codon', codons)
    def test_Codon(self, codon):
        _types.Codon(codon)

    @pytest.mark.parametrize('codon', bad_codons)
    def test_bad_Codon(self, codon):
        with pytest.raises(TypeError):
            _types.Codon(codon)


class TestValidate:
    def test_good_validate_assembly(self, good_validate_assembly):
        parser, subparsers = cli.get_parsers()
        parsed = parser.parse_args(good_validate_assembly)
        cli.validate_assembly(parsed, subparsers['assembly'])

    def test_bad_validate_assembly(self, bad_validate_assembly):
        parser, subparsers = cli.get_parsers()
        parsed = parser.parse_args(bad_validate_assembly)
        with pytest.raises(SystemExit) as excinfo:
            cli.validate_assembly(parsed, subparsers['assembly'])
        assert excinfo.value.code == 2

    def test_good_validate_database(self, good_validate_database):
        parser, subparsers = cli.get_parsers()
        parsed = parser.parse_args(good_validate_database)
        cli.validate_database(parsed, subparsers['database'])

    def test_bad_validate_database(self, bad_validate_database):
        parser, subparsers = cli.get_parsers()
        parsed = parser.parse_args(bad_validate_database)
        with pytest.raises(SystemExit) as excinfo:
            cli.validate_database(parsed, subparsers['database'])
        assert excinfo.value.code == 2
