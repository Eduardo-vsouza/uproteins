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


import argparse
import itertools
import pathlib

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
    ['--reads1', 'tmp_file'],
    ['--reads2', 'tmp_file'],
    ['--single', 'tmp_file', '--reads1', 'tmp_file'],
    ['--single', 'tmp_file', '--reads2', 'tmp_file'],
    ['--single', 'tmp_file', '--reads1', 'tmp_file', '--reads2', 'tmp_file']
])
def bad_validate_assembly(request, assembly_args, tmp_file):
    args = assembly_args
    tmp_file = str(tmp_file)

    for arg in request.param:
        if arg == 'tmp_file':
            args.append(tmp_file)
        else:
            args.append(arg)
    return args


@pytest.fixture(params=[
    ['--reads1', 'tmp_file', '--reads2', 'tmp_file'],
    ['--single', 'tmp_file']
])
def good_validate_assembly(request, assembly_args, tmp_file):
    args = assembly_args
    tmp_file = str(tmp_file)

    for arg in request.param:
        if arg == 'tmp_file':
            args.append(tmp_file)
        else:
            args.append(arg)
    return args


@pytest.fixture(params=[
    ['--external_transcriptome', 'tmp_file'],
    ['--external_gtf', 'tmp_file'],
    ['--starts', 'ATC,ATG', '--stops', 'ATA,ATG']
])
def bad_validate_database(request, database_args, tmp_file):
    args = database_args
    tmp_file = str(tmp_file)

    for arg in request.param:
        if arg == 'tmp_file':
            args.append(tmp_file)
        else:
            args.append(arg)
    return args


@pytest.fixture(params=[
    [],
    ['--external_gtf', 'tmp_file', '--external_transcriptome', 'tmp_file'],
    ['--starts', 'ATC,ATA', '--stops', 'ATG,CTC'],
])
def good_validate_database(request, database_args, tmp_file):
    args = database_args
    tmp_file = str(tmp_file)

    for arg in request.param:
        if arg == 'tmp_file':
            args.append(tmp_file)
        else:
            args.append(arg)
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
    yes = ['yes'[:i] for i in range(1, 4)] + ['YES'[:i] for i in range(1, 4)]
    no = ['no'[:i] for i in range(1, 3)] + ['NO'[:i] for i in range(1, 3)]

    def test_Executable(self, inexistent_path):
        exec = _types.Executable('exec')
        assert exec('exec') == pathlib.Path('exec')
        with pytest.raises(argparse.ArgumentTypeError):
            exec(str(inexistent_path))

    def test_DirectoryName(self, tmp_path, tmp_file, inexistent_path):
        assert _types.DirectoryName(str(tmp_path)) == tmp_path.absolute()
        assert (
            _types.DirectoryName(str(inexistent_path)) ==
            inexistent_path.absolute()
        )
        with pytest.raises(TypeError):
            _types.DirectoryName(str(tmp_file))

    def test_DirectoryPath(self, tmp_path, inexistent_path, tmp_file):
        assert _types.DirectoryPath(str(tmp_path)) == tmp_path.absolute()
        with pytest.raises(TypeError):
            _types.DirectoryPath(str(inexistent_path))
        with pytest.raises(TypeError):
            _types.DirectoryPath(str(tmp_file))

    def test_FileName(self, tmp_path, inexistent_path, tmp_file):
        assert (
            _types.FileName(str(inexistent_path)) == inexistent_path.absolute()
        )
        assert _types.FileName(str(tmp_file)) == tmp_file.absolute()
        with pytest.raises(TypeError):
            _types.FileName(str(tmp_path))

    def test_FilePath(self, tmp_path, inexistent_path, tmp_file):
        assert _types.FilePath(str(tmp_file)) == tmp_file.absolute()
        with pytest.raises(TypeError):
            _types.FilePath(str(inexistent_path))
        with pytest.raises(TypeError):
            _types.FilePath(str(tmp_path))

    @pytest.mark.parametrize('codon', codons)
    def test_Codon(self, codon):
        assert _types.Codon(codon) == codon

    @pytest.mark.parametrize('codon', bad_codons)
    def test_bad_Codon(self, codon):
        with pytest.raises(TypeError):
            _types.Codon(codon)

    def test_CommaList(self, tmp_path, inexistent_path):
        DirList = _types.CommaList(_types.DirectoryPath)
        path_list = DirList(str(tmp_path) + ',' + str(tmp_path))
        assert path_list is not None
        for path in path_list:
            assert path == tmp_path.absolute()
        with pytest.raises(argparse.ArgumentTypeError):
            DirList(str(tmp_path) + ',' + str(inexistent_path))

    @pytest.mark.parametrize('yes', yes)
    def test_yes_YesOrNoBooleanAction(self, yes):
        parser = argparse.ArgumentParser()
        namespace = argparse.Namespace()
        yes_or_no = _types.YesOrNoBooleanAction('act', 'act')
        yes_or_no(parser, namespace, yes)
        assert namespace.act

    @pytest.mark.parametrize('no', no)
    def test_no_YesOrNoBooleanAction(self, no):
        parser = argparse.ArgumentParser()
        namespace = argparse.Namespace()
        yes_or_no = _types.YesOrNoBooleanAction('act', 'act')
        yes_or_no(parser, namespace, no)
        assert not namespace.act

    def test_bad_YesOrNoBooleanAction(self):
        parser = argparse.ArgumentParser()
        namespace = argparse.Namespace()
        yes_or_no = _types.YesOrNoBooleanAction('act', 'act')
        with pytest.raises(SystemExit) as excinfo:
            yes_or_no(parser, namespace, 'B')
        assert excinfo.value.code == 2


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
