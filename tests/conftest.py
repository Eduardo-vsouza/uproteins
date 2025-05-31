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
import random
import typing as t

import pytest


@pytest.fixture
def tmp_file(tmp_path: pathlib.Path) -> pathlib.Path:
    file = (tmp_path / 'file.txt')
    file.write_text('')
    return file


@pytest.fixture
def inexistent_path() -> pathlib.Path:
    while True:
        random_nums = [str(random.randint(1, 10)) for _ in range(10)]
        path = pathlib.Path(''.join(random_nums))
        if not path.exists():
            return path


@pytest.fixture
def assembly_args(tmp_path, tmp_file) -> t.Generator[list[str], None, None]:
    """Return the smallest number of args necessary for the assembly parser.

    Note: Doesn't include `--single ...` or  `--reads1 ... --reads2 ...`, so
    it fails :func:`validate_assembly`.

    <b>Equivalent to the command-line input:</b>
    .. code-block:: text
        assembly --outdir tmp_path --genome tmp_file --gtf tmp_file
    """
    args = [
        'assembly',
        '--outdir', str(tmp_path),
        '--genome', str(tmp_file),
        '--gtf', str(tmp_file),
    ]

    def reset():
        nonlocal args
        args = [
            'assembly',
            '--outdir', str(tmp_path),
            '--genome', str(tmp_file),
            '--gtf', str(tmp_file),
        ]

    yield args
    reset()


@pytest.fixture
def database_args(tmp_path, tmp_file) -> t.Generator[list[str], None, None]:
    """Return the smallest number of args necessary for the database parser.

    <b>Equivalent to the command-line input:</b>
    .. code-block:: text
        database --outdir tmp_path --genome tmp_file --proteome tmp_file
    """
    args = [
        'database',
        '--outdir', str(tmp_path),
        '--genome', str(tmp_file),
        '--proteome', str(tmp_file),
    ]

    def reset():
        nonlocal args
        args = [
            'database',
            '--outdir', str(tmp_path),
            '--genome', str(tmp_file),
            '--proteome', str(tmp_file),
        ]

    yield args
    reset()
