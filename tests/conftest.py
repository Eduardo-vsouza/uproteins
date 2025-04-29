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
        --outdir tmp_path assembly --genome tmp_file --gtf tmp_file
        --gffcompare py --gffread py
    """
    tmp_path = str(tmp_path)
    tmp_file = str(tmp_file)

    args = [
        '--outdir', tmp_path,
        'assembly',
        '--genome', tmp_file,
        '--gtf', tmp_file,
        '--gffcompare', 'py',
        '--gffread', 'py',
    ]

    def reset():
        nonlocal args
        args = [
            '--outdir', tmp_path,
            'assembly',
            '--genome', tmp_file,
            '--gtf', tmp_file,
            '--gffcompare', 'py',
            '--gffread', 'py',
        ]

    yield args
    reset()


@pytest.fixture
def database_args(tmp_path, tmp_file) -> t.Generator[list[str], None, None]:
    """Return the smallest number of args necessary for the database parser.

    <b>Equivalent to the command-line input:</b>
    .. code-block:: text
        --outdir tmp_path database --genome tmp_file --proteome tmp_file
    """
    tmp_path = str(tmp_path)
    tmp_file = str(tmp_file)

    args = [
        '--outdir', tmp_path,
        'database',
        '--genome', tmp_file,
        '--proteome', tmp_file,
    ]

    def reset():
        nonlocal args
        args = [
            '--outdir', tmp_path,
            'database',
            '--genome', tmp_file,
            '--proteome', tmp_file,
        ]

    yield args
    reset()
