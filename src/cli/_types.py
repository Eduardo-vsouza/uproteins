"""Utilitarian module containing a number of helper functions designed to serve
as values for the :param:`type` parameter of the :mod:`argparse` module.

Those functions raise :exc:`ArgumentError` when the str value passed to them
fails to conform to their expected format, else convert the value to the
pertinent type and return it. Pathes are returned in their absolute forms.

All functions can receive a None value, in which case they just return None.
"""

import argparse
import pathlib
import shutil
import re


def Executable(val: str | None) -> pathlib.Path | None:
    """Receive a str path and raise an :exc:`ArgumentError` if the value is not
    a valid executable, else return a :obj:`Path` object.
    """
    if val is None:
        return val

    which = shutil.which(val)
    if which is None:
        raise argparse.ArgumentError
    return pathlib.Path(which).absolute()


def FilePath(val: str | None) -> pathlib.Path | None:
    """Receive a str path and raise an :exc:`ArgumentError` if the value is not
    the path to a file, else return a :obj:`Path` object.
    """
    if val is None:
        return val

    path = pathlib.Path(val)
    if not (path.exists() and path.is_file()):
        raise argparse.ArgumentError
    return path.absolute()


def FileName(val: str | None) -> pathlib.Path | None:
    """Receive a str path and raise an :exc:`ArgumentError` if the value is not
    a filename, else return a :obj:`Path` object.

    A filename is considered to be a path to a file or a path that doesn't yet
    exist.
    """
    if val is None:
        return val

    path = pathlib.Path(val)
    if path.is_file() or not path.exists():
        raise argparse.ArgumentError
    return path.absolute()


def DirectoryPath(val: str | None) -> pathlib.Path | None:
    """Receive a str path and raise an :exc:`ArgumentError` if the value is not
    a valid directory, else return a :obj:`Path` object.
    """
    if val is None:
        return val

    path = pathlib.Path(val)
    if not (path.exists() and path.is_dir()):
        raise argparse.ArgumentError
    return path.absolute()


def DirectoryName(val: str | None) -> pathlib.Path | None:
    """Receive a str path and raise an :exc:`ArgumentError` if the value is not
    a valid directoryname, else return a :obj:`Path` object.

    A valid directoryname is considered to be a path to a directory or a path
    that doesn't yet exist.
    """
    if val is None:
        return val

    path = pathlib.Path(val)
    if path.exists() and not path.is_dir():
        raise argparse.ArgumentError
    return path.absolute()


def Codon(val: str | None) -> str | None:
    """Receive a str val and raise an :exc:`ArgumentError` if the value is not
    a valid codon, else return the str untouched.

    A valid codon is considered to be any three letter str composed only of
    the `ATCG` characteres.
    """
    if val is None:
        return val

    if re.fullmatch(r'[ATCG]{3}', val) is None:
        raise argparse.ArgumentError
    return val
