"""Utilitarian module containing a number of helper functions designed to serve
as values for the :param:`type` parameter of the :mod:`argparse` module.

Those functions raise :exc:`ArgumentError` when the str value passed to them
fails to conform to their expected format, else convert the value to the
pertinent type and return it. Pathes are returned in their absolute forms.

All functions can receive a None value, in which case they just return None.
"""

import pathlib
import shutil
import re
import typing as t


class Executable:
    def __init__(self, default: t.Optional[str]) -> None:
        """Executable type. This type is a callable class. Initialize it at the
        type field with the default value. The type checking won't fail even
        if it's not an Executable.

        When called
        -----------
        Receive a str path and raise an :exc:`ArgumentError` if the value
        is not a valid executable, else return a :obj:`Path` object.

        Create an :obj:`Executable` with `'exec'` as default.
        Passing 'exec' won't raise even if `'exec'` is not an executable.
        >>> exec = Executable('exec')
        >>> exec('exec')
        Path('path/to/exec')

        But passing another invalid executable raises.
        >>> exec('foo')
        TypeError

        If the arg is a valid executable, it doesn't raise.
        >>> exec('valid')
        Path('path/to/valid.exe')
        """
        self.default = default

    def __call__(self, val: t.Optional[str]) -> t.Optional[pathlib.Path]:
        """Receive a str path and raise an :exc:`ArgumentError` if the value
        is not a valid executable, else return a :obj:`Path` object.
        """
        if val is None:
            return val

        if val == self.default:
            return pathlib.Path(val)

        which = shutil.which(val)
        if which is None:
            raise TypeError

        return pathlib.Path(which).absolute()


def FilePath(val: t.Optional[str]) -> t.Optional[pathlib.Path]:
    """Receive a str path and raise an :exc:`ArgumentError` if the value is not
    the path to a file, else return a :obj:`Path` object.
    """
    if val is None:
        return val

    path = pathlib.Path(val)
    if not (path.exists() and path.is_file()):
        raise TypeError
    return path.absolute()


def FileName(val: t.Optional[str]) -> t.Optional[pathlib.Path]:
    """Receive a str path and raise an :exc:`ArgumentError` if the value is not
    a filename, else return a :obj:`Path` object.

    A filename is considered to be a path to a file or a path that doesn't yet
    exist.
    """
    if val is None:
        return val

    path = pathlib.Path(val)
    if path.exists() and not path.is_file():
        raise TypeError
    return path.absolute()


def DirectoryPath(val: t.Optional[str]) -> t.Optional[pathlib.Path]:
    """Receive a str path and raise an :exc:`ArgumentError` if the value is not
    a valid directory, else return a :obj:`Path` object.
    """
    if val is None:
        return val

    path = pathlib.Path(val)
    if not (path.exists() and path.is_dir()):
        raise TypeError
    return path.absolute()


def DirectoryName(val: t.Optional[str]) -> t.Optional[pathlib.Path]:
    """Receive a str path and raise an :exc:`ArgumentError` if the value is not
    a valid directoryname, else return a :obj:`Path` object.

    A valid directoryname is considered to be a path to a directory or a path
    that doesn't yet exist.
    """
    if val is None:
        return val

    path = pathlib.Path(val)
    if path.exists() and not path.is_dir():
        raise TypeError
    return path.absolute()


def Codon(val: t.Optional[str]) -> t.Optional[str]:
    """Receive a str val and raise an :exc:`ArgumentError` if the value is not
    a valid codon, else return the str untouched.

    A valid codon is considered to be any three letter str composed only of
    the `ATCG` characteres.
    """
    if val is None:
        return val

    if re.fullmatch(r'[ATCG]{3}', val) is None:
        raise TypeError
    return val
