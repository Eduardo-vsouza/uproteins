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


"""Those functions perform more complex validations of command line arguments
that can't be performed by simply configuring the parser.

If some arg fails to be validated, the function calls
:class:`ArgumentParser`.:func:`error`.
"""

import argparse


def validate_assembly(
    args: argparse.Namespace,
    parser: argparse.ArgumentParser
) -> None:
    """Make sure that the assembly args were correctly given.

    This function exists the cli with an error message if no read is present,
    both a single read and any of --reads1 or --reads2 aer present, or if any
    of --reads1 or --reads2 is present without the other.

    Arguments
    ---------
    args : Namespace
        The namespace containing the args, generated from the assembly parser.
    parser : ArgumentParser
        A reference to the assembly parser to invoke
            :class:`ArgumentParser`.:func:`error` from.
    """
    single_present = args.single is not None
    reads1_present = args.reads1 is not None
    reads2_present = args.reads2 is not None

    if not (single_present or reads1_present or reads2_present):
        parser.error(
            'either --single or --reads1 and --reads2 must be given'
        )

    if not single_present and not (reads1_present and reads2_present):
        parser.error(
            '--reads1 and --reads2 must be given together'
        )

    elif single_present and (reads1_present or reads2_present):
        parser.error(
            'argument --single: '
            'not allowed with arguments --reads1 and --reads2'
        )


def validate_database(
    args: argparse.Namespace,
    parser: argparse.ArgumentParser
) -> None:
    """Make sure that the datanase args were correctly given.

    This function exists the cli with an error message if either --external-gtf
    or --external-transcriptome is given without the other.

    This function exists the cli with an error message a codons appears as
    both a start codon and an end codon.

    Arguments
    ---------
    args : Namespace
        The namespace containing the args, generated from the database parser.
    parser : ArgumentParser
        A reference to the database parser to invoke
            :class:`ArgumentParser`.:func:`error` from.
    """
    # This check only pass if one is None and the other isn't
    if (args.external_gtf is None) != (args.external_transcriptome is None):
        parser.error(
            '--external-transcriptome and '
            '--external-gtf must be given together'
        )

    # We check if there is an intersection and capture it to use on the error
    # message
    if shared := set(args.starts).intersection(set(args.stops)):
        parser.error(
            "argument --starts: not allowed to share codons with argument "
            f"--stops: '{','.join(shared)}'"
        )


def validate_postms(
    args: argparse.Namespace,
    parser: argparse.ArgumentParser
) -> None:
    """Make sure that the postms args were correctly given.

    This function exists the cli with an error message a codons appears as
    both a start codon and an end codon.

    Arguments
    ---------
    args : Namespace
        The namespace containing the args, generated from the database parser.
    parser : ArgumentParser
        A reference to the database parser to invoke
            :class:`ArgumentParser`.:func:`error` from.
    """
    # We check if there is an intersection and capture it to use on the error
    # message
    if shared := set(args.starts).intersection(set(args.stops)):
        parser.error(
            "argument --starts: not allowed to share codons with argument "
            f"--stops: '{','.join(shared)}'"
        )
