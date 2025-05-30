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
