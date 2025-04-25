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

    elif single_present and reads1_present or reads2_present:
        parser.error(
            'argument --single: '
            'not allowed with arguments --reads1 and --reads2'
        )


def validate_database(
    args: argparse.Namespace,
    parser: argparse.ArgumentParser
) -> None:
    # This check only pass if one is None and the other isn't
    if (args.external_gtf is None) != (args.external_transcriptome is None):
        parser.error(
            '--external-transcriptome and '
            '--external-gtf must be given together'
        )
