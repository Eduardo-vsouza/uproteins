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
import sys
import os
import typing as t

from src.cli import _types


# TODO: Better help messages for each arg, likely removing the formatter_class
_parser = argparse.ArgumentParser(
    prog='uProteInS',
    description="Run ORF identification pipeline",
    fromfile_prefix_chars='@',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

# Reserving '-v' for a future '--verbose'
_parser.add_argument(
    '-V', '--version',
    action='version',
    version='%(prog)s 1.4.0 (GPL-3.0-or-later)'
)

# Subparsers
_modes = _parser.add_subparsers(
    title='modes',
    dest='mode',
    help='Pipeline mode to be run',
    required=True
)

# =============
# ASSEMBLY MODE
# =============
_assembly_parser = _modes.add_parser('assembly')
_assembly_parser.add_argument(
    '-o', '--outdir',
    help='Inform output directory',
    type=_types.DirectoryName,
    required=True
)
_assembly_parser.add_argument(
    "--genome", "-g",
    help="Path to the reference genome.",
    required=True,
    type=_types.FilePath
)
_assembly_parser.add_argument(
    "--gtf", "-G",
    help="Path to the GTF file. It does not accept a GFF3 format. If you "
    "only have a GFF file, convert it to GTF using gffread.",
    required=True,
    type=_types.FilePath
)
# Reserving '-s' for a future '--silently'
_assembly_parser.add_argument(
    '--single', '-S',
    type=_types.CommaList(_types.FilePath),
    help="Single-end reads. Comma separated list."
)
_assembly_parser.add_argument(
    '--reads1', '-1',
    type=_types.CommaList(_types.FilePath),
    help="Paired-end reads. Inform your reads2 as well. Comma separated "
    "list."
)
_assembly_parser.add_argument(
    '--reads2', '-2',
    type=_types.CommaList(_types.FilePath),
    help="Paired-end reads. Inform your reads1 as well. Comma separated "
    "list."
)
_assembly_parser.add_argument(
    "--strandness", "-L",
    choices=('RF', 'FR', 'F', 'R'),
    help="If your RNA-Seq experiment is strand-specific, specify the "
    "orientation. If paired-end reads: RF or FR. If single-end reads: F "
    "or R.",
)
_assembly_parser.add_argument(
    "--threads",
    help="Number of threads for samtools and stringtie to use.",
    type=_types.PositiveInt
)
_assembly_parser.add_argument(
    '--memory',
    help="Maximum memory per thread for samtools sort to use.",
    type=_types.Memory

)
_assembly_parser.add_argument(
    "--gffcompare_path",
    help="Path to gffcompare executable. If not provided, be sure to "
    "have it on your PATH, else the pipeline raises an error.",
    default='gffcompare',
    type=_types.Executable('gffcompare')
)
_assembly_parser.add_argument(
    "--gffread_path",
    help="Path to gffread executable. If not provided, be sure to "
    "have it on your PATH, else the pipeline raises an error.",
    default='gffread',
    type=_types.Executable('gffread')
)
_assembly_parser.add_argument(
    "--step",
    help="Assembly step to skip to. Developer parameter.",
    type=int,
    default=0
)

# =============
# DATABASE MODE
# =============
_database_parser = _modes.add_parser('database')
_database_parser.add_argument(
    '-o', '--outdir',
    help='Inform output directory',
    type=_types.DirectoryName,
    required=True
)
_database_parser.add_argument(
    "--genome", "-g",
    help="Genome fasta file to be translated to the six reading frames",
    required=True,
    type=_types.FilePath
)
_database_parser.add_argument(
    "--proteome", "-p",
    help="Proteome fasta file, i.e. 'Uniprot.fasta'",
    required=True,
    type=_types.FilePath
)
_database_parser.add_argument(
    "--starts",
    help="The start codons that should be used for the six and three "
    "frame translation. By default, the NCBI bacterial genetic code is "
    "used. Space separated list.",
    type=_types.CommaList(_types.Codon),
    default='TTG,CTG,ATT,ATC,ATA,ATG,GTG'
)
_database_parser.add_argument(
    "--stops",
    help="The stop codons that should be used for the six and three frame "
    "translation. Space separated list.",
    type=_types.CommaList(_types.Codon),
    default='TAA,TAG,TGA'
)
_database_parser.add_argument(
    "--Transcriptome", "-T",
    action=_types.YesOrNoBooleanAction,
    dest='transcriptome',
    help="Generate transcriptome database. A YES or NO action. Default: NO."
)
_database_parser.add_argument(
    "--maxsize",
    help="Max size for the predicted ORFs (in nucleotides)",
    default=1_000_000,
    type=int
)
_database_parser.add_argument(
    "--minsize",
    help="Min size for the predicted ORFs (in nucleotides)",
    default=30,
    type=int
)
_database_parser.add_argument(
    "--external_transcriptome",
    help="Fasta file containing the transcripts sequences from a "
    "transcriptome assembly that was done elsewhere. Also implies "
    "--external_gtf.",
    type=_types.FilePath
)
_database_parser.add_argument(
    "--external_gtf",
    help="If --external_transcriptome is provided, this argument should "
    "contain the path to a GTF file originated from the same "
    "transcriptome assembly",
    type=_types.FilePath
)

# =======
# MS MODE
# =======
# TODO: I have to check each of these args to add validation
# The problem is they are from another CLI application
_ms_parser = _modes.add_parser('ms')
_ms_parser.add_argument(
    '-o', '--outdir',
    help='Inform output directory',
    type=_types.DirectoryName,
    required=True
)
_ms_parser.add_argument(
    "--Transcriptome", "-T",
    dest='transcriptome',
    action=_types.YesOrNoBooleanAction,
    help="Whether transcriptome database was generated or not. If the "
    "transcriptome database was not generated, ignore this."
)
_ms_parser.add_argument(
    "--Mass_spec", "-M",
    dest='mass_spec',
    help="Inform the directory containing all your .mzML files.",
    required=True,
    type=os.path.abspath
)
_ms_parser.add_argument(
    "--t",
    help="The precursor mass tolerance. (e.g. 2.5Da, 20ppm or "
    "0.5Da,2.5Da; Default: 20ppm) Use a comma to define asymmetric "
    "values. E.g. --t 0.5Da,2.5Da will set 0.5Da to the left (ObsMass < "
    "TheoMass) and 2.5Da to the right (ObsMass > TheoMass)"
)
_ms_parser.add_argument(
    "--ti",
    help="The isotope error range. (Range of allowed isotope peak errors; "
    "Default: 0,1). Takes into account the error introduced by choosing a "
    "non-monoisotopic peak for fragmentation. The combination of --t and "
    "--ti determines the precursor mass tolerance. E.g. --t 20ppm --ti "
    "-1,2 tests abs(ObservedPepMass - TheoreticalPepMass - n * 1.00335Da) "
    "< 20ppm for n = -1, 0, 1, 2."
)
_ms_parser.add_argument(
    "--thread",
    help="Number of concurrent threads to be executed; Default: Number of "
    "available cores. This is best set to the number of physical cores in "
    "a single NUMA node. Generally a single NUMA node is 1 physical "
    "processor. The default will try to use hyperthreading cores, which "
    "can increase the amount of time this process will take. This is "
    "because the part of Scoring param generation that is multithreaded "
    "is also I/O intensive."
)
_ms_parser.add_argument(
    "--tasks",
    help="(Override the number of tasks to use on the threads; Default: "
    "(internally calculated based on inputs)) More tasks than threads "
    "will reduce the memory requirements of the search, but will be "
    "slower (how much depends on the inputs). 1 <= tasks <= numThreads: "
    "will create one task per thread, which is the original behavior. "
    "tasks = 0: use default calculation - minimum of: (threads*3) and "
    "(numSpectra/250). tasks < 0: multiply number of threads by "
    "abs(tasks) to determine number of tasks (i.e., -2 means 2 * "
    "numThreads tasks). One task per thread will use the most memory, but "
    "will usually finish the fastest. 2-3 tasks per thread will use "
    "comparably less memory, but may cause the search to take 1.5 to 2 "
    "times as long."
)
_ms_parser.add_argument(
    "--m",
    help="[-m FragmentationMethodID] (Fragmentation Method, Default: 0) 0 "
    "means as written in the spectrum or CID if no info (Default) 1 means "
    "CID, 2 means ETD, 3 means HCD"
)
_ms_parser.add_argument(
    "--inst",
    help="The instrument used for the experiment. (0: Low-res LCQ/LTQ "
    "(Default), 1: Orbitrap/FTICR/Lumos, 2: TOF, 3: Q-Exactive)"
)
_ms_parser.add_argument(
    "--e",
    help="The enzyme used for protein digestion. (0: unspecific cleavage, "
    "1: Trypsin (Default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: "
    "glutamyl endopeptidase, 6: Arg-C, 7: Asp-N, 8: alphaLP, 9: no "
    "cleavage)"
)
_ms_parser.add_argument(
    "--protocol",
    help="0: Automatic (Default), 1: Phosphorylation, 2: iTRAQ, 3: "
    "iTRAQPhospho, 4: TMT, 5: Standard"
)
_ms_parser.add_argument(
    "--ntt",
    help="Number of Tolerable Termini, Default: 2) E.g. For trypsin, 0: "
    "non-tryptic, 1: semi-tryptic, 2: fully-tryptic peptides only."
)
_ms_parser.add_argument(
    "--mod",
    help="Modification file; Default: standard amino acids with fixed "
    "C+57; only if -mod is not specified)"
)
_ms_parser.add_argument(
    "--minLength",
    help="Minimum peptide length to consider; Default: 6"
)
_ms_parser.add_argument(
    "--maxLength",
    help="Maximum peptide length to consider; Default: 40"
)
_ms_parser.add_argument(
    "--minCharge",
    help="Minimum precursor charge to consider if charges are not "
    "specified in the spectrum file; Default: 2"
)
_ms_parser.add_argument(
    "--maxCharge",
    help="Maximum precursor charge to consider if charges are not "
    "specified in the spectrum file; Default: 3"
)
_ms_parser.add_argument(
    "--n",
    help="Number of matches per spectrum to be reported; Default: 1"
)
_ms_parser.add_argument(
    "--ccm",
    help="Mass of charge carrier; Default: mass of proton (1.00727649)"
)
_ms_parser.add_argument(
    "--maxMissedCleavages",
    help="Exclude peptides with more than this number of missed cleavages "
    "from the search; Default: -1 (no limit)"
)
_ms_parser.add_argument(
    "--numMods",
    help="Maximum number of dynamic (variable) modifications per peptide; "
    "Default: 3"
)

# ===========
# POSTMS MODE
# ===========
_postms_parser = _modes.add_parser('postms')
_postms_parser.add_argument(
    '-o', '--outdir',
    help='Inform output directory',
    type=_types.DirectoryName,
    required=True
)
_postms_parser.add_argument(
    "--Mass_spec", "-M",
    dest='mass_spec',
    help="Inform the directory informed in the Peptide search step",
    type=_types.DirectoryPath
)
_postms_parser.add_argument(
    "--Transcriptome", "-T",
    dest='transcriptome',
    action=_types.YesOrNoBooleanAction,
    help="Whether transcriptome database was generated or not. If the "
    "transcriptome database was not generated, ignore this."
)
_postms_parser.add_argument(
    "--proteome", "-p",
    help="Proteome fasta file, i.e. 'Uniprot.fasta'",
    required=True,
    type=_types.FilePath
)
_postms_parser.add_argument(
    "--genome", "-g",
    help="Path to the reference genome.",
    required=True,
    type=_types.FilePath
)
_postms_parser.add_argument(
    "--genbank",
    help="Path to the genbank file. If this is abscent, search for ORFs "
    "inside ncRNAs is not going to be performed.",
    type=_types.FilePath
)
_postms_parser.add_argument(
    "--gff",
    help="Path to a GFF file",
    type=_types.FilePath
)
# TODO: Check appropiate type
_postms_parser.add_argument(
    "--pep",
    help="Posterior Error Probability (PEP) cutoff for percolator output.",
    default=1
)
_postms_parser.add_argument(
    "--qvalue",
    help="Q-Value cutoff for percolator output",
    type=float,
    default=0.01
)
_postms_parser.add_argument(
    "--starts",
    help="The start codons that should be used for the six and three "
    "frame translation. By default, the NCBI bacterial genetic code is "
    "used. Space separated list.",
    type=_types.CommaList(_types.Codon),
    default='TTG,CTG,ATT,ATC,ATA,ATG,GTG'
    )
_postms_parser.add_argument(
    "--stops",
    help="The stop codons that should be used for the six and three frame "
    "translation. Space separated list.",
    type=_types.CommaList(_types.Codon),
    default='TAA,TAG,TGA'
)
_postms_parser.add_argument(
    "--rrna",
    help="Path to the fasta file containing the sequence for the 16S rRNA",
    type=_types.FilePath
)
_postms_parser.add_argument(
    "--maxsize",
    help="Maximum ORF size (in nucleotides)",
    type=int,
    default=300
)

# =============
# VALIDATE MODE
# =============
_validate_parser = _modes.add_parser('validate')
_validate_parser.add_argument(
    '-o', '--outdir',
    help='Inform output directory',
    type=_types.DirectoryName,
    required=True
)
_validate_parser.add_argument(
    "--Transcriptome", "-T",
    dest='transcriptome',
    action=_types.YesOrNoBooleanAction,
    help="Whether transcriptome database was generated or not. If the "
    "transcriptome database was not generated, ignore this."
)

# ============
# TESTING MODE
# ============
_testing_parser = _modes.add_parser('testing')
_testing_parser.add_argument(
    '-o', '--outdir',
    help='Inform output directory',
    type=_types.DirectoryName,
    required=True
)
_testing_parser.add_argument(
    "--skip_assembly",
    help="Write TRUE if you want to skip the assembly checking.",
    default="FALSE"
)
_testing_parser.add_argument(
    "--skip_db",
    help="Write TRUE if you want to skip the database checking.",
    default="FALSE"
)
_testing_parser.add_argument(
    "--skip_ms",
    help="Write TRUE if you want to skip the peptide search checking.",
    default="FALSE"
)
_testing_parser.add_argument(
    "--skip_postms",
    help="Write TRUE if you want to skip the results processingchecking.",
    default="FALSE"
)
_testing_parser.add_argument(
    "--skip_validation",
    help="Write TRUE if you want to skip the validation checking",
    default="FALSE"
)
_testing_parser.add_argument(
    "--threads",
    help="Number of threads for samtools and stringtie to use.",
    type=_types.PositiveInt
)
_testing_parser.add_argument(
    '--memory',
    help="Maximum memory per thread for samtools sort to use.",
    type=_types.Memory
)

# ============
# METRICS MODE
# ============
_metrics_parser = _modes.add_parser('metrics')
_metrics_parser.add_argument(
    '-o', '--outdir',
    help='Inform output directory',
    type=_types.DirectoryName,
    required=True
)
_metrics_parser.add_argument("--non-redundant", action="store_true")
_metrics_parser.add_argument("--maxsize", type=int, default=100)


def get_parsers() -> tuple[
    argparse.ArgumentParser,
    dict[str, argparse.ArgumentParser],
]:
    """Get the main parser and a dict of the subparsers for the pipeline.

    The key to each subparser is its respective mode, stored in the `mode`
    variable of the args :class:`Namespace`.

    Returns
    -------
    tuple[argparse.ArgumentParser, list[argparse.ArgumentParser]]
        The first item is the main parser, while the second item is a dict
        of `mode: subparser`.
    """
    subparsers = {
        'assembly': _assembly_parser,
        'database': _database_parser,
        'ms': _ms_parser,
        'postms': _postms_parser,
        'validate': _validate_parser,
        'metrics': _metrics_parser,
    }

    return _parser, subparsers


def exit(
    status: int = 0,
    message: t.Optional[str] = None,
    end: t.Optional[str] = "\n"
) -> t.NoReturn:
    if end is not None and message is not None:
        message += end
    _parser.exit(status, message)


def error(message: str, end: t.Optional[str] = "\n") -> t.NoReturn:
    """Print to stderr, then exit with return value 2."""
    if end is not None:
        message += end
    _parser.error(message)


def stderr(message: str, end: t.Optional[str] = "\n") -> None:
    """Print to stderr."""
    print(message, file=sys.stderr, end=end)
