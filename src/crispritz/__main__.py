"""
CRISPRitz {version}

Copyright (C) 2026 Pinello Lab & InfOmics Lab <lpinello@mgh.harvard.edu> <rosalba.giugno@univr.it>

CRISPRitz: High-Throughput and Variant-Aware In Silico Off-Target Sites Identification
For CRISPR Genome Editing

CRISPRitz is a software package containing five tools to perform variant-aware
off-target prediction and result assessement on CRISPR/Cas experiments.

Usage:
    crispritz add-variants <options>

Run 'crispritz -h/--help' to display the complete help
"""

from .utils import TOOLNAME, SUBCOMMANDS
from .crispritz_argparse import CrispritzArgumentParser, CrispritzEnrichmentInputArgs
from .version import __version__
from .exception_handlers import sigint_handler
from .enrichment import add_variants_cli
from .ternary_search_tree import index_genome_cli

from argparse import _SubParsersAction
from time import time

import sys
import os


def create_parser_crispritz() -> CrispritzArgumentParser:
    # force displaying docstring at each usage display and force
    # the default help to not being shown
    parser = CrispritzArgumentParser(usage=__doc__, add_help=False)  # type: ignore
    group = parser.add_argument_group("Options")  # arguments group
    # add help and version arguments
    group.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit"
    )
    group.add_argument(
        "--version",
        action="version",
        help=f"Show {TOOLNAME} version and exit",
        version=__version__,
    )
    # create subparsers for different functionalities
    subparsers = parser.add_subparsers(
        dest="command",
        title="Available commands",
        metavar="",  # needed for help formatting (avoid <command to be displayed>)
        description=None,
    )
    # crispritz enrichment
    create_enrichment_parser(subparsers)
    create_indexing_parser(subparsers)
    return parser


def create_enrichment_parser(subparser: _SubParsersAction) -> _SubParsersAction:
    """Create and configure the argument parser for the enrichment subcommand.

    Defines required and optional arguments for running the genome enrichment
    pipeline, including input directories, enrichment options, and runtime
    settings.

    Args:
        subparser (_SubParsersAction): The subparsers collection to which the
            enrichment parser will be added.

    Returns:
        _SubParsersAction: The configured enrichment subparser.
    """
    parser_enrichment = subparser.add_parser(
        SUBCOMMANDS[0],
        usage="CRISPRitz add-variants {version}\n\nUsage:\n"
        "\tcrisprhawk add-variants --vcf <vcf> --genome <fasta>\n\n",
        description="Genome enrichment pipeline: parses input VCF files to "
        "integrate sequence variants (SNPs and indels) into the reference FASTA "
        "files",
        help="Genome enrichment pipeline that processes input VCF files and adds "
        "sequence variants (SNPs and indels) to the corresponding reference FASTA "
        "files. For each input FASTA, an enriched FASTA is generated. SNPs are "
        "represented using IUPAC ambiguity codes to encode both reference and "
        "alternative alleles. When enabled, indels are incorporated individually "
        "by inserting or removing bases in the regions surrounding the variant "
        "position",
        add_help=False,
    )
    general_group = parser_enrichment.add_argument_group("General options")
    general_group.add_argument(
        "-h", "--help", action="help", help="show this help message and exit"
    )
    required_group = parser_enrichment.add_argument_group("Options")
    required_group.add_argument(
        "--vcf",
        type=str,
        metavar="VCF-DIR",
        dest="vcf",
        required=True,
        help="directory containing the VCF files used for genome enrichment. "
        "Each chromosome must be stored in a separate VCF file "
        "(e.g., chr1.vcf.gz, chr2.vcf.gz)",
    )
    required_group.add_argument(
        "--genome",
        type=str,
        metavar="FASTA-DIR",
        dest="genome",
        required=True,
        help="directory containing the reference genome FASTA files. "
        "Each chromosome must be stored in a separate FASTA file "
        "(e.g., chr1.fa, chr2.fa). All FASTA files in this directory "
        "will be used as the reference genome",
    )
    optional_group = parser_enrichment.add_argument_group("Optional arguments")
    optional_group.add_argument(
        "--indels",
        action="store_true",
        dest="indels",
        default=False,
        help="include indels during genome enrichment. "
        "If enabled, insertions and deletions are applied to the reference "
        "sequence individually (default: disabled)",
    )
    optional_group.add_argument(
        "--keep",
        action="store_true",
        dest="keep",
        default=False,
        help="include all variants during genome enrichment, regardless of their "
        "FILTER status. By default, only variants with FILTER=PASS are included "
        "(default: disabled)",
    )
    optional_group.add_argument(
        "--outdir",
        type=str,
        metavar="OUTDIR",
        dest="outdir",
        required=False,
        default=os.getcwd(),
        help="directory where output folder will be written. "
        "(default: a `variants_genome` folder will be created in the current "
        "working directory)",
    )
    optional_group.add_argument(
        "--threads",
        type=int,
        metavar="THREADS",
        dest="threads",
        required=False,
        default=1,
        help="number of threads. Use 0 for using all available cores (default: 1)",
    )
    optional_group.add_argument(
        "--verbosity",
        type=int,
        metavar="VERBOSITY",
        dest="verbosity",
        required=False,
        default=1,  # minimal output
        help="verbosity level of output messages: 0 = Silent, 1 = Normal, 2 = "
        "Verbose, 3 = Debug (default: 1)",
    )
    optional_group.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="enter debug mode and trace the full error stack",
    )
    return parser_enrichment

def create_indexing_parser(subparser: _SubParsersAction) -> _SubParsersAction:
    parser_enrichment = subparser.add_parser(
        SUBCOMMANDS[1],
        usage="CRISPRitz index-genome {version}\n\nUsage:\n"
        "\tcrisprhawk index-genome\n\n",
        description="",
        help="",
        add_help=False,
    )
    general_group = parser_enrichment.add_argument_group("General options")
    general_group.add_argument(
        "-h", "--help", action="help", help="show this help message and exit"
    )
    # required_group = parser_enrichment.add_argument_group("Options")
    # required_group.add_argument(
    #     "--vcf",
    #     type=str,
    #     metavar="VCF-DIR",
    #     dest="vcf",
    #     required=True,
    #     help="directory containing the VCF files used for genome enrichment. "
    #     "Each chromosome must be stored in a separate VCF file "
    #     "(e.g., chr1.vcf.gz, chr2.vcf.gz)",
    # )
    # required_group.add_argument(
    #     "--genome",
    #     type=str,
    #     metavar="FASTA-DIR",
    #     dest="genome",
    #     required=True,
    #     help="directory containing the reference genome FASTA files. "
    #     "Each chromosome must be stored in a separate FASTA file "
    #     "(e.g., chr1.fa, chr2.fa). All FASTA files in this directory "
    #     "will be used as the reference genome",
    # )
    # optional_group = parser_enrichment.add_argument_group("Optional arguments")
    # optional_group.add_argument(
    #     "--indels",
    #     action="store_true",
    #     dest="indels",
    #     default=False,
    #     help="include indels during genome enrichment. "
    #     "If enabled, insertions and deletions are applied to the reference "
    #     "sequence individually (default: disabled)",
    # )
    # optional_group.add_argument(
    #     "--keep",
    #     action="store_true",
    #     dest="keep",
    #     default=False,
    #     help="include all variants during genome enrichment, regardless of their "
    #     "FILTER status. By default, only variants with FILTER=PASS are included "
    #     "(default: disabled)",
    # )
    # optional_group.add_argument(
    #     "--outdir",
    #     type=str,
    #     metavar="OUTDIR",
    #     dest="outdir",
    #     required=False,
    #     default=os.getcwd(),
    #     help="directory where output folder will be written. "
    #     "(default: a `variants_genome` folder will be created in the current "
    #     "working directory)",
    # )
    # optional_group.add_argument(
    #     "--threads",
    #     type=int,
    #     metavar="THREADS",
    #     dest="threads",
    #     required=False,
    #     default=1,
    #     help="number of threads. Use 0 for using all available cores (default: 1)",
    # )
    # optional_group.add_argument(
    #     "--verbosity",
    #     type=int,
    #     metavar="VERBOSITY",
    #     dest="verbosity",
    #     required=False,
    #     default=1,  # minimal output
    #     help="verbosity level of output messages: 0 = Silent, 1 = Normal, 2 = "
    #     "Verbose, 3 = Debug (default: 1)",
    # )
    # optional_group.add_argument(
    #     "--debug",
    #     action="store_true",
    #     default=False,
    #     help="enter debug mode and trace the full error stack",
    # )
    return parser_enrichment



def main():
    start = time()  # track eleapsed time
    try:
        parser = create_parser_crispritz()  # parse input argument using custom parser
        if not sys.argv[1:]:  # no input args -> print help and exit
            parser.error_noargs()
        args = parser.parse_args(sys.argv[1:])  # parse input args
        if args.command == SUBCOMMANDS[0]:  # add-variants
            add_variants_cli(CrispritzEnrichmentInputArgs(args, parser))
        if args.command == SUBCOMMANDS[1]:  # index-genome
            index_genome_cli()
    except KeyboardInterrupt:
        sigint_handler()  # catch SIGINT and exit gracefully
    sys.stdout.write(f"{TOOLNAME} - Elapsed time {time() - start:.2f}s\n")


# --------------------------------> ENTRY POINT <--------------------------------
if __name__ == "__main__":
    main()
