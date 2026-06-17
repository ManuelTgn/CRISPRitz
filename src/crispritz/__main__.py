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

from argparse import _SubParsersAction
from time import time

import os
import sys


from .crispritz_argparse import CrispritzArgumentParser
from .crispritz_inputargs import (
    CrispritzEnrichmentInputArgs,
    CrispritzIndexingInputArgs,
    CrispritzSearchInputArgs,
    OUTPUT_MODES,
)
from .enrichment import add_variants_cli
from .exception_handlers import sigint_handler
from .indexing import index_genome_cli
from .search import search_offtargets_cli
from .utils import TOOLNAME, SUBCOMMANDS
from .version import __version__


def _create_parser_crispritz() -> CrispritzArgumentParser:
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
    _create_enrichment_parser(subparsers)
    _create_indexing_parser(subparsers)
    _create_search_parser(subparsers)
    return parser


def _create_enrichment_parser(subparser: _SubParsersAction) -> _SubParsersAction:
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
        "\tcrispritz add-variants --vcf <vcf> --genome <fasta>\n\n",
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


def _create_indexing_parser(subparser: _SubParsersAction) -> _SubParsersAction:
    parser_indexing = subparser.add_parser(
        SUBCOMMANDS[1],
        usage="CRISPRitz index-genome {version}\n\nUsage:\n"
        "\tcrispritz index-genome --genome <genome-dir> --genome-name "
        "<genome-name> --pam <pam-file>\n\n",
        description="Create a genome TST index for fast (optionally bulge-aware) "
        "off-target candidate searches.",
        help="Build a genome index (TST: Ternary Search Tree) used for fast "
        "off-target candidate retrieval. This command scans the input FASTA "
        "files, extracts all candidate targets matching the provided PAM, and "
        "stores them in a compact TST-based index. The resulting index enables "
        "rapid searches and supports bulge-aware matching (DNA/RNA bulges)",
        add_help=False,
    )
    general_group = parser_indexing.add_argument_group("General options")
    general_group.add_argument(
        "-h", "--help", action="help", help="show this help message and exit"
    )
    required_group = parser_indexing.add_argument_group("Options")
    required_group.add_argument(
        "--genome",
        type=str,
        metavar="FASTA-DIR",
        dest="genome",
        required=True,
        help="path to a directory containing the reference/enriched genome in "
        "FASTA format. Each chromosome must be provided as a separate FASTA file "
        "(e.g., chr1.fa, chr2.fa, chrX.fa). All FASTA files in this directory "
        "will be scanned to extract target candidates and to build the TST index",
    )
    required_group.add_argument(
        "--genome-name",
        type=str,
        metavar="GENOME-NAME",
        dest="genome_name",
        required=True,
        help="identifier used to name the generated TST index. A folder with "
        "this name will be created to store the index associated with the input "
        "genome",
    )
    required_group.add_argument(
        "--pam",
        type=str,
        metavar="PAM-FILE",
        dest="pam_file",
        required=True,
        help="path to a text file specifying the PAM model. The file must "
        "contain: (1) the full pattern including a number of 'N' characters "
        "equal to the guide length, and (2) a space-separated integer indicating "
        "the PAM length. Example format: NNNNNNNNNNNNNNNNNNNNGG 3",
    )
    optional_group = parser_indexing.add_argument_group("Optional arguments")
    optional_group.add_argument(
        "--bmax",
        type=int,
        metavar="BMAX",
        dest="bmax",
        default=0,
        help="maximum number of bulges allowed during index construction and "
        "off-target search. Larger bulges increase search sensitivity but also "
        "computational cost (default: 0)",
    )
    optional_group.add_argument(
        "--outdir",
        type=str,
        metavar="OUTDIR",
        dest="outdir",
        required=False,
        default=os.getcwd(),
        help="directory where output files will be written. If not specified, a "
        "folder named '<GENOME-NAME>_<PAM>_<BMAX>' will be created in the "
        "current working directory.",
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
    return parser_indexing


def _create_search_parser(subparser: _SubParsersAction) -> _SubParsersAction:
    """Create and configure the argument parser for the search subcommand.

    Defines required and optional arguments for the off-target search pipeline:
    the pre-computed TST genome index, the PAM model, the guide RNAs, the edit
    budget (mismatches plus optional DNA/RNA bulges), and the output controls
    (which files to produce and in what serialization format).

    Args:
        subparser (_SubParsersAction): The subparsers collection to which the
            search parser will be added.

    Returns:
        _SubParsersAction: The configured search subparser.
    """
    parser_search = subparser.add_parser(
        SUBCOMMANDS[2],
        usage="CRISPRitz search {version}\n\nUsage:\n"
        "\tcrispritz search --index-genome <index-genome-dir> --pam <pam-file> "
        "--guides <guides-file> --mm <mismatches> [--bdna <n>] [--brna <n>] "
        "[--output-mode <targets|profile|both>] "
        "[--output-format <tsv|targets>] [--outdir <dir>] [--threads <n>]\n\n",
        description="Search for candidate CRISPR off-target sites for each input "
        "guide RNA against a pre-computed genome index based on a Ternary Search "
        "Tree (TST). The index already encodes the PAM and guide geometry fixed "
        "during the 'index-genome' step; this command supplies only the search "
        "tolerances and the desired output.",
        help="Identify candidate off-target sites within a user-defined edit "
        "distance, given as the number of mismatches (required) plus optional "
        "DNA and RNA bulges. The search traverses a pre-computed TST genome "
        "index, enumerating every site whose alignment to a guide stays within "
        "the edit budget. Results are reported as a tab-separated targets table "
        "and/or per-guide mismatch profiles, ready for downstream scoring and "
        "annotation.",
        add_help=False,
    )
    general_group = parser_search.add_argument_group("General options")
    general_group.add_argument(
        "-h", "--help", action="help", help="show this help message and exit"
    )
    required_group = parser_search.add_argument_group("Options")
    required_group.add_argument(
        "--index-genome",
        type=str,
        metavar="INDEX-GENOME",
        dest="genome_index",
        required=True,
        help="path to the genome index directory produced by the 'index-genome' "
        "command. It must contain the pre-computed TST partitions (*.bin) used "
        "for the off-target search",
    )
    required_group.add_argument(
        "--pam",
        type=str,
        metavar="PAM-FILE",
        dest="pam_file",  # fixed: was "pamfile"; CrispritzSearchInputArgs reads pam_file
        required=True,
        help="path to a text file specifying the PAM model. The file must "
        "contain: (1) the full pattern, with a run of 'N' characters equal to "
        "the guide length, and (2) a space-separated integer giving the PAM "
        "length. Must match the PAM used to build the index. "
        "Example: NNNNNNNNNNNNNNNNNNNNGG 3",
    )
    required_group.add_argument(
        "--guides",
        type=str,
        metavar="GUIDES-FILE",
        dest="guides_file",
        required=True,
        help="path to a text file with one guide RNA per line. Each guide must "
        "match the length and PAM structure implied by the index. "
        "Example: CTAACAGTTGCTTTTATCACNNN",
    )
    required_group.add_argument(
        "--mm",
        type=int,
        metavar="MISMATCHES",
        dest="mm",
        required=True,
        help="maximum number of mismatches (substitutions) allowed between a "
        "guide RNA and a candidate off-target site",
    )
    optional_group = parser_search.add_argument_group("Optional arguments")
    optional_group.add_argument(
        "--bdna",
        type=int,
        metavar="DNA-BULGES",
        dest="bdna",
        required=False,
        default=0,
        help="maximum number of DNA bulges allowed in the guide/off-target "
        "alignment. Cannot exceed the bulge budget the index was built with "
        "(default: 0)",
    )
    optional_group.add_argument(
        "--brna",
        type=int,
        metavar="RNA-BULGES",
        dest="brna",
        required=False,
        default=0,
        help="maximum number of RNA bulges allowed in the guide/off-target "
        "alignment. Cannot exceed the bulge budget the index was built with "
        "(default: 0)",
    )
    optional_group.add_argument(
        "--output-mode",
        type=str,
        metavar="MODE",
        dest="output_mode",
        choices=OUTPUT_MODES,
        required=False,
        default=OUTPUT_MODES[0],  # "both"
        help="which result files to produce: 'targets' writes only the "
        "off-target sites table; 'profile' writes only the per-guide mismatch "
        "profile summaries; 'both' writes targets and profiles "
        "(choices: targets, profile, both; default: both)",
    )
    optional_group.add_argument(
        "--outdir",
        type=str,
        metavar="OUTDIR",
        dest="outdir",
        required=False,
        default=os.getcwd(),
        help="directory where output files will be written. If not specified, "
        "results are written to the current working directory",
    )
    optional_group.add_argument(
        "--threads",
        type=int,
        metavar="THREADS",
        dest="threads",
        required=False,
        default=1,
        help="number of threads. Use 0 for all available cores (default: 1)",
    )
    optional_group.add_argument(
        "--verbosity",
        type=int,
        metavar="VERBOSITY",
        dest="verbosity",
        required=False,
        default=1,
        help="verbosity level of output messages: 0 = Silent, 1 = Normal, "
        "2 = Verbose, 3 = Debug (default: 1)",
    )
    optional_group.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="enter debug mode and trace the full error stack",
    )
    return parser_search

def _parse_input_args():
    parser = _create_parser_crispritz()  # parse input argument using custom parser
    if not sys.argv[1:]:  # no input args -> print help and exit
        parser.error_noargs()
    args = parser.parse_args(sys.argv[1:])  # parse input args
    return parser, args


def main():
    start = time()  # track eleapsed time
    try:
        parser, args = _parse_input_args()  # parse input arguments
        if args.command == SUBCOMMANDS[0]:  # add-variants
            add_variants_cli(CrispritzEnrichmentInputArgs(args, parser))
        elif args.command == SUBCOMMANDS[1]:  # index-genome
            index_genome_cli(CrispritzIndexingInputArgs(args, parser))
        elif args.command == SUBCOMMANDS[2]:  # search
            search_offtargets_cli(CrispritzSearchInputArgs(args, parser))
    except KeyboardInterrupt:
        sigint_handler()  # catch SIGINT and exit gracefully
    sys.stdout.write(f"{TOOLNAME} - Elapsed time {time() - start:.2f}s\n")


# --------------------------------> ENTRY POINT <--------------------------------
if __name__ == "__main__":
    main()
