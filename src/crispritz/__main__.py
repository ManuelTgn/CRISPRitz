"""
CRISPRitz {version}

Copyright (C) 2026 Pinellolab <lpinello@mgh.harvard.edu>

CRISPRitz: High-Throughput and Variant-Aware In Silico Off-Target Sites Identification 
For CRISPR Genome Editing

CRISPRitz is a software package containing five tools to perform variant-aware 
off-target prediction and result assessement on CRISPR/Cas experiments.

Usage:
    crispritz add-variants <options>

Run 'crispritz -h/--help' to display the complete help
"""

from .utils import TOOLNAME
from .crispritz_argparse import CrispritzArgumentParser
from .version import __version__
from .exception_handlers import sigint_handler

from time import time

import sys

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
    return parser


def main():
    start = time()  # track eleapsed time
    try:
        parser = create_parser_crispritz()  # parse input argument using custom parser
        if not sys.argv[1:]:  # no input args -> print help and exit
            parser.error_noargs()
        args = parser.parse_args(sys.argv[1:])  # parse input args
    except KeyboardInterrupt:
        sigint_handler()  # catch SIGINT and exit gracefully
    sys.stdout.write(f"{TOOLNAME} - Elapsed time {time() - start:.2f}s\n")



# --------------------------------> ENTRY POINT <--------------------------------
if __name__ == "__main__":
    main()

