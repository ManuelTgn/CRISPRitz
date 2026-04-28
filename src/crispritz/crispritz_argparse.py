""" """

from .utils import COMMAND, VERBOSITYLVL
from .version import __version__

from argparse import (
    SUPPRESS,
    ArgumentParser,
    HelpFormatter,
    Action,
    _MutuallyExclusiveGroup,
    Namespace,
)
from typing import Iterable, Optional, TypeVar, Tuple, Dict, NoReturn, List
from colorama import Fore
from glob import glob

import multiprocessing
import sys
import os

# define abstract generic types for typing
_D = TypeVar("_D")
_V = TypeVar("_V")


class CrispritzArgumentParser(ArgumentParser):

    class CrispritzHelpFormatter(HelpFormatter):

        def add_usage(  # type: ignore
            self,
            usage: str,
            actions: Iterable[Action],
            groups: Iterable[_MutuallyExclusiveGroup],
            prefix: Optional[str] = None,
        ) -> None:
            # add usage description for help only if the set action is not to
            # suppress the display of the help formatter
            if usage != SUPPRESS:
                args = (usage, actions, groups, "")
                self._add_item(self._format_usage, args)  # initialize the formatter

    def __init__(self, *args: Tuple[_D], **kwargs: Dict[_D, _V]) -> None:
        # set custom help formatter defined as
        kwargs["formatter_class"] = self.CrispritzHelpFormatter  # type: ignore
        # replace the default version display in usage help with a custom
        # version display formatter
        if "usage" in kwargs:
            kwargs["usage"] = kwargs["usage"].replace("{version}", __version__)  # type: ignore
        # initialize argument parser object with input parameters for
        # usage display
        super().__init__(*args, **kwargs)  # type: ignore

    def error(self, error: str) -> NoReturn:  # type: ignore
        # display error messages raised by argparse in red
        errormsg = (
            f"{Fore.RED}\nERROR: {error}.{Fore.RESET}"
            + f"\n\nRun {COMMAND} -h for usage\n\n"
        )
        sys.stderr.write(errormsg)  # write error to stderr
        sys.exit(os.EX_USAGE)  # exit execution -> usage error

    def error_noargs(self) -> None:
        self.print_help()  # if no input argument, print help
        sys.exit(os.EX_NOINPUT)  # exit with no input code


class CrispritzInputArgs:
    def __init__(self, args: Namespace, parser: CrispritzArgumentParser) -> None:
        self._args = args
        self._parser = parser

    def _validate_output_folder(self) -> None:
        parentdir = os.path.dirname(self._args.outdir)
        _check_folder(
            parentdir, self._parser, f"Cannot create output folder {self._args.outdir}"
        )
        outdir = os.path.abspath(self._args.outdir)
        if not os.path.exists(outdir) or not os.path.isdir(outdir):
            os.makedirs(outdir)  # create output folder
        self._outdir = outdir
        assert os.path.isdir(self._outdir)

    def _validate_threads(self) -> None:
        self._threads = _validate_threads_num(self._args.threads, self._parser)

    def _validate_verbosity(self) -> None:
        self._verbosity = _validate_verbosity_value(self._args.verbosity, self._parser)

    @property
    def outdir(self) -> str:
        return self._args.outdir

    @property
    def threads(self) -> int:
        return self._args.threads

    @property
    def verbosity(self) -> int:
        return self._args.verbosity

    @property
    def debug(self) -> bool:
        return self._args.debug


class CrispritzEnrichmentInputArgs(CrispritzInputArgs):

    def __init__(self, args: Namespace, parser: CrispritzArgumentParser) -> None:
        super().__init__(args, parser)
        self._check_consistency()

    def _validate_vcf_folder(self) -> None:
        _check_folder(
            self._args.vcf, self._parser, f"Cannot find VCF folder {self._args.vcf}"
        )
        self._vcfs = glob(os.path.join(self._args.vcf, "*.vcf.gz"))
        _check_retrieved_files(
            self._vcfs, self._parser, f"No VCF file found in {self._args.vcf}"
        )

    def _validate_genome_folder(self) -> None:
        _check_folder(
            self._args.genome,
            self._parser,
            f"Cannot find input genome folder {self._args.genome}",
        )
        self._fastas = glob(os.path.join(self._args.genome, "*.fa")) + glob(
            os.path.join(self._args.genome, "*.fasta")
        )
        _check_retrieved_files(
            self._fastas, self._parser, f"No FASTA file found in {self._args.genome}"
        )

    def _check_consistency(self) -> None:
        self._validate_vcf_folder()  # check vcf folder
        self._validate_genome_folder()  # check genome folder
        self._validate_output_folder()  # check output folder
        self._validate_threads()  # check threads number
        self._validate_verbosity()  # check verbosity

    @property
    def vcfs(self) -> List[str]:
        return self._vcfs

    @property
    def fastas(self) -> List[str]:
        return self._fastas

    @property
    def indels(self) -> bool:
        return self._args.indels

    @property
    def keep(self) -> bool:
        return self._args.keep


class CrispritzIndexingInputArgs(CrispritzInputArgs):

    def __init__(self, args: Namespace, parser: CrispritzArgumentParser) -> None:
        super().__init__(args, parser)
        self._check_consistency()  # check input args consistency

    def _validate_genome_folder(self) -> None:
        _check_folder(
            self._args.genome,
            self._parser,
            f"Cannot find input genome folder {self._args.genome}",
        )
        self._fastas = glob(os.path.join(self._args.genome, "*.fa")) + glob(
            os.path.join(self._args.genome, "*.fasta")
        )
        _check_retrieved_files(
            self._fastas, self._parser, f"No FASTA file found in {self._args.genome}"
        )

    def _validate_pam_file(self) -> None:
        _check_file(
            self._args.pam_file,
            self._parser,
            f"Cannot find input PAM file {self._args.pam_file}",
        )

    def _validate_bmax(self) -> None:
        if self._args.bmax < 0:
            self._parser.error(f"Invalid max bulge value: {self._args.bmax}")

    def _check_consistency(self) -> None:
        self._validate_genome_folder()  # check genome folder
        self._validate_pam_file()  # check pam file
        self._validate_bmax()  # check max bulge
        self._validate_output_folder()  # check output folder
        self._validate_threads()  # check threads number
        self._validate_verbosity()  # check verbosity

    @property
    def fastas(self) -> List[str]:
        return self._fastas

    @property
    def pam_file(self) -> str:
        return self._args.pam_file

    @property
    def bmax(self) -> int:
        return self._args.bmax


class CrispritzSearchInputArgs(CrispritzInputArgs):

    def __init__(self, args: Namespace, parser: CrispritzArgumentParser) -> None:
        super().__init__(args, parser)
        self._check_consistency()

    def _check_consistency(self) -> None:
        self._validate_index_genome()  # check genome index folder
        self._validate_pam_file()  # check pam file
        self._validate_guides_file()  # check guides file
        self._validate_mm()  # check mismatch
        self._validate_bdna()  # check dna bulge
        self._validate_brna()  # check rna bulge
        self._validate_output_folder()  # check output folder
        self._validate_threads()  # check threads number
        self._validate_verbosity()  # check verbosity

    def _validate_index_genome(self) -> None:
        _check_folder(
            self._args.genome_index,
            self._parser,
            f"Cannot find genome index folder {self._args.genome_index}",
        )
        self._indexes = glob(os.path.join(self._args.genome_index, "*.bin"))
        _check_retrieved_files(
            self._indexes,
            self._parser,
            f"No TST index found in {self._args.genome_index}",
        )

    def _validate_pam_file(self) -> None:
        _check_file(
            self._args.pam_file,
            self._parser,
            f"Cannot find input PAM file {self._args.pam_file}",
        )

    def _validate_guides_file(self) -> None:
        _check_file(
            self._args.guides_file,
            self._parser,
            f"Cannot find input guides file {self._args.guides_file}",
        )

    def _validate_mm(self) -> None:
        if self._args.mm < 0:
            self._parser.error(f"Invalid max mismatch value: {self._args.mm}")

    def _validate_bdna(self) -> None:
        if self._args.bdna < 0:
            self._parser.error(f"Invalid max DNA bulge value: {self._args.bdna}")

    def _validate_brna(self) -> None:
        if self._args.brna < 0:
            self._parser.error(f"Invalid max RNA bulge value: {self._args.brna}")

    @property
    def indexes(self) -> List[str]:
        return self._indexes

    @property
    def pam_file(self) -> str:
        return self._args.pam_file

    @property
    def guides_file(self) -> str:
        return self._args.guides_file

    @property
    def mm(self) -> int:
        return self._args.mm

    @property
    def bdna(self) -> int:
        return self._args.bdna

    @property
    def brna(self) -> int:
        return self._args.brna


def _check_folder(dirname: str, parser: CrispritzArgumentParser, msg: str) -> None:
    if not os.path.exists(dirname) or not os.path.isdir(dirname):
        parser.error(msg)


def _check_file(fname: str, parser: CrispritzArgumentParser, msg: str) -> None:
    if not os.path.exists(fname) or not os.path.isfile(fname):
        parser.error(msg)


def _check_retrieved_files(
    fnames: List[str], parser: CrispritzArgumentParser, msg: str
) -> None:
    if not fnames:
        parser.error(msg)


def _validate_threads_num(threads: int, parser: CrispritzArgumentParser) -> int:
    max_threads = multiprocessing.cpu_count()
    if threads < 0 or threads > max_threads:
        parser.error(
            f"Forbidden number of threads provided ({threads}). "
            f"Max number of available cores: {max_threads}"
        )
    return max_threads if threads == 0 else threads


def _validate_verbosity_value(verbosity: int, parser: CrispritzArgumentParser) -> int:
    if verbosity not in VERBOSITYLVL:
        parser.error(f"Forbidden verbosity level selected ({verbosity})")
    return verbosity
