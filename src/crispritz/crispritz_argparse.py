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
    """Custom argument parser for CRISPRitz command-line interface.

    This class extends argparse.ArgumentParser to provide custom help formatting,
    error handling, and version display for the CRISPRitz tool.

    Attributes:
        usage (str): The usage string for the parser, with version information.
        formatter_class (type): The custom help formatter class.
    """

    class CrispritzHelpFormatter(HelpFormatter):
        """Custom help formatter for CRISPRitz argument parser.

        This formatter customizes the usage message display for the help output.

        Attributes:
            None
        """

        def add_usage(  # type: ignore
            self,
            usage: str,
            actions: Iterable[Action],
            groups: Iterable[_MutuallyExclusiveGroup],
            prefix: Optional[str] = None,
        ) -> None:
            """Add a usage message to the help output.

            Displays the usage description unless suppressed.

            Args:
                usage (str): The usage string to display.
                actions (Iterable[Action]): The actions associated with the parser.
                groups (Iterable[_MutuallyExclusiveGroup]): Mutually exclusive
                    groups.
                prefix (Optional[str]): Optional prefix for the usage message.
            """
            # add usage description for help only if the set action is not to
            # suppress the display of the help formatter
            if usage != SUPPRESS:
                args = (usage, actions, groups, "")
                self._add_item(self._format_usage, args)  # initialize the formatter

    def __init__(self, *args: Tuple[_D], **kwargs: Dict[_D, _V]) -> None:
        """Initialize the CRISPR-HAWK argument parser.

        Sets up the parser with a custom help formatter and version display.

        Args:
            *args: Positional arguments for ArgumentParser.
            **kwargs: Keyword arguments for ArgumentParser.
        """
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
        """Display an error message and exit.

        Shows the error in red and suggests running the help command.

        Args:
            error (str): The error message to display.

        Raises:
            SystemExit: Exits the program with a usage error code.
        """
        # display error messages raised by argparse in red
        errormsg = (
            f"{Fore.RED}\nERROR: {error}.{Fore.RESET}"
            + f"\n\nRun {COMMAND} -h for usage\n\n"
        )
        sys.stderr.write(errormsg)  # write error to stderr
        sys.exit(os.EX_USAGE)  # exit execution -> usage error

    def error_noargs(self) -> None:
        """Display help and exit when no arguments are provided.

        Prints the help message and exits with a no input code.

        Raises:
            SystemExit: Exits the program with a no input error code.
        """
        self.print_help()  # if no input argument, print help
        sys.exit(os.EX_NOINPUT)  # exit with no input code


class CrispritzEnrichmentInputArgs:
    """Manage and validate enrichment analysis input arguments.

    This class validates input paths, numeric parameters, and verbosity settings,
    and exposes normalized values through convenient properties.
    """

    def __init__(self, args: Namespace, parser: CrispritzArgumentParser) -> None:
        """Initializes the CrispritzEnrichmentInputArgs with parsed arguments and
        parser.

        Stores the parsed arguments and parser, then checks argument consistency.
        """
        self._args = args
        self._parser = parser
        self._check_consistency()  # check input args consistency

    def _validate_vcf_folder(self) -> None:
        """Validate the VCF input folder and discover VCF files.

        Ensures the specified VCF directory exists and contains at least one
        compressed VCF file, storing the discovered files for later use.

        Returns:
            None
        """
        if not os.path.isdir(self._args.vcf):
            self._parser.error(f"Cannot find VCF folder {self._args.vcf}")
        self._vcfs = glob(os.path.join(self._args.vcf, "*.vcf.gz"))
        if self._args.vcf and not self._vcfs:
            self._parser.error(f"No VCF file found in {self._args.vcf}")

    def _validate_genome_folder(self) -> None:
        """Validate the genome input folder and discover FASTA files.

        Ensures the specified genome directory exists and contains at least one
        FASTA file, storing the discovered files for later use.

        Returns:
            None
        """
        if not os.path.exists(self._args.genome) or not os.path.isdir(
            self._args.genome
        ):
            self._parser.error(f"Cannot find input genome folder {self._args.genome}")
        self._fastas = glob(os.path.join(self._args.genome, "*.fa")) + glob(
            os.path.join(self._args.genome, "*.fasta")
        )
        if not self._fastas:
            self._parser.error(f"No FASTA file found in {self._args.genome}")

    def _validate_output_folder(self) -> None:
        """Validate the genome input folder and discover FASTA files.

        Ensures the specified genome directory exists and contains at least one
        FASTA file, storing the discovered files for later use.

        Returns:
            None
        """
        if not os.path.exists(self._args.outdir) or not os.path.isdir(
            self._args.outdir
        ):
            if not os.path.isdir(
                os.path.dirname(self._args.outdir)
            ):  # parent doesn't exist
                self._parser.error(f"Cannot find output folder {self._args.outdir}")
            os.makedirs(self._args.outdir)  # create output folder
        self._outdir = os.path.abspath(self._args.outdir)
        assert os.path.isdir(self._outdir)

    def _validate_threads(self) -> None:
        """Validates the thread count argument for allowed range.

        This function checks that the number of threads is non-negative and does
        not exceed the number of available CPU cores.

        Returns:
            None
        """
        max_threads = multiprocessing.cpu_count()
        if self._args.threads < 0 or self._args.threads > max_threads:
            self._parser.error(
                f"Forbidden number of threads provided ({self._args.threads}). "
                f"Max number of available cores: {max_threads}"
            )
        self._threads = max_threads if self._args.threads == 0 else self._args.threads

    def _validate_verbosity(self) -> None:
        """Validates the verbosity level argument for allowed values.

        This function checks that the verbosity level is one of the accepted values.

        Returns:
            None
        """
        if self._args.verbosity not in VERBOSITYLVL:
            self._parser.error(
                f"Forbidden verbosity level selected ({self._args.verbosity})"
            )

    def _check_consistency(self) -> None:
        """Validate all enrichment input arguments for consistency.

        Runs a sequence of validation steps to ensure all required input folders,
        thread settings, and verbosity options are correctly configured.

        Returns:
            None
        """
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
