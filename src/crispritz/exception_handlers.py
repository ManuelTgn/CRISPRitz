""" """

import sys
import os


def sigint_handler() -> None:
    """Handle SIGINT (interrupt signal) to exit the program gracefully.

    Prints a message to standard error and exits the program with an OS error
    code when SIGINT is received.
    """
    # print message when SIGINT is caught to exit gracefully from the execution
    sys.stderr.write(f"\nCaught SIGINT. Exit CRISPRitz\n")
    sys.exit(os.EX_OSERR)  # mark as os error code
