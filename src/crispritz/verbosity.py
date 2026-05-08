""" """

import sys


# ------------------------------------------------------------------------------
# verbosity constant variables
# ------------------------------------------------------------------------------

# enabled verbosity levels:
# 0 = silent
# 1 = normal
# 2 = detail
# 3 = debug
VERBOSITY_LVL = [0, 1, 2, 3]  


# ------------------------------------------------------------------------------
# verbosity functions
# ------------------------------------------------------------------------------

def print_verbosity(message: str, verbosity: int, verbosity_threshold: int) -> None:
    """Print a message if the verbosity level meets the threshold.

    Writes the message to standard output if the current verbosity is greater
    than or equal to the specified threshold.

    Parameters
    ----------
    message: 
        The message to print.
    verbosity: 
        The current verbosity level.
    verbosity_threshold: 
        The minimum verbosity level required to print the message.
    """
    if verbosity >= verbosity_threshold:
        sys.stdout.write(f"{message}\n")
    return








