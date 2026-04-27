""" """

from ..crispritz_error import CrispritzTstError, CrispritzSearchError
from ..exception_handlers import exception_handler
from ..utils import print_verbosity, VERBOSITYLVL
from ..pam import PAM

def search_offtargets_tst(genome_index: str, pam_file: str, guides, mm: int, bdna: int, brna: int, outdir: str, threads: int, verbosity: int, debug: bool) -> None:
    pam = PAM(pam_file, debug)  # initialize pam_file




