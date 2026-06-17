""" """

from ..crispritz_errors import CrispritzTstError, CrispritzSearchError
from ..exception_handlers import exception_handler
from ..guide import GuideList
from ..pam import PAM
from ..verbosity import VERBOSITY_LVL, print_verbosity

from typing import List

def search_offtargets_tst(indexes: List[str], pam_file: str, guides_file: str, mm: int, bdna: int, brna: int, outdir: str, threads: int, verbosity: int, debug: bool) -> None:
    pam = PAM(pam_file, debug)  # initialize pam
    guides = GuideList(guides_file, pam, debug)  # initialize guides

    # NOTE: each guide must be given as both fwd and reverse strand!!!




