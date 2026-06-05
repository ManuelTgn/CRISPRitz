""" """

from typing import List


from ..offtarget import OffTarget

from .cfd import load_mismatch_pam_scores, compute_cfd


def cfd_score(offtargets: List[OffTarget], debug: bool) -> List[float]:
    mmscores, pamscores = load_mismatch_pam_scores(debug)  # load scoring models
    return [compute_cfd(ot.grna, ot.spacer, ot.spacer[-2:], mmscores, pamscores) for ot in offtargets]


