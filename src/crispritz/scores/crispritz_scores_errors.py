""" """

from ..crispritz_errors import CrispritzError


class CfdScoreError(CrispritzError):
    def __init__(self, value: str):
        # initialize exception object when raised
        super().__init__(value)  # error message or error related info

    def __str__(self):
        return super().__str__()  # string representation for the exception