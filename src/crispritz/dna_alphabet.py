""" """

from itertools import permutations


# ------------------------------------------------------------------------------
# Define DNA-related constant variables
# ------------------------------------------------------------------------------

# define dna alphabet
DNA = ["A", "C", "G", "T", "N"]

# define complete iupac alphabet
IUPAC = DNA + ["R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"]

# define reverse complement dictionary
RC = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "U": "A",
    "R": "Y",
    "Y": "R",
    "M": "K",
    "K": "M",
    "H": "D",
    "D": "H",
    "B": "V",
    "V": "B",
    "N": "N",
    "S": "S",
    "W": "W",
    "a": "t",
    "c": "g",
    "g": "c",
    "t": "a",
    "u": "a",
    "r": "y",
    "y": "r",
    "m": "k",
    "k": "m",
    "h": "d",
    "d": "h",
    "b": "v",
    "v": "b",
    "n": "n",
    "s": "s",
    "w": "w",
}

# define dictionary to encode nucleotides combinations as iupac characters
IUPACTABLE = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "M": "AC",
    "K": "GT",
    "S": "CG",
    "W": "AT",
    "H": "ACT",
    "B": "CGT",
    "V": "ACG",
    "D": "AGT",
    "N": "ACGT",
}

# define dictionary to encode nucleotide strings as iupac characters
IUPAC_ENCODER = {
    perm: k
    for k, v in IUPACTABLE.items()
    for perm in {"".join(p) for p in permutations(v)}
}

# ------------------------------------------------------------------------------
# Define DNA-related functions
# ------------------------------------------------------------------------------

def reverse_complement(sequence: str) -> str:
    return "".join(RC[nt] for nt in sequence[::-1])


def dna2rna(sequence: str) -> str:
    """Convert a DNA sequence to its RNA equivalent.

    Replaces all occurrences of 'T' with 'U' and 't' with 'u' in the input sequence.

    Args:
        sequence (str): The DNA sequence to convert.

    Returns:
        str: The RNA sequence.
    """
    return sequence.replace("T", "U").replace("t", "u")


