""" """

from typing import Dict, FrozenSet


class IUPACCodes:
    
    # Mapping from IUPAC code to all possible nucleotide combinations
    CODE_TO_BASES: Dict[str, FrozenSet[str]] = {
        'R': frozenset(['AG', 'GA']),
        'Y': frozenset(['CT', 'TC']),
        'S': frozenset(['GC', 'CG']),
        'W': frozenset(['AT', 'TA']),
        'K': frozenset(['GT', 'TG']),
        'M': frozenset(['AC', 'CA']),
        'B': frozenset(['CGT', 'GCT', 'TGC', 'GTC', 'CTG', 'TCG']),
        'D': frozenset(['AGT', 'GAT', 'TAG', 'ATG', 'GTA', 'TGA']),
        'H': frozenset(['ACT', 'CAT', 'TCA', 'ATC', 'CTA', 'TAC']),
        'V': frozenset(['ACG', 'CAG', 'GAC', 'AGC', 'CGA', 'GCA']),
        'N': frozenset([
            'ACGT', 'CAGT', 'GACT', 'AGCT', 'CGAT', 'GCAT',
            'GCTA', 'CGTA', 'TGCA', 'GTCA', 'CTGA', 'TCGA',
            'TAGC', 'ATGC', 'GTAC', 'TGAC', 'AGTC', 'GATC',
            'CATG', 'ACTG', 'TCAG', 'CTAG', 'ATCG', 'TACG'
        ])
    }
    
    # Simple decomposition mapping
    DECOMPOSITION: Dict[str, str] = {
        "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
        "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
        "H": "ACT", "V": "ACG", "N": "ATGC",
        "A": "A", "T": "T", "C": "C", "G": "G"
    }
    
    @classmethod
    def get_iupac_code(cls, nucleotides: str) -> str:
        nucleotides_sorted = ''.join(sorted(set(nucleotides.upper())))        
        for code, combinations in cls.CODE_TO_BASES.items():
            if nucleotides_sorted in combinations:
                return code
        return 'N'  # Default to 'N' if no match found

