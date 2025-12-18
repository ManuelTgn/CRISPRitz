""" """

from .utils import find_fasta_index

from pysam import faidx, FastaFile
from typing import Tuple, List
from pathlib import Path

import pysam

class GenomeReader:

    # class attributes
    _sequence: List[str]
    _contig: str
    
    # class functions
    def __init__(self, fasta_path: str):
        self._fasta_path = fasta_path
    
    def read(self) -> None:
        if not find_fasta_index(self._fasta_path):  # Search for FASTA index
            faidx(self._fasta_path)  # index not found, compute it
        fa = FastaFile(self._fasta_path)  # load fasta with pysam
        self._contig = fa.references[0]  # assume single contig fasta
        self._sequence = list(fa.fetch(self._contig))

    @property
    def contig(self) -> str:
        return self._contig
    
    @property
    def sequence(self) -> List[str]:
        return self._sequence
               

class GenomeWriter:
    
    @staticmethod
    def write_enriched_genome(
        output_path: Path,
        header: str,
        sequence_list: List[str]
    ) -> None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        sequence = ''.join(sequence_list)
        with open(output_path, 'w') as f:
            if not header.startswith('>'):
                header = '>' + header
            if not header.endswith('\n'):
                header += '\n'
            
            f.write(header)
            f.write(sequence)
            f.write('\n')
