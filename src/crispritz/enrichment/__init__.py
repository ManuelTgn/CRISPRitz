from ..crispritz_argparse import CrispritzEnrichmentInputArgs
from .enricher import enrich_genome

from typing import List

def add_variants_cli(args: CrispritzEnrichmentInputArgs) -> None:
    enrich_genome(args.fastas, args.vcfs, args.verbosity, args.debug)

def add_variants(fasta_files: List[str], vcf_files: List[str], verbosity: int = 1, debug: bool = False) -> None:
    enrich_genome(fasta_files, vcf_files, verbosity, debug)