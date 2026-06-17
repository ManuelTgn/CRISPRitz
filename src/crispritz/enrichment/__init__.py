from ..crispritz_inputargs import CrispritzEnrichmentInputArgs
from .enricher import enrich_genome

from typing import List

import os


def add_variants_cli(args: CrispritzEnrichmentInputArgs) -> None:
    enrich_genome(
        args.fastas,
        args.vcfs,
        args.keep,
        args.indels,
        False,  # option only for crisprme-integration
        args.outdir,
        args.threads,
        args.verbosity,
        args.debug,
    )


def add_variants(
    fasta_files: List[str],
    vcf_files: List[str],
    keep: bool = False,
    process_indels: bool = False,
    store_dictionary: bool = False,
    outdir: str = os.getcwd(),
    threads: int = 1,
    verbosity: int = 1,
    debug: bool = False,
) -> None:
    enrich_genome(
        fasta_files,
        vcf_files,
        keep,
        process_indels,
        store_dictionary,
        outdir,
        threads,
        verbosity,
        debug,
    )
