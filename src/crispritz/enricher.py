""" """

from .crispritz_error import CrispritzEnrichmentError, EnrichmentPairError
from .exception_handlers import exception_handler
from .genome_io import GenomeReader, GenomeWriter
from .crispritz_argparse import CrispritzEnrichmentInputArgs
from .utils import print_verbosity, create_folder, VERBOSITYLVL

from typing import List, Dict, Set, Optional, Tuple
from pysam import FastaFile, VariantFile
from io import TextIOWrapper
from time import time

import gzip
import os


# output folders
VARIANTGENOMEDIR = "variants_genome"  # root folder
SNPDIR = os.path.join(VARIANTGENOMEDIR, "SNPs_genome")  # snps genome
INDELSDIR = os.path.join(VARIANTGENOMEDIR, "INDELs_genome")  # indels genome


class EnrichPair:

    def __init__(self, debug: bool) -> None:
        self._debug = debug  # store debug flag
        self._fasta: Optional[str] = None
        self._vcf: Optional[str] = None

    def __repr__(self) -> str:
        """Return a detailed string representation for debugging."""
        return f"EnrichPair(fasta={self._fasta!r}, vcf={self._vcf!r})"

    def __str__(self) -> str:
        """Return a human-readable string representation."""
        fasta_str = self._fasta if self._fasta is not None else "not set"
        vcf_str = self._vcf if self._vcf is not None else "not set"
        return f"EnrichPair: FASTA={fasta_str}, VCF={vcf_str}"

    @property
    def fasta(self) -> Optional[str]:
        return self._fasta

    @fasta.setter
    def fasta(self, value: str) -> None:
        if not isinstance(value, str) or not value:
            exception_handler(
                EnrichmentPairError,
                "FASTA file must be a non-empty str, got "
                f"{type(value).__name__} instead",
                os.EX_DATAERR,
                self._debug,
            )
        self._fasta = value

    @property
    def vcf(self) -> Optional[str]:
        return self._vcf

    @vcf.setter
    def vcf(self, value: str) -> None:
        if not isinstance(value, str):
            exception_handler(
                EnrichmentPairError,
                "VCF file must be a non-empty str, got "
                f"{type(value).__name__} instead",
                os.EX_DATAERR,
                self._debug,
            )
        self._vcf = value


def _retrieve_contig_name(fasta: FastaFile, debug: bool) -> str:
    if len(fasta.references) != 1:  # assumes each fasta is chromosome separated
        contigs = ", ".join(fasta.references)
        exception_handler(
            CrispritzEnrichmentError,
            f"FASTA {fasta.filename} contains multiple contigs: {contigs}. Each "
            "FASTA is expected to contain exactly one contig",
            os.EX_DATAERR,
            debug,
        )
    contig = fasta.references[0]  # assumes single contig in fasta
    return contig if contig.startswith("chr") else f"chr{contig}"


def retrieve_contig_names(
    fasta_files: List[str], verbosity: int, debug: bool
) -> Set[str]:
    # retrieve contig names for each fasta file in genome folder
    print_verbosity(
        "Retrieving contig names from FASTA in genome folder",
        verbosity,
        VERBOSITYLVL[3],
    )
    return {_retrieve_contig_name(FastaFile(f), debug) for f in fasta_files}


def initialize_fasta(
    fasta_vcf_map: Dict[str, EnrichPair], fasta_files: List[str], debug: bool
) -> Dict[str, EnrichPair]:
    for f in fasta_files:
        contig = FastaFile(f).references[0]  # retrieve contig name
        if not contig.startswith("chr"):
            contig = f"chr{contig}"  # avoid mismatch (see 1000G)
        if (
            fasta_vcf_map[contig].fasta is not None
        ):  # multiple fasta pointing to same contig
            exception_handler(
                CrispritzEnrichmentError,
                f"Multiple FASTA file pointing to contig {contig}: "
                f"{f} - {fasta_vcf_map[contig].fasta}",
                os.EX_DATAERR,
                debug,
            )
        fasta_vcf_map[contig].fasta = f  # assign fasta slot
    return fasta_vcf_map


def _retrieve_contig_vcf(vcf: VariantFile, debug: bool) -> str:
    contigs = list(map(str, vcf.header.contigs))  # get contig names in vcf
    if len(contigs) != 1:
        contigs = ", ".join(contigs)
        exception_handler(
            CrispritzEnrichmentError,
            f"VCF file {vcf.filename} contains multiple contigs: {contigs}. Each "
            "VCF is expected to contain exactly one contig",
            os.EX_DATAERR,
            debug,
        )
    contig = contigs[0]  # assumes single contig in vcf
    return contig if contig.startswith("chr") else f"chr{contig}"


def initialize_vcf(
    fasta_vcf_map: Dict[str, EnrichPair], vcf_files: List[str], debug: bool
) -> Dict[str, EnrichPair]:
    for f in vcf_files:
        contig = _retrieve_contig_vcf(
            VariantFile(f, mode="r"), debug
        )  # retrieve vcf contig
        if (
            fasta_vcf_map[contig].vcf is not None
        ):  # multiple vcf pointing to same contig
            exception_handler(
                CrispritzEnrichmentError,
                f"Multiple VCF file pointing to conting {contig}: "
                f"{f} - {fasta_vcf_map[contig].vcf}",
                os.EX_DATAERR,
                debug,
            )
        fasta_vcf_map[contig].vcf = f  # assign vcf slot
    return fasta_vcf_map


def construct_fasta_vcf_map(
    fasta_files: List[str], vcf_files: List[str], verbosity: int, debug: bool
) -> Dict[str, EnrichPair]:
    # retrieve contig names in fasta from genome folder and initialize map
    fasta_vcf_map = {
        contig: EnrichPair(debug)
        for contig in retrieve_contig_names(fasta_files, verbosity, debug)
    }
    # initialize fasta elements in the map
    fasta_vcf_map = initialize_fasta(fasta_vcf_map, fasta_files, debug)
    # initialize vcf elements in the map
    fasta_vcf_map = initialize_vcf(fasta_vcf_map, vcf_files, debug)
    return fasta_vcf_map


def split_contigs(fasta_vcf_map: Dict[str, EnrichPair]) -> Tuple[List[str], List[str]]:
    # divide contigs with variants from those without for different processing
    contigs_vcf = [contig for contig, p in fasta_vcf_map.items() if p.vcf is not None]
    contigs_wo_vcf = [contig for contig, p in fasta_vcf_map.items() if p.vcf is None]
    return contigs_vcf, contigs_wo_vcf


def prepare_output_dir(outdir: str) -> Tuple[str, str]:
    return create_folder(os.path.join(outdir, SNPDIR)), create_folder(
        os.path.join(outdir, INDELSDIR)
    )


def enrich_no_variants(
    fasta_vcf_map: Dict[str, EnrichPair], contigs: List[str], outdir: str, debug: bool
) -> None:
    for contig in contigs:  # just copy fasta without variants for enrichment
        assert fasta_vcf_map[contig].vcf is None
        reader = GenomeReader(fasta_vcf_map[contig].fasta, debug)  # type: ignore
        reader.read()  # read contig sequence
        # define output fasta filename
        prefix = os.path.splitext(os.path.basename(fasta_vcf_map[contig].fasta))[0]  # type: ignore
        fasta_enr = os.path.join(outdir, f"{prefix}.enriched.fa")
        writer = GenomeWriter(fasta_enr, debug)  # write contig sequence
        writer.write(reader.header, reader.sequence)  # type: ignore

def _extract_samples(vcf_header: List[str], vcf_fname: str, debug: bool) -> List[str]:
    try:
        return vcf_header[9:]
    except Exception as e:
        exception_handler(CrispritzEnrichmentError, f"Failed retrieving samples from VCF header: {vcf_fname}", os.EX_IOERR, debug, e)

def retrieve_samples(vcfin: TextIOWrapper, vcf_fname: str, debug: bool) -> List[str]:
    for line in vcfin:  # parse VCF header
        if "#CHROM" in line:  # end of header reached
            header = line.strip().split()  # store header to retrieve samples and af data
            # retrieve samples and allele frequency information from vcfs
            return _extract_samples(header, vcf_fname, debug)
    # header not found?
    exception_handler(CrispritzEnrichmentError, f"VCF header parsing failed: {vcf_fname}", os.EX_IOERR, debug)

def _extract_af(info: str, vcf_fname: str, debug: bool) -> int:
    for i, e in enumerate(info.split(";")):  # look in INFO field
        if e[:2] == "AF":
            return i
    exception_handler(CrispritzEnrichmentError, f"Failed retrieving AF index from VCF: {vcf_fname}", os.EX_IOERR, debug)

def _skip_variant(variant_filter: str) -> bool:
    return variant_filter != "PASS"

def insert_snps(vcfin: TextIOWrapper, vcf_fname: str, debug: bool):
    for line in vcfin:  # iterate over variants
        variant = line.strip().split("\t")  # split variant in its fields
        if _skip_variant(variant[6]):  # filter != PASS
            continue
        afidx = _extract_af(variant[7], vcf_fname, debug)

    

def enrich_variants(fasta_vcf_map: Dict[str, EnrichPair], contigs: List[str], outdir: str, debug: bool) -> None:
    for contig in contigs:
        reader = GenomeReader(fasta_vcf_map[contig].fasta, debug)  # type: ignore
        reader.read()  # read contig sequence
        with gzip.open(fasta_vcf_map[contig].vcf, mode="rt") as fin:  # type: ignore
            samples = retrieve_samples(fin, fasta_vcf_map[contig].vcf, debug) # type: ignore
            #TODO: indels






def enrich_genome(fasta_vcf_map: Dict[str, EnrichPair], outdir: str, debug: bool):
    # retrieve contig to enrich with variants and those without variants associated
    contigs_vcf, contigs_wo_vcf = split_contigs(fasta_vcf_map)
    snpsdir, indelsdir = prepare_output_dir(outdir)  # prepare enrichment output folder
    # copy content of original fasta for contig without variants
    enrich_no_variants(fasta_vcf_map, contigs_wo_vcf, snpsdir, debug)


def add_variants(args: CrispritzEnrichmentInputArgs) -> None:
    # construct a fasta-vcf files map
    fasta_vcf_map = construct_fasta_vcf_map(
        args.fastas, args.vcfs, args.verbosity, args.debug
    )
    start = time()  # genome enrichment start point
    print_verbosity(
        "Enriching genome with input variants", args.verbosity, VERBOSITYLVL[1]
    )
    enrich_genome(fasta_vcf_map, args.outdir, args.debug)  # genome enrichment
