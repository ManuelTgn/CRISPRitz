""" """

from ..crispritz_error import CrispritzEnrichmentError, EnrichmentPairError
from ..exception_handlers import exception_handler
from .genome_io import GenomeReader, GenomeWriter, INDELOFFSET
from ..dna_alphabet import IUPAC_ENCODER, IUPACTABLE
from ..crispritz_argparse import CrispritzEnrichmentInputArgs
from ..utils import print_verbosity, create_folder, find_tabix_index, VERBOSITYLVL
from .enrichment_pair import EnrichPair
from .variants import Snp, Snps, Indel, Indels, IndelsSet, IndelPair, IndelInfo

from typing import List, Dict, Set, Tuple, Union
from pysam import FastaFile, VariantFile, tabix_index
from pysam.utils import SamtoolsError
from io import TextIOWrapper
from dataclasses import dataclass
from time import time

import json
import gzip
import os


# output folders
VARIANTGENOMEDIR = "variants_genome"  # root folder
SNPDIR = os.path.join(VARIANTGENOMEDIR, "SNPs_genome")  # snps genome
INDELSDIR = os.path.join(VARIANTGENOMEDIR, "INDELs_genome")  # indels genome




def _retrieve_contig_name(fasta: FastaFile, debug: bool) -> str:
    """Retrieve the single contig name from a FASTA file and normalize it. The
    contig name is validated for uniqueness and adjusted to start with 'chr'
    when needed.

    This function checks that the given FASTA file contains exactly one contig,
    raises an error otherwise, and returns the standardized contig identifier.

    Args:
        fasta: An open pysam.FastaFile object from which to read contig names.
        debug: Flag indicating whether to use debug-aware error handling.

    Returns:
        The normalized contig name extracted from the FASTA file.
    """
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
    """Retrieve standardized contig names from a list of FASTA files. Each FASTA
    file is expected to contain exactly one contig.

    This function reads the headers of the provided FASTA files, extracts their
    contig names, normalizes them to start with 'chr', and returns the set of
    unique contig identifiers found.

    Args:
        fasta_files: List of paths to FASTA files representing genome contigs.
        verbosity: Verbosity level controlling printed progress information.
        debug: Flag indicating whether to use debug-aware error handling.

    Returns:
        A set of normalized contig names extracted from the FASTA files.
    """
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
    """Populate the FASTA entries of a contig-to-file mapping. This associates each
    contig key with exactly one corresponding FASTA file.

    This function inspects each input FASTA, derives its contig name in a
    normalized 'chr' form, checks for duplicate assignments and records the file
    path in the mapping.

    Args:
        fasta_vcf_map: Dictionary mapping contig names to `EnrichPair` objects 
            to be updated with FASTA paths.
        fasta_files: List of FASTA file paths to register in the mapping.
        debug: Flag indicating whether to use debug-aware error handling.

    Returns:
        The updated contig-to-file mapping including the assigned FASTA paths.
    """
    for f in fasta_files:
        contig = FastaFile(f).references[0]  # retrieve contig name
        if not contig.startswith("chr"):
            contig = f"chr{contig}"  # avoid mismatch (see 1000G)
        # multiple fasta pointing to same contig
        if (fasta_vcf_map[contig].fasta is not None):
            exception_handler(
                CrispritzEnrichmentError,
                f"Multiple FASTA file pointing to contig {contig}: "
                f"{f} - {fasta_vcf_map[contig].fasta}",
                os.EX_DATAERR,
                debug,
            )
        fasta_vcf_map[contig].fasta = f  # assign fasta slot
    return fasta_vcf_map


def _tabix_index(vcf_fname: str, verbosity: int, debug: bool) -> None:
    """Ensure that a VCF file has an associated tabix index. This prepares the VCF
        for random access during downstream enrichment.

        The function checks for an existing index, creates one if missing, and uses
        debug-aware error handling to report failures.

        Args:
            vcf_fname: Path to the VCF file to be indexed.
            verbosity: Verbosity level controlling printed progress information.
            debug: Flag indicating whether to use debug-aware error handling.

        Returns:
            None
    """
    if find_tabix_index(vcf_fname):  # index found, do nothing
        return
    try:  # tabix index not found, compute index
        print_verbosity(
            f"Index not found, indexing (VCF: {vcf_fname})", verbosity, VERBOSITYLVL[3]
        )
        tabix_index(vcf_fname, preset="vcf")
    except (SamtoolsError, Exception) as e:
        exception_handler(
            CrispritzEnrichmentError,
            f"Failed indexing VCF: {vcf_fname}",
            os.EX_DATAERR,
            debug,
            e,
        )


def _retrieve_contig_vcf(vcf: VariantFile, debug: bool) -> str:
    """Retrieve the single contig name from a VCF file and normalize it. The
    contig name is validated for uniqueness and adjusted to start with 'chr'
    when needed.

    This function checks that the given VCF file declares exactly one contig in
    its header, raises an error otherwise, and returns the standardized contig
    identifier.

    Args:
        vcf: An open pysam.VariantFile object from which to read contig names.
        debug: Flag indicating whether to use debug-aware error handling.

    Returns:
        The normalized contig name extracted from the VCF header.
    """
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
    fasta_vcf_map: Dict[str, EnrichPair],
    vcf_files: List[str],
    verbosity: int,
    debug: bool,
) -> Dict[str, EnrichPair]:
    """Populate the VCF entries of a contig-to-file mapping. This associates each
    contig key with exactly one corresponding VCF file containing its variants.

    This function indexes each VCF if needed, derives its contig name in a
    normalized 'chr' form, checks for duplicate assignments and records the file
    path in the mapping.

    Args:
        fasta_vcf_map: Dictionary mapping contig names to `EnrichPair` objects 
            to be updated with VCF paths.
        vcf_files: List of VCF file paths to register in the mapping.
        verbosity: Verbosity level controlling printed progress information.
        debug: Flag indicating whether to use debug-aware error handling.

    Returns:
        The updated contig-to-file mapping including the assigned VCF paths.
    """
    for f in vcf_files:
        # retrieve vcf contig
        _tabix_index(f, verbosity, debug)
        contig = _retrieve_contig_vcf(VariantFile(f, mode="r"), debug)
        # multiple vcf pointing to same contig
        if fasta_vcf_map[contig].vcf is not None:
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
    """Build a mapping between contig names and their FASTA/VCF file pairs. This
    prepares per-contig inputs needed for downstream genome enrichment.

    The function discovers contigs from the FASTA files, initializes an
    `EnrichPair` for each, then populates the mapping with corresponding FASTA
    and VCF paths.

    Args:
        fasta_files: List of FASTA file paths representing contig sequences.
        vcf_files: List of VCF file paths providing variant calls per contig.
        verbosity: Verbosity level controlling printed progress information.
        debug: Flag indicating whether to use debug-aware error handling.

    Returns:
        A dictionary mapping normalized contig names to their associated
        `EnrichPair` instances.
    """
    # retrieve contig names in fasta from genome folder and initialize map
    fasta_vcf_map = {
        contig: EnrichPair(debug)
        for contig in retrieve_contig_names(fasta_files, verbosity, debug)
    }
    # initialize fasta elements in the map
    fasta_vcf_map = initialize_fasta(fasta_vcf_map, fasta_files, debug)
    # initialize vcf elements in the map
    fasta_vcf_map = initialize_vcf(fasta_vcf_map, vcf_files, verbosity, debug)
    return fasta_vcf_map


def split_contigs(
    fasta_vcf_map: Dict[str, EnrichPair], verbosity: int
) -> Tuple[List[str], List[str]]:
    """Split contigs into those with and without associated VCF files. This helps
    route contigs through variant-aware or copy-only enrichment workflows.

    This function inspects the contig-to-file mapping, groups contigs based on
    whether a VCF path is present, logs the group sizes, and returns the two
    resulting contig lists.

    Args:
        fasta_vcf_map: Mapping from contig names to `EnrichPair` objects
            containing FASTA and VCF file associations.
        verbosity: Verbosity level controlling printed progress information.

    Returns:
        A tuple containing a list of contigs with VCFs and a list of contigs
        without VCFs.
    """
    # divide contigs with variants from those without for different processing
    print_verbosity(
        "Retrieving contigs with associated VCFs", verbosity, VERBOSITYLVL[3]
    )
    contigs_vcf = [contig for contig, p in fasta_vcf_map.items() if p.vcf is not None]
    contigs_wo_vcf = [contig for contig, p in fasta_vcf_map.items() if p.vcf is None]
    print_verbosity(
        f"Contigs with VCFs: {len(contigs_vcf)}, contigs without VCFs: {len(contigs_wo_vcf)}",
        verbosity,
        VERBOSITYLVL[3],
    )
    return contigs_vcf, contigs_wo_vcf


def prepare_output_dir(outdir: str) -> Tuple[str, str]:
    """Prepare the output directory structure for enrichment results. This ensures
    separate subfolders exist for SNP and INDEL enriched genomes.

    The function creates (or reuses) the SNP and INDEL output folders inside the
    given root directory and returns their paths.

    Args:
        outdir: Root output directory where variant genome subfolders are created.

    Returns:
        A tuple containing the SNP output directory path and the INDEL output
        directory path, in that order.
    """
    # create snps and indels out directory
    return create_folder(os.path.join(outdir, SNPDIR)), create_folder(
        os.path.join(outdir, INDELSDIR)
    )


def enrich_no_variants(
    fasta_vcf_map: Dict[str, EnrichPair],
    contigs: List[str],
    outdir: str,
    verbosity: int,
    debug: bool,
) -> None:
    """Enrich contigs without variant data by copying their reference sequences. 
    This prepares consistent enriched FASTA outputs for contigs lacking VCFs.

    The function iterates over the provided contigs, reads each reference FASTA,
    and writes an 'enriched' copy to the output directory while reporting
    progress.

    Args:
        fasta_vcf_map: Mapping from contig names to `EnrichPair` objects containing 
            FASTA paths and optional VCF paths.
        contigs: List of contig names to process that have no associated VCF files.
        outdir: Directory where enriched FASTA files will be written.
        verbosity: Verbosity level controlling printed progress information.
        debug: Flag indicating whether to use debug-aware error handling.

    Returns:
        None
    """
    for contig in contigs:  # just copy fasta without variants for enrichment
        print_verbosity(f"Enriching contig {contig}", verbosity, VERBOSITYLVL[3])
        start = time()  # track enrichment running time
        assert fasta_vcf_map[contig].vcf is None
        reader = GenomeReader(fasta_vcf_map[contig].fasta, debug)  # type: ignore
        reader.read()  # read contig sequence
        # define output fasta filename
        prefix = os.path.splitext(os.path.basename(fasta_vcf_map[contig].fasta))[0]  # type: ignore
        fasta_enr = os.path.join(outdir, f"{prefix}.enriched.fa")
        writer = GenomeWriter(fasta_enr, debug)  # write contig sequence
        writer.write(reader.header, reader.sequence)  # type: ignore
        print_verbosity(
            f"Enrichment on contig  {contig} completed in {time() - start:.2f}s",
            verbosity,
            VERBOSITYLVL[3],
        )


def _extract_samples(vcf_header: List[str], vcf_fname: str, debug: bool) -> List[str]:
    """Extract sample identifiers from a parsed VCF header line. This assumes the
    header is already split into fields and contains standard VCF columns.

    The function returns all columns after the fixed FORMAT column, or raises a
    debug-aware error if extraction fails.

    Args:
        vcf_header: List of header fields from the '#CHROM' VCF header line.
        vcf_fname: Path to the VCF file being parsed, used for error reporting.
        debug: Flag indicating whether to use debug-aware error handling.

    Returns:
        A list of sample names defined in the VCF header.
    """
    try:
        return vcf_header[9:]
    except Exception as e:
        exception_handler(
            CrispritzEnrichmentError,
            f"Failed retrieving samples from VCF header: {vcf_fname}",
            os.EX_IOERR,
            debug,
            e,
        )


def retrieve_samples(vcfin: TextIOWrapper, vcf_fname: str, debug: bool) -> List[str]:
    """Read a VCF stream and extract the list of sample names. This scans header
    lines until the '#CHROM' line is found and then parses its fields.

    The function delegates to `_extract_samples` to pull out sample identifiers
    and raises a debug-aware error if the header cannot be located.

    Args:
        vcfin: Open text stream for the VCF file, positioned at the start.
        vcf_fname: Path to the VCF file being parsed, used for error reporting.
        debug: Flag indicating whether to use debug-aware error handling.

    Returns:
        A list of sample names defined in the VCF header.
    """
    for line in vcfin:  # parse VCF header
        if "#CHROM" in line:  # end of header reached
            header = (
                line.strip().split()
            )  # store header to retrieve samples and af data
            # retrieve samples and allele frequency information from vcfs
            return _extract_samples(header, vcf_fname, debug)
    # header not found?
    exception_handler(
        CrispritzEnrichmentError,
        f"VCF header parsing failed: {vcf_fname}",
        os.EX_IOERR,
        debug,
    )



def _skip_variant(variant_filter: str) -> bool:
    """Decide whether a VCF variant should be skipped based on its FILTER field. 
    This treats only 'PASS' variants as eligible for enrichment.

    The function returns True for filtered-out variants and False for variants
    that passed quality control.

    Args:
        variant_filter: The FILTER column value from a VCF record.

    Returns:
        True if the variant should be skipped, False if it should be processed.
    """
    return variant_filter != "PASS"


def _extract_af_idx(info: str, debug: bool) -> int:
    """Locate the index of the allele-frequency (AF) field in a VCF INFO string. 
    This identifies which semicolon-separated entry encodes AF values.

    The function scans the INFO components, returns the position of the first
    entry starting with 'AF', and raises an error if no such entry is found.

    Args:
        info: The INFO column string from a VCF record.
        debug: Flag indicating whether to use debug-aware error handling.

    Returns:
        The zero-based index of the AF entry within the semicolon-separated INFO 
            fields.
    """
    for i, e in enumerate(info.split(";")):  # look in INFO field
        if e[:2] == "AF":
            return i
    exception_handler(
        CrispritzEnrichmentError, "Failed retrieving AF index", os.EX_IOERR, debug
    )

def _split_snps_indels(pos: int, ref: str, alts: str) -> Tuple[Snps, Indels]:
    """Separate alternate alleles into SNPs and indels relative to a reference base. 
    This prepares per-type containers that drive downstream enrichment logic.

    The function walks over all alternate alleles, classifies single-base
    substitutions as SNPs and all other length-changing alleles as indels, and
    records their shared position and per-allele genotype index.

    Args:
        pos: Zero-based genomic position of the variant.
        ref: Reference allele sequence from the VCF record.
        alts: Comma-separated string of alternate alleles from the VCF record.

    Returns:
        A tuple containing a `Snps` collection and an `Indels` collection built
        from the provided alleles.
    """
    # retrieve reference, snps and indels for current variant
    snps, indels = Snps(), Indels()  # snps and indels containers
    for i, alt in enumerate(alts.strip().split(",")):
        if len(alt) == len(ref) == 1:  # snp found
            snps.add(Snp(pos, ref, alt, i))
        else:  # indel found
            indels.add(Indel(pos, ref, alt, i))
    return snps, indels


def _compute_vid(chrom: str, pos: Union[int, str], ref: str, alt: str) -> str:
    """Construct a stable identifier string for a variant. This encodes the
    chromosome, position and allele change in a compact, comparable form.

    The function normalizes chromosome names to start with 'chr' and joins all
    components into a single 'chr-pos-ref/alt' label.

    Args:
        chrom: Chromosome name from the variant record.
        pos: Genomic position of the variant, either as int or string.
        ref: Reference allele sequence.
        alt: Alternate allele sequence.

    Returns:
        A variant identifier string of the form 'chrX-pos-ref/alt'.
    """
    chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
    return f"{chrom}-{pos}-{ref}/{alt}"


def _retrieve_carriers(genotypes: List[str], samples: List[str], gtidx: str, indels: bool = False) -> str:
    """Identify which samples carry a specific allele genotype. This can optionally 
    suppress genotype details for indel reporting.

    The function scans all genotype entries, selects those whose leading allele
    index matches the requested value, and returns either 'sample' or
    'sample:genotype' labels depending on the indels flag.

    Args:
        genotypes: List of genotype strings (e.g. '0/1:...') for each sample.
        samples: List of sample names aligned with the genotype list.
        gtidx: Genotype index (as a string) that identifies the allele of interest.
        indels: Flag indicating whether to omit genotype strings from the output.

    Returns:
        A comma-separated string of carrier labels, or an empty string if no
        carriers are found.
    """
    if indels:  # do not report genotype associated to samples
        carriers = [
            f"{samples[i]}"
            for i, gt in enumerate(genotypes)
            if gtidx in (g := gt.split(":")[0])
        ]
    else:  # report genotype associated to samples
        carriers = [
            f"{samples[i]}:{g}"
            for i, gt in enumerate(genotypes)
            if gtidx in (g := gt.split(":")[0])
        ]
    return ",".join(sorted(carriers))

def _retrieve_af(info: str, idx: int, gtidx: int) -> str:
    """Retrieve the allele-frequency value for a specific alternate allele. This
    operates on a parsed INFO field and an index previously identified for AF.

    The function selects the AF entry at the given position, strips its 'AF='
    prefix, and returns the value corresponding to the requested genotype index.

    Args:
        info: The INFO column string from a VCF record.
        idx: Zero-based index of the AF entry within the semicolon-separated INFO 
            fields.
        gtidx: One-based genotype index pointing to the target alternate allele.

    Returns:
        The allele-frequency string for the specified alternate allele.
    """
    return info.split(";")[idx][3:].split(",")[gtidx - 1]


def _create_snp_dict_entry(carriers: str, alleles: str, vid: str, af: str) -> str:
    """Assemble a compact dictionary entry string for a SNP. This encodes carrier 
    samples, alleles, a variant identifier and allele frequency in a single record.

    The function conditionally prefixes the entry with carrier information and
    always includes alleles, variant ID and AF separated by semicolons.

    Args:
        carriers: Comma-separated 'sample:genotype' labels for carrier samples, 
            or an empty string.
        alleles: Reference and alternate alleles encoded as 'REF,ALT'.
        vid: Stable variant identifier string (e.g. from `_compute_vid`).
        af: Allele-frequency value for the SNP.

    Returns:
        A semicolon-delimited string representing the SNP dictionary entry.
    """
    if carriers:
        return f"{carriers};{alleles};{vid};{af}"
    return f";{alleles};{vid};{af}"


def insert_snp_in_dict(
    chrom_snps_dict: Dict,
    contig: str,
    info: str,
    genotypes: List[str],
    samples: List[str],
    afidx: int,
    snps: Snps,
) -> Dict[str, str]:
    """Insert SNP information for a single genomic position into the chromosome 
    SNP dictionary. This consolidates multiallelic SNP data into a compact, 
    per-position entry.

    The function builds per-allele records containing carriers, alleles, a
    variant identifier and allele frequency, joins them when multiple alleles
    share the same position, and updates the dictionary in place.

    Args:
        chrom_snps_dict: Dictionary storing SNP annotations keyed by 'contig,pos'.
        contig: Normalized contig name for the SNPs being recorded.
        info: INFO column string from the VCF record supplying allele-frequency 
            values.
        genotypes: List of genotype strings for all samples at this record.
        samples: List of sample names aligned with the genotype list.
        afidx: Zero-based index of the AF entry within the semicolon-separated 
            INFO fields.
        snps: Collection of `Snp` objects representing all SNP alleles at this 
            position.

    Returns:
        The updated chromosome SNP dictionary containing the new position entry.
    """
    snpkey = f"{contig},{snps.pos()}"  # retrieve snp key
    # compute dictionary entry for each snp (multiallelic sites)
    entries = []
    for snp in snps.items:
        # retrieve snp carriers
        carriers = _retrieve_carriers(genotypes, samples, str(snp.gtidx))
        af = _retrieve_af(info, afidx, snp.gtidx)  # snp af
        alleles = f"{snp.ref},{snp.alt}"  # snp alleles
        vid = _compute_vid(contig, snp.pos, snp.ref, snp.alt)  # compute id
        entries.append(_create_snp_dict_entry(carriers, alleles, vid, af))
    # join multiallelic snps on same dictionary entry
    chrom_snps_dict[snpkey] = "$".join(entries) if len(entries) > 1 else entries[0]
    return chrom_snps_dict


def _process_snp(
    variant: List[str],
    snps: Snps,
    contig: str,
    reader: GenomeReader,
    chrom_snps_dict: Dict[str, str],
    samples: List[str],
    afidx: int,
    store_dictionary: bool,
    debug: bool,
) -> Dict[str, str]:
    """Apply SNP alleles from a single VCF record to the in-memory contig sequence. 
    This both updates the enriched sequence and optionally records dictionary 
    metadata for the SNPs.

    The function validates reference bases against the FASTA, encodes the
    combined alleles using an IUPAC symbol, inserts that symbol into the
    sequence, and delegates SNP dictionary insertion when requested.

    Args:
        variant: Full list of VCF fields for the current record.
        snps: Collection of `Snp` objects representing all SNP alleles at this 
            position.
        contig: Normalized contig name for the current record.
        reader: GenomeReader instance holding the contig sequence to be enriched.
        chrom_snps_dict: Dictionary storing SNP annotations keyed by 'contig,pos'.
        samples: List of sample names aligned with genotype fields in the VCF record.
        afidx: Zero-based index of the AF entry within the semicolon-separated 
            INFO fields.
        store_dictionary: Flag indicating whether SNP metadata should be stored 
            in the dictionary.
        debug: Flag indicating whether to use debug-aware error handling.

    Returns:
        The updated chromosome SNP dictionary, potentially with a new entry for 
            this position.
    """
    # retrieve ref allele from contig sequence
    pos = snps.pos()  # snp position
    ref_nt = reader.sequence[pos]
    ref = snps.ref()  # snp reference allele
    if snps.ref() not in IUPACTABLE[ref_nt]:  # mismatch between VCF and contig FASTA data
        vid = _compute_vid(contig, pos, ref, ",".join(snps.alts()))
        exception_handler(
            CrispritzEnrichmentError,
            f"Mismatching REF alleles in VCF and FASTA: {ref} - {ref_nt} (variant: {vid})",
            os.EX_DATAERR,
            debug,
        )
    # enrich contig sequence with iupac character
    reader.insert_snp(IUPAC_ENCODER["".join(snps.alts() + [ref])], pos)
    if store_dictionary:  # insert snp in dictionary
        insert_snp_in_dict(chrom_snps_dict, contig, variant[7], variant[9:],samples, afidx, snps)
    return chrom_snps_dict


def _initialize_samples_dict_indels(
    indels: Indels, genotypes: List[str], samples: List[str]
) -> Dict[str, str]:
    """Build a lookup table mapping each indel allele to its carrier samples. 
    This prepares per-allele sample annotations used when logging indel 
    information.

    The function iterates over all stored indels, computes carriers for each
    allele from the genotype data, and returns a dictionary keyed by ALT
    sequence.

    Args:
        indels: Collection of `Indels` objects representing all indel alleles at 
            a position.
        genotypes: List of genotype strings for all samples at the current record.
        samples: List of sample names aligned with the genotype list.

    Returns:
        A dictionary mapping each indel ALT sequence to a comma-separated string
        of 'sample:genotype' carrier labels.
    """
    samples_dict: Dict[str, str] = {indel.alt: "" for indel in indels.items}
    for indel in indels.items:
        samples_dict[indel.alt] = _retrieve_carriers(genotypes, samples, str(indel.gtidx))
    return samples_dict


def _insert_indel(reader: GenomeReader, indel: str, pos: int, offset: int, indels_set: IndelsSet) -> Tuple[IndelPair, IndelInfo]:
    """Insert a single indel into the enriched contig sequence and register it. 
    This both updates the synthetic sequence and records bookkeeping information 
    for later use.

    The function delegates to the genome reader to construct the reference and
    indel-flanked sequences, pushes the new indel into the `IndelsSet`, and
    returns both the sequence pair and its assigned metadata.

    Args:
        reader: GenomeReader instance holding the contig sequence to be enriched.
        indel: Alternate indel sequence to be inserted.
        pos: Zero-based genomic position at which the indel is anchored.
        offset: Length of the reference allele being replaced.
        indels_set: IndelsSet collection that tracks all inserted indel sequences.

    Returns:
        A tuple containing the `IndelPair` with reference and indel sequences,
        and the `IndelInfo` describing the indel's synthetic coordinates.
    """
    # insert indel in reference sequence
    indel_pair = reader.insert_indel(indel, pos, offset)
    # compute indel info (fake start/stop, idx)
    indel_info = indels_set.push(indel_pair.indelseq)
    return indel_pair, indel_info

def _compute_indel_coordinates(ref: str, pos: int) -> Tuple[int, int]:
    """Compute the flanking coordinate window used to describe an indel. 
    This expands around the variant site by a fixed offset on both sides.

    The function subtracts and adds a constant number of bases around the
    provided position, taking the reference allele length into account to define
    the stop coordinate.

    Args:
        ref: Reference allele sequence for the indel.
        pos: Zero-based genomic position at which the indel is anchored.

    Returns:
        A tuple containing the start and stop coordinates of the flanking window.
    """
    # compute start/stop coordinates for indel
    start = pos - INDELOFFSET
    stop = pos + INDELOFFSET + len(ref)
    return start, stop


def insert_indel_in_dict(
    logindels: List[List[str]],
    contig: str,
    info: str,
    indel: Indel,
    indel_info: IndelInfo,
    indel_pair: IndelPair,
    samples: str,
    afidx: int,
) -> List[List[str]]:
    """Append a fully annotated indel entry to the logging structure. This 
    captures coordinate, allele, frequency and sample information for downstream 
    reporting.

    The function derives descriptive and extended identifiers, computes a
    flanking coordinate window, retrieves the allele-specific frequency, and
    appends a tab-ready row describing the indel to the log list.

    Args:
        logindels: Accumulated list of indel log rows to be written later.
        contig: Normalized contig name on which the indel is located.
        info: INFO column string from the VCF record supplying allele-frequency 
            values.
        indel: `Indel` object describing the variant allele and its position.
        indel_info: `IndelInfo` describing the synthetic coordinates and index of 
            the inserted sequence.
        indel_pair: `IndelPair` containing reference and indel-flanked sequences.
        samples: Comma-separated 'sample:genotype' labels for carriers of this indel.
        afidx: Zero-based index of the AF entry within the semicolon-separated 
            INFO fields.

    Returns:
        The updated list of indel log rows including the new entry.
    """
    # compute indel desc and extended id
    indel_desc = f"{contig}_{indel.pos}_{indel.ref}_{indel.alt}"
    indel_start, indel_stop = _compute_indel_coordinates(indel.ref, indel.pos)
    indel_id_ext = f"{contig}_{indel_start}-{indel_stop}_{indel_info.idx}"
    # retrieve allele-specific af
    af = _retrieve_af(info, afidx, indel.gtidx)
    vid = _compute_vid(contig, indel.pos, indel.ref, indel.alt)  # compute id
    fakepos = f"{indel_info.start},{indel_info.stop}"  # position key
    # fill indels log file with current indel data
    logindels.append([indel_id_ext, samples, vid, af, indel_desc, fakepos, "".join(indel_pair.refseq)])
    return logindels


def _process_indel(
    variant: List[str],
    indels: Indels,
    contig: str,
    reader: GenomeReader,
    indels_set: IndelsSet,
    logindels: List[List[str]],
    samples: List[str],
    afidx: int,
    store_dictionary: bool,
    debug: bool,
) -> Tuple[IndelsSet, List[List[str]]]:
    """Apply indel alleles from a single VCF record to the enriched contig sequence.
    This both reconstructs synthetic indel sequences and optionally records detailed 
    log metadata.

    The function validates reference alleles against the FASTA, builds per-allele
    carrier information, inserts supported indels into the `IndelsSet`, and,
    when requested, appends descriptive rows to the indel log.

    Args:
        variant: Full list of VCF fields for the current record.
        indels: Collection of `Indels` objects representing all indel alleles at 
            this position.
        contig: Normalized contig name for the current record.
        reader: GenomeReader instance holding the contig sequence to be enriched.
        indels_set: IndelsSet collection tracking all inserted indel sequences for 
            the contig.
        logindels: Accumulated list of indel log rows to be written later.
        samples: List of sample names aligned with genotype fields in the VCF record.
        afidx: Zero-based index of the AF entry within the semicolon-separated 
            INFO fields.
        store_dictionary: Flag indicating whether indel metadata should be stored 
            in the log list.
        debug: Flag indicating whether to use debug-aware error handling.

    Returns:
        A tuple containing the updated `IndelsSet` and the updated list of indel 
            log rows.
    """
    # retrieve ref allele from contig sequence
    pos = indels.pos()  # indel position
    ref_nt = "".join(reader.sequence[pos : pos + len(variant[3])])
    ref = indels.ref()  # indel reference allele
    if ref != ref_nt:  # mismatch between VCF and contig FASTA data
        vid = _compute_vid(variant[0], pos, ref, ",".join(indels.alts()))
        exception_handler(
            CrispritzEnrichmentError,
            f"Mismatching REF alleles in VCF and FASTA: {ref} - {ref_nt} (variant: {vid})",
            os.EX_DATAERR,
            debug,
        )
    # TODO: no 1:1 comparison with old (shift be -1 bp)
    #TODO: remove genotyping from samples in log
    # initialize samples dictionary for indels
    samples_dict = _initialize_samples_dict_indels(indels, variant[9:], samples)
    for indel in indels.items:
        if samples_dict[indel.alt]:  # carriers found for indel
            # reconstruct indel sequence
            indel_pair, indel_info = _insert_indel(reader, indel.alt, indel.pos, len(indel.ref), indels_set)
            if store_dictionary:
                logindels = insert_indel_in_dict(
                    logindels,
                    contig,
                    variant[7],
                    indel,
                    indel_info,
                    indel_pair,
                    samples_dict[indel.alt],
                    afidx,
                )
    return indels_set, logindels


def insert_variants(
    vcfin: TextIOWrapper,
    reader: GenomeReader,
    samples: List[str],
    contig: str,
    indels_analysis: bool,
    chrom_snps_dict: Dict[str, str],
    logindels: List[List[str]],
    store_dictionary: bool,
    debug: bool,
) -> Tuple[Dict[str, str], IndelsSet, List[List[str]]]:
    """Parse a VCF stream and insert all supported variants into a contig sequence. 
    This coordinates SNP enrichment, optional indel processing and dictionary/log 
    construction for a single contig.

    The function walks over VCF records, filters by quality, classifies alleles
    into SNPs and indels, updates the in-memory enriched sequence, and returns
    both SNP dictionaries and indel tracking structures for downstream output.

    Args:
        vcfin: Open text stream for the VCF file, positioned at the first record.
        reader: GenomeReader instance holding the contig sequence to be enriched.
        samples: List of sample names defined in the VCF header.
        contig: Normalized contig name being enriched.
        indels_analysis: Flag indicating whether indel variants should be processed 
            in addition to SNPs.
        chrom_snps_dict: Dictionary storing SNP annotations keyed by 'contig,pos'.
        logindels: Accumulated list of indel log rows to be written later.
        store_dictionary: Flag indicating whether SNP/indel metadata structures 
            should be populated.
        debug: Flag indicating whether to use debug-aware error handling.

    Returns:
        A tuple containing the updated SNP dictionary, the populated `IndelsSet`
        for this contig, and the updated list of indel log rows.
    """
    # initialize indel-specific variables (updated and only used for indels)
    indels_set = IndelsSet(debug)
    # allele frequency position in info field
    afidx = -1
    for line in vcfin:  # iterate over variants
        variant = line.strip().split("\t")  # split variant in its fields
        if _skip_variant(variant[6]):  # filter != PASS
            continue
        if afidx == -1:  # retrieve AF position in info field (done once)
            afidx = _extract_af_idx(variant[7], debug)
        assert afidx > -1
        # retrieve ref, snps, and indel alleles
        snps, indels = _split_snps_indels(int(variant[1]), variant[3], variant[4])
        if snps.items:   # insert snp in contig sequence
            chrom_snps_dict = _process_snp(
                variant,
                snps,
                contig,
                reader,
                chrom_snps_dict,
                samples,
                afidx,
                store_dictionary,
                debug,
            )
        if indels_analysis and indels.items:
            indels_set, logindels = _process_indel(variant, indels, contig, reader, indels_set, logindels, samples, afidx, store_dictionary, debug)
    return chrom_snps_dict, indels_set, logindels


def save_enriched_contig(reader: GenomeReader, outdir: str, debug: bool):
    # retrieve contig fasta prefix
    prefix = os.path.splitext(os.path.basename(reader.fname))[0]
    writer = GenomeWriter(os.path.join(outdir, f"{prefix}.enriched.fa"), debug)
    writer.write(reader.header, reader.sequence_enr)  # write enriched contig sequence


def store_dictionary_json(
    chrom_snps_dict: Dict, contig: str, outdir: str, verbosity: int, debug: bool
) -> None:
    # store dictionary in json file
    fname = os.path.join(outdir, f"snps_dict_{contig}.json")
    print_verbosity(
        f"Storing SNPs on conting {contig} in JSON dictionary",
        verbosity,
        VERBOSITYLVL[3],
    )
    start = time()  # track json dumping run time
    try:
        with open(fname, mode="w") as fout:
            json.dump(chrom_snps_dict, fout)
    except Exception as e:
        exception_handler(
            CrispritzEnrichmentError,
            f"Failed JSON dump on {contig}",
            os.EX_IOERR,
            debug,
            e,
        )
    print_verbosity(
        f"Storing SNPs on conting {contig} in JSON dictionary completed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )


def save_indels_fasta(
    indelsdir: str, contig: str, indels_contig: IndelsSet, debug: bool
) -> None:
    fasta_path = os.path.join(indelsdir, f"fake{contig}.fa")
    try:
        with open(fasta_path, "w") as fout:
            fout.write(f">fake{contig}\n")
            for seq in indels_contig.sequences:
                fout.write("".join(seq))
                fout.write("\nN\n")
    except Exception as e:
        exception_handler(
            CrispritzEnrichmentError,
            f"Failed writing fake indels FASTA for contig {contig}",
            os.EX_IOERR,
            debug,
            e,
        )


def store_indels_log(
    indelsdir: str, contig: str, logindels: List[List[str]], debug: bool
) -> None:
    log_path = os.path.join(indelsdir, f"log{contig}.txt")
    header = ["CHR", "SAMPLES", "rsID", "AF", "indel", "FAKEPOS", "refseq"]
    try:
        with open(log_path, "w") as fout:
            fout.write("\t".join(header) + "\n")
            for row in logindels:
                fout.write("\t".join(map(str, row)) + "\n")
    except Exception as e:
        exception_handler(
            CrispritzEnrichmentError,
            f"Failed writing indels log for contig {contig}",
            os.EX_IOERR,
            debug,
            e,
        )


def enrich_variants(
    fasta_vcf_map: Dict[str, EnrichPair],
    contigs: List[str],
    indels_analysis: bool,
    snpsdir: str,
    indelsdir: str,
    store_dictionary: bool,
    verbosity: int,
    debug: bool,
) -> None:
    for contig in contigs:
        chrom_snps_dict: Dict[str, str] = {}
        logindels: List[List[str]] = []  # initialize variants dictionaries
        print_verbosity(f"Enriching contig {contig}", verbosity, VERBOSITYLVL[3])
        start = time()  # track enrichment running time
        reader = GenomeReader(fasta_vcf_map[contig].fasta, debug)  # type: ignore
        reader.read()  # read contig sequence
        try:
            with gzip.open(fasta_vcf_map[contig].vcf, mode="rt") as fin:  # type: ignore
                samples = retrieve_samples(fin, fasta_vcf_map[contig].vcf, debug)  # type: ignore
                # enrich contig sequences with snps and indels
                chrom_snps_dict, indels_contig, logindels = insert_variants(fin, reader, fasta_vcf_map[contig].vcf, samples, indels_analysis, chrom_snps_dict, logindels, store_dictionary, debug)  # type: ignore
                # TODO: indels
                # store enriched contig sequence
                save_enriched_contig(reader, snpsdir, debug)
                if indels_analysis:
                    save_indels_fasta(indelsdir, contig, indels_contig, debug)
                if store_dictionary:  # dump snp dictionary in json
                    store_dictionary_json(
                        chrom_snps_dict, contig, snpsdir, verbosity, debug
                    )
                    if indels_analysis:
                        store_indels_log(indelsdir, contig, logindels, debug)
        except (IOError, Exception) as e:
            exception_handler(CrispritzEnrichmentError, f"Failed parsing VCF: {fasta_vcf_map[contig].vcf}", os.EX_IOERR, debug, e)
        print_verbosity(
            f"Enrichment on contig  {contig} completed in {time() - start:.2f}s",
            verbosity,
            VERBOSITYLVL[3],
        )


def run_enrich_genome(
    fasta_vcf_map: Dict[str, EnrichPair],
    indels_analysis: bool,
    outdir: str,
    store_dictionary: bool,
    verbosity: int,
    debug: bool,
) -> None:
    # retrieve contig to enrich with variants and those without variants associated
    contigs_vcf, contigs_wo_vcf = split_contigs(fasta_vcf_map, verbosity)
    snpsdir, indelsdir = prepare_output_dir(outdir)  # prepare enrichment output folder
    # copy content of original fasta for contig without variants
    enrich_no_variants(fasta_vcf_map, contigs_wo_vcf, snpsdir, verbosity, debug)
    # enrich contig fasta with vcf variants
    enrich_variants(
        fasta_vcf_map,
        contigs_vcf,
        indels_analysis,
        snpsdir,
        indelsdir,
        store_dictionary,
        verbosity,
        debug,
    )


def enrich_genome(fastas: List[str], vcfs: List[str], process_indels: bool, store_dictionary: bool, outdir: str, verbosity: int, debug: bool) -> None:
    # construct a fasta-vcf files map
    fasta_vcf_map = construct_fasta_vcf_map(fastas, vcfs, verbosity, debug)
    print_verbosity(
        "Enriching genome with input variants", verbosity, VERBOSITYLVL[1]
    )
    start = time()  # genome enrichment start point
    run_enrich_genome(
        fasta_vcf_map, process_indels, outdir, store_dictionary, verbosity, debug
    )  # genome enrichment
    print_verbosity(
        f"Genome enrichment on {len(fasta_vcf_map)} contigs completed in "
        f"{time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[2],
    )
