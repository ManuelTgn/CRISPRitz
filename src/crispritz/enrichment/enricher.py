""" """

from .crispritz_enrichment_error import CrispritzEnrichmentError

from ..exception_handlers import exception_handler
from ..genome_io import GenomeReader, GenomeWriter, INDELOFFSET
from ..dna_alphabet import IUPAC_ENCODER, IUPACTABLE
from ..utils import create_folder, find_tabix_index, set_processes
from ..verbosity import VERBOSITY_LVL, print_verbosity
from ..progress import progress_bar, progress_bar_parallel
from .enrichment_pair import EnrichPair
from .variants import Snp, Snps, Indel, Indels, IndelsSet, IndelPair, IndelInfo

from typing import List, Dict, Set, Tuple, Union
from pysam import FastaFile, VariantFile, tabix_index
from pysam.utils import SamtoolsError
from io import TextIOWrapper
from multiprocessing.pool import Pool
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
        fasta (FastaFile): An open pysam.FastaFile object from which to read
            contig names.
        debug (bool): Flag indicating whether to use debug-aware error handling.

    Returns:
        str: The normalized contig name extracted from the FASTA file.
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


def _retrieve_contig_names(
    fasta_files: List[str], verbosity: int, debug: bool
) -> Set[str]:
    """Retrieve standardized contig names from a list of FASTA files. Each FASTA
    file is expected to contain exactly one contig.

    This function reads the headers of the provided FASTA files, extracts their
    contig names, normalizes them to start with 'chr', and returns the set of
    unique contig identifiers found.

    Args:
        fasta_files (List[str]): List of paths to FASTA files representing
            genome contigs.
        verbosity (int): Verbosity level controlling printed progress information.
        debug (bool): Flag indicating whether to use debug-aware error handling.

    Returns:
        Set[str]: A set of normalized contig names extracted from the FASTA files.
    """
    # retrieve contig names for each fasta file in genome folder
    print_verbosity(
        "Retrieving contig names from FASTA in genome folder",
        verbosity,
        VERBOSITY_LVL[3],
    )
    return {_retrieve_contig_name(FastaFile(f), debug) for f in fasta_files}


def _initialize_fasta(
    fasta_vcf_map: Dict[str, EnrichPair], fasta_files: List[str], debug: bool
) -> Dict[str, EnrichPair]:
    """Populate the FASTA entries of a contig-to-file mapping. This associates each
    contig key with exactly one corresponding FASTA file.

    This function inspects each input FASTA, derives its contig name in a
    normalized 'chr' form, checks for duplicate assignments and records the file
    path in the mapping.

    Args:
        fasta_vcf_map (Dict[str,EnrichPair]): Dictionary mapping contig names
            to `EnrichPair` objects to be updated with FASTA paths.
        fasta_files (List[str]): List of FASTA file paths to register in the mapping.
        debug (bool): Flag indicating whether to use debug-aware error handling.

    Returns:
        Dict[str,EnrichPair]: The updated contig-to-file mapping including the assigned FASTA paths.
    """
    for f in fasta_files:
        contig = FastaFile(f).references[0]  # retrieve contig name
        if not contig.startswith("chr"):
            contig = f"chr{contig}"  # avoid mismatch (see 1000G)
        # multiple fasta pointing to same contig
        if fasta_vcf_map[contig].fasta:
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
        vcf_fname (str): Path to the VCF file to be indexed.
        verbosity (int): Verbosity level controlling printed progress information.
        debug (bool): Flag indicating whether to use debug-aware error handling.
    """
    if find_tabix_index(vcf_fname):  # index found, do nothing
        return
    try:  # tabix index not found, compute index
        print_verbosity(
            f"Index not found, indexing (VCF: {vcf_fname})", verbosity, VERBOSITY_LVL[3]
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
        vcf (VariantFile): An open pysam.VariantFile object from which to read
            contig names.
        debug (bool): Flag indicating whether to use debug-aware error handling.

    Returns:
        str: The normalized contig name extracted from the VCF header.
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


def _initialize_vcf(
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
        fasta_vcf_map (Dict[str, EnrichPair]): Dictionary mapping contig names to
            `EnrichPair` objects to be updated with VCF paths.
        vcf_files (List[str]): List of VCF file paths to register in the mapping.
        verbosity (int): Verbosity level controlling printed progress information.
        debug (bool): Flag indicating whether to use debug-aware error handling.

    Returns:
        Dict[str,EnrichPair]: The updated contig-to-file mapping including the
            assigned VCF paths.
    """
    for f in vcf_files:
        # retrieve vcf contig
        _tabix_index(f, verbosity, debug)
        contig = _retrieve_contig_vcf(VariantFile(f, mode="r"), debug)
        # multiple vcf pointing to same contig
        if fasta_vcf_map[contig].vcf:
            exception_handler(
                CrispritzEnrichmentError,
                f"Multiple VCF file pointing to conting {contig}: "
                f"{f} - {fasta_vcf_map[contig].vcf}",
                os.EX_DATAERR,
                debug,
            )
        fasta_vcf_map[contig].vcf = f  # assign vcf slot
    return fasta_vcf_map


def _construct_fasta_vcf_map(
    fasta_files: List[str], vcf_files: List[str], verbosity: int, debug: bool
) -> Dict[str, EnrichPair]:
    """Build a mapping between contig names and their FASTA/VCF file pairs. This
    prepares per-contig inputs needed for downstream genome enrichment.

    The function discovers contigs from the FASTA files, initializes an
    `EnrichPair` for each, then populates the mapping with corresponding FASTA
    and VCF paths.

    Args:
        fasta_files (List[str]): List of FASTA file paths representing contig sequences.
        vcf_files (List[str]): List of VCF file paths providing variant calls per contig.
        verbosity (int): Verbosity level controlling printed progress information.
        debug (bool): Flag indicating whether to use debug-aware error handling.

    Returns:
        Dict[str,EnrichPair]: A dictionary mapping normalized contig names to
            their associated `EnrichPair` instances.
    """
    # retrieve contig names in fasta from genome folder and initialize map
    fasta_vcf_map = {
        contig: EnrichPair(debug)
        for contig in _retrieve_contig_names(fasta_files, verbosity, debug)
    }
    # initialize fasta elements in the map
    fasta_vcf_map = _initialize_fasta(fasta_vcf_map, fasta_files, debug)
    # initialize vcf elements in the map
    fasta_vcf_map = _initialize_vcf(fasta_vcf_map, vcf_files, verbosity, debug)
    return fasta_vcf_map


def _split_contigs(
    fasta_vcf_map: Dict[str, EnrichPair], verbosity: int
) -> Tuple[List[str], List[str]]:
    """Split contigs into those with and without associated VCF files. This helps
    route contigs through variant-aware or copy-only enrichment workflows.

    This function inspects the contig-to-file mapping, groups contigs based on
    whether a VCF path is present, logs the group sizes, and returns the two
    resulting contig lists.

    Args:
        fasta_vcf_map (Dict[str, EnrichPair]): Mapping from contig names to
            `EnrichPair` objects containing FASTA and VCF file associations.
        verbosity (int): Verbosity level controlling printed progress information.

    Returns:
        Tuple[List[str],List[str]]: A tuple containing a list of contigs with VCFs
            and a list of contigs without VCFs.
    """
    # divide contigs with variants from those without for different processing
    print_verbosity(
        "Retrieving contigs with associated VCFs", verbosity, VERBOSITY_LVL[3]
    )
    contigs_vcf = [contig for contig, p in fasta_vcf_map.items() if p.vcf is not None]
    contigs_wo_vcf = [contig for contig, p in fasta_vcf_map.items() if p.vcf is None]
    print_verbosity(
        f"Contigs with VCFs: {len(contigs_vcf)}, contigs without VCFs: {len(contigs_wo_vcf)}",
        verbosity,
        VERBOSITY_LVL[3],
    )
    return contigs_vcf, contigs_wo_vcf


def _prepare_output_dir(outdir: str) -> Tuple[str, str]:
    """Prepare the output directory structure for enrichment results. This ensures
    separate subfolders exist for SNP and INDEL enriched genomes.

    The function creates (or reuses) the SNP and INDEL output folders inside the
    given root directory and returns their paths.

    Args:
        outdir (str): Root output directory where variant genome subfolders are created.

    Returns:
        Tuple[str,str]: A tuple containing the SNP output directory path and the
            INDEL output directory path, in that order.
    """
    # create snps and indels out directory
    return create_folder(os.path.join(outdir, SNPDIR)), create_folder(
        os.path.join(outdir, INDELSDIR)
    )


def _enrich_no_variants(
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
        fasta_vcf_map (Dict[str, EnrichPair]): Mapping from contig names to
            `EnrichPair` objects containing FASTA paths and optional VCF paths.
        contigs (List[str]): List of contig names to process that have no associated
            VCF files.
        outdir (str): Directory where enriched FASTA files will be written.
        verbosity (int): Verbosity level controlling printed progress information.
        debug (bool): Flag indicating whether to use debug-aware error handling.
    """
    for contig in contigs:  # just copy fasta without variants for enrichment
        print_verbosity(f"Enriching contig {contig}", verbosity, VERBOSITY_LVL[3])
        start = time()  # track enrichment running time
        assert not fasta_vcf_map[contig].vcf
        reader = GenomeReader(fasta_vcf_map[contig].fasta, debug)
        reader.read()  # read contig sequence
        # define output fasta filename
        prefix = os.path.splitext(os.path.basename(fasta_vcf_map[contig].fasta))[0]
        fasta_enr = os.path.join(outdir, f"{prefix}.enriched.fa")
        writer = GenomeWriter(fasta_enr, debug)  # write contig sequence
        writer.write(reader.header, reader.sequence)
        print_verbosity(
            f"Enrichment on contig  {contig} completed in {time() - start:.2f}s",
            verbosity,
            VERBOSITY_LVL[3],
        )


def _extract_samples(vcf_header: List[str], vcf_fname: str, debug: bool) -> List[str]:
    """Extract sample identifiers from a parsed VCF header line. This assumes the
    header is already split into fields and contains standard VCF columns.

    The function returns all columns after the fixed FORMAT column, or raises a
    debug-aware error if extraction fails.

    Args:
        vcf_header (List[str]): List of header fields from the '#CHROM' VCF header line.
        vcf_fname (str): Path to the VCF file being parsed, used for error reporting.
        debug (bool): Flag indicating whether to use debug-aware error handling.

    Returns:
        List[str]: A list of sample names defined in the VCF header.
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


def _retrieve_samples(vcfin: TextIOWrapper, vcf_fname: str, debug: bool) -> List[str]:
    """Read a VCF stream and extract the list of sample names. This scans header
    lines until the '#CHROM' line is found and then parses its fields.

    The function delegates to `_extract_samples` to pull out sample identifiers
    and raises a debug-aware error if the header cannot be located.

    Args:
        vcfin (TextIOWrapper
            the start.
        vcf_fname (str): Path to the VCF file being parsed, used for error reporting.
        debug (bool): Flag indicating whether to use debug-aware error handling.

    Returns:
        List[strA:  list of sample names defined in the VCF header.
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


def _skip_variant(variant_filter: str, keep: bool) -> bool:
    """Decide whether a variant should be skipped based on its FILTER status.
    This encapsulates the logic for retaining or discarding non-PASS records.

    The function respects a global keep flag that forces retention of all
    variants, otherwise it skips those whose FILTER value differs from 'PASS'.

    Args:
        variant_filter (str): FILTER field value from the VCF record.
        keep (bool): Flag indicating whether to keep all variants regardless of
            their FILTER status.

    Returns:
        bool: True if the variant should be skipped, False otherwise.
    """
    return False if keep else variant_filter != "PASS"


def _extract_af_idx(info: str, debug: bool) -> int:
    """Locate the index of the allele-frequency (AF) field in a VCF INFO string.
    This identifies which semicolon-separated entry encodes AF values.

    The function scans the INFO components, returns the position of the first
    entry starting with 'AF', and raises an error if no such entry is found.

    Args:
        info (str): The INFO column string from a VCF record.
        debug (bool): Flag indicating whether to use debug-aware error handling.

    Returns:
        int: The zero-based index of the AF entry within the semicolon-separated
            INFO fields.
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
        pos (int): Zero-based genomic position of the variant.
        ref (str): Reference allele sequence from the VCF record.
        alts (str): Comma-separated string of alternate alleles from the VCF record.

    Returns:
        Tuple[Snps,Indels]: A tuple containing a `Snps` collection and an `Indels`
            collection built from the provided alleles.
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
        chrom (str): Chromosome name from the variant record.
        pos (int): Genomic position of the variant, either as int or string.
        ref (str): Reference allele sequence.
        alt (str): Alternate allele sequence.

    Returns:
        str: A variant identifier string of the form 'chrX-pos-ref/alt'.
    """
    chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
    return f"{chrom}-{pos}-{ref}/{alt}"


def _retrieve_carriers(
    genotypes: List[str], samples: List[str], gtidx: str, indels: bool = False
) -> str:
    """Identify which samples carry a specific allele genotype. This can optionally
    suppress genotype details for indel reporting.

    The function scans all genotype entries, selects those whose leading allele
    index matches the requested value, and returns either 'sample' or
    'sample:genotype' labels depending on the indels flag.

    Args:
        genotypes (List[str]): List of genotype strings (e.g. '0/1:...') for each sample.
        samples (List[str]): List of sample names aligned with the genotype list.
        gtidx (str): Genotype index (as a string) that identifies the allele of interest.
        indels (bool): Flag indicating whether to omit genotype strings from the output.

    Returns:
        str: A comma-separated string of carrier labels, or an empty string if no
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
        info (str): The INFO column string from a VCF record.
        idx (int): Zero-based index of the AF entry within the semicolon-separated
            INFO fields.
        gtidx (int): One-based genotype index pointing to the target alternate allele.

    Returns:
        str: The allele-frequency string for the specified alternate allele.
    """
    return info.split(";")[idx][3:].split(",")[gtidx - 1]


def _create_snp_dict_entry(carriers: str, alleles: str, vid: str, af: str) -> str:
    """Assemble a compact dictionary entry string for a SNP. This encodes carrier
    samples, alleles, a variant identifier and allele frequency in a single record.

    The function conditionally prefixes the entry with carrier information and
    always includes alleles, variant ID and AF separated by semicolons.

    Args:
        carriers (str): Comma-separated 'sample:genotype' labels for carrier samples,
            or an empty string.
        alleles (str): Reference and alternate alleles encoded as 'REF,ALT'.
        vid (str): Stable variant identifier string (e.g. from `_compute_vid`).
        af (str): Allele-frequency value for the SNP.

    Returns:
        A semicolon-delimited string representing the SNP dictionary entry.
    """
    if carriers:
        return f"{carriers};{alleles};{vid};{af}"
    return f";{alleles};{vid};{af}"


def insert_snp_in_dict(
    chrom_snps_dict: Dict[str, str],
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
        chrom_snps_dict (Dict[str, str]): Dictionary storing SNP annotations
            keyed by 'contig,pos'.
        contig (str): Normalized contig name for the SNPs being recorded.
        info (str): INFO column string from the VCF record supplying allele-frequency
            values.
        genotypes (List[str]): List of genotype strings for all samples at this record.
        samples (List[str]): List of sample names aligned with the genotype list.
        afidx (int): Zero-based index of the AF entry within the semicolon-separated
            INFO fields.
        snps (Snps): Collection of `Snp` objects representing all SNP alleles at this
            position.

    Returns:
        Dict[str,str]: The updated chromosome SNP dictionary containing the new
            position entry.
    """
    snpkey = f"{contig},{snps.pos}"  # retrieve snp key
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
        variant List[str]: Full list of VCF fields for the current record.
        snps (Snps): Collection of `Snp` objects representing all SNP alleles at
            this position.
        contig (str): Normalized contig name for the current record.
        reader (GenomeReader): GenomeReader instance holding the contig sequence
            to be enriched.
        chrom_snps_dict (Dict[str, str]): Dictionary storing SNP annotations
            keyed by 'contig,pos'.
        samples (List[str]): List of sample names aligned with genotype fields in
            the VCF record.
        afidx (int): Zero-based index of the AF entry within the semicolon-separated
            INFO fields.
        store_dictionary (bool): Flag indicating whether SNP metadata should be stored
            in the dictionary.
        debug (bool): Flag indicating whether to use debug-aware error handling.

    Returns:
        Dict[str,str]: The updated chromosome SNP dictionary, potentially with a
            new entry for this position.
    """
    # retrieve ref allele from contig sequence
    pos = snps.pos - 1  # snp position
    ref_nt = reader.sequence[pos]
    ref = snps.ref  # snp reference allele
    if snps.ref not in IUPACTABLE[ref_nt]:  # mismatch between VCF and contig FASTA data
        vid = _compute_vid(contig, snps.pos, ref, ",".join(snps.alts))
        exception_handler(
            CrispritzEnrichmentError,
            f"Mismatching REF alleles in VCF and FASTA: {ref} - {ref_nt} (variant: {vid})",
            os.EX_DATAERR,
            debug,
        )
    # enrich contig sequence with iupac character
    reader.insert_snp(IUPAC_ENCODER["".join(snps.alts + [ref])], pos)
    if store_dictionary:  # insert snp in dictionary
        insert_snp_in_dict(
            chrom_snps_dict, contig, variant[7], variant[9:], samples, afidx, snps
        )
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
        indels (Indels): Collection of `Indels` objects representing all indel
            alleles at a position.
        genotypes (List[str]): List of genotype strings for all samples at the
            current record.
        samples (List[str]): List of sample names aligned with the genotype list.

    Returns:
        Dict[str,str]: A dictionary mapping each indel ALT sequence to a comma-separated
            string of 'sample:genotype' carrier labels.
    """
    samples_dict: Dict[str, str] = {indel.alt: "" for indel in indels.items}
    for indel in indels.items:
        samples_dict[indel.alt] = _retrieve_carriers(
            genotypes, samples, str(indel.gtidx)
        )
    return samples_dict


def _insert_indel(
    reader: GenomeReader, indel: str, pos: int, offset: int, indels_set: IndelsSet
) -> Tuple[IndelPair, IndelInfo]:
    """Insert a single indel into the enriched contig sequence and register it.
    This both updates the synthetic sequence and records bookkeeping information
    for later use.

    The function delegates to the genome reader to construct the reference and
    indel-flanked sequences, pushes the new indel into the `IndelsSet`, and
    returns both the sequence pair and its assigned metadata.

    Args:
        reader (GenomeReader): GenomeReader instance holding the contig sequence t
            o be enriched.
        indel (str): Alternate indel sequence to be inserted.
        pos (int): Zero-based genomic position at which the indel is anchored.
        offset (int): Length of the reference allele being replaced.
        indels_set (IndelsSet): IndelsSet collection that tracks all inserted
            indel sequences.

    Returns:
        Tuple[IndelPair,IndelInfo]: A tuple containing the `IndelPair` with reference
            and indel sequences, and the `IndelInfo` describing the indel's synthetic
            coordinates.
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
        ref (str): Reference allele sequence for the indel.
        pos (int): Zero-based genomic position at which the indel is anchored.

    Returns:
        Tuple[int,int]: A tuple containing the start and stop coordinates of the
            flanking window.
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
        logindels (List[List[str]]): Accumulated list of indel log rows to be
            written later.
        contig (str): Normalized contig name on which the indel is located.
        info (str): INFO column string from the VCF record supplying allele-frequency
            values.
        indel (Indel): `Indel` object describing the variant allele and its position.
        indel_info (IndelInfo): `IndelInfo` describing the synthetic coordinates
            and index of the inserted sequence.
        indel_pair (IndelPair): `IndelPair` containing reference and indel-flanked
            sequences.
        samples (str): Comma-separated 'sample:genotype' labels for carriers of
            this indel.
        afidx (int): Zero-based index of the AF entry within the semicolon-separated
            INFO fields.

    Returns:
        List[List[str]]: The updated list of indel log rows including the new entry.
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
    logindels.append(
        [
            indel_id_ext,
            samples,
            vid,
            af,
            indel_desc,
            fakepos,
            "".join(indel_pair.refseq),
        ]
    )
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
        variant (List[str]): Full list of VCF fields for the current record.
        indels (Indels): Collection of `Indels` objects representing all indel
            alleles at this position.
        contig (str): Normalized contig name for the current record.
        reader (GenomeReader): GenomeReader instance holding the contig sequence
            to be enriched.
        indels_set (IndelsSet): IndelsSet collection tracking all inserted indel
            sequences for the contig.
        logindels (list[List[str]]): Accumulated list of indel log rows to be written later.
        samples (List[str]): List of sample names aligned with genotype fields in
            the VCF record.
        afidx (int): Zero-based index of the AF entry within the semicolon-separated
            INFO fields.
        store_dictionary (bool): Flag indicating whether indel metadata should be stored
            in the log list.
        debug (bool): Flag indicating whether to use debug-aware error handling.

    Returns:
        Tuple[IndelsSet,List[List[str]]]: A tuple containing the updated `IndelsSet`
            and the updated list of indel log rows.
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
    # initialize samples dictionary for indels
    samples_dict = _initialize_samples_dict_indels(indels, variant[9:], samples)
    for indel in indels.items:
        if samples_dict[indel.alt]:  # carriers found for indel
            # reconstruct indel sequence
            indel_pair, indel_info = _insert_indel(
                reader, indel.alt, indel.pos, len(indel.ref), indels_set
            )
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


def _insert_variants(
    vcfin: TextIOWrapper,
    reader: GenomeReader,
    samples: List[str],
    contig: str,
    keep: bool,
    indels_analysis: bool,
    chrom_snps_dict: Dict[str, str],
    logindels: List[List[str]],
    store_dictionary: bool,
    debug: bool,
) -> Tuple[Dict[str, str], IndelsSet, List[List[str]]]:
    """Stream all variants from a VCF handle and insert them into a contig sequence.
    This coordinates SNP and optional indel handling while maintaining per-contig
    dictionaries and logs.

    The function iterates over VCF records, applies optional FILTER-based skipping,
    lazily discovers the AF field index, dispatches SNPs and indels to their
    respective processors, and returns the updated SNP dictionary, indel set and
    indel log.

    Args:
        vcfin (TextIOWrapper): Open, text-mode VCF file handle positioned at the
            first variant record.
        reader (GenomeReader): GenomeReader instance holding the contig sequence to
            be enriched.
        samples (List[str]): List of sample names aligned with genotype fields in
            the VCF records.
        contig (str): Normalized contig name associated with the VCF records.
        keep (bool): Flag indicating whether to keep all variants (True) or skip
            non-PASS variants based on the FILTER field (False).
        indels_analysis (bool): Flag indicating whether indel variants should be
            processed in addition to SNPs.
        chrom_snps_dict (Dict[str,str]): Dictionary storing SNP annotations keyed by
            'contig,pos' that will be updated in place.
        logindels (List[List[str]]): Accumulated list of indel log rows to be updated.
        store_dictionary (bool): Flag indicating whether SNP and indel metadata
            should be recorded into the dictionary and log structures.
        debug (bool): Flag indicating whether to use debug-aware error handling.

    Returns:
        Tuple[Dict[str,str],IndelsSet,List[List[str]]]: A tuple containing the updated
            SNP dictionary, the per-contig `IndelsSet`, and the updated indel log.
    """
    # initialize indel-specific variables (updated and only used for indels)
    indels_set = IndelsSet(debug)
    # allele frequency position in info field
    afidx = -1
    for line in vcfin:  # iterate over variants
        variant = line.strip().split("\t")  # split variant in its fields
        if _skip_variant(variant[6], keep):  # filter != PASS
            continue
        if afidx == -1:  # retrieve AF position in info field (done once)
            afidx = _extract_af_idx(variant[7], debug)
        assert afidx > -1
        # retrieve ref, snps, and indel alleles
        snps, indels = _split_snps_indels(int(variant[1]), variant[3], variant[4])
        if snps.items:  # insert snp in contig sequence
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
            indels_set, logindels = _process_indel(
                variant,
                indels,
                contig,
                reader,
                indels_set,
                logindels,
                samples,
                afidx,
                store_dictionary,
                debug,
            )
    return chrom_snps_dict, indels_set, logindels


def _save_enriched_contig(reader: GenomeReader, outdir: str, debug: bool) -> None:
    """Write an enriched contig sequence to a FASTA file on disk. This function
    derives a file name from the original contig FASTA and delegates the actual
    writing to a `GenomeWriter`.

    The function preserves the original contig header, appends an `.enriched.fa`
    suffix to the base FASTA name, and writes the in-memory enriched sequence
    so it can be used for downstream analyses.

    Args:
        reader (GenomeReader): GenomeReader instance providing the contig header
            and enriched sequence to be written.
        outdir (str): Output directory where the enriched FASTA file will be created.
        debug (bool): Flag indicating whether to use debug-aware error handling during
            file writing.
    """
    # retrieve contig fasta prefix
    prefix = os.path.splitext(os.path.basename(reader.fname))[0]
    writer = GenomeWriter(os.path.join(outdir, f"{prefix}.enriched.fa"), debug)
    writer.write(reader.header, reader.sequence_enr)  # write enriched contig sequence


def _store_dictionary_json(
    chrom_snps_dict: Dict[str, str],
    contig: str,
    outdir: str,
    verbosity: int,
    debug: bool,
) -> None:
    """Persist the SNP dictionary for a contig as a JSON file on disk. This enables
    downstream tools to reload compact SNP annotations without reparsing the VCF.

    The function builds a contig-specific JSON path, reports progress according
    to the verbosity level, attempts to dump the dictionary, and delegates error
    handling to the enrichment exception handler.

    Args:
        chrom_snps_dict (Dict[str,str]): Dictionary storing SNP annotations keyed
            by 'contig,pos'.
        contig (str): Normalized contig name whose SNP dictionary is being saved.
        outdir (str): Output directory where the JSON file will be written.
        verbosity (int): Verbosity level controlling progress messages.
        debug (bool): Flag indicating whether to include debug details in error
            handling.
    """
    # store dictionary in json file
    fname = os.path.join(outdir, f"snps_dict_{contig}.json")
    print_verbosity(
        f"Storing SNPs on conting {contig} in JSON dictionary",
        verbosity,
        VERBOSITY_LVL[3],
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
        VERBOSITY_LVL[3],
    )


def _save_indels_fasta(
    indelsdir: str, contig: str, indels_contig: IndelsSet, debug: bool
) -> None:
    """Write all synthetic indel sequences for a contig to a FASTA file. This
    produces a dedicated 'fake' contig file that aggregates every reconstructed
    indel sequence for downstream analyses.

    The function opens a contig-specific FASTA path, writes a single header,
    streams each indel sequence separated by sentinel 'N' lines, and reports
    any I/O failures through the enrichment exception handler.

    Args:
        indelsdir (str): Output directory where the indel FASTA file will be written.
        contig (str): Normalized contig name whose indel sequences are being exported.
        indels_contig (IndelSet): `IndelsSet` collection providing the synthetic
            indel sequences.
        debug (bool): Flag indicating whether to include debug details in error handling.
    """
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


def _store_indels_log(
    indelsdir: str, contig: str, logindels: List[List[str]], debug: bool
) -> None:
    """Write the collected indel log entries for a contig to a tab-delimited file.
    This creates a human-readable summary of indel coordinates, alleles and
    carrier information for downstream inspection.

    The function opens a contig-specific log path, writes a fixed header, streams
    each log row as a tab-separated line, and reports any I/O failures through
    the enrichment exception handler.

    Args:
        indelsdir (str): Output directory where the indel log file will be written.
        contig (str): Normalized contig name whose indel log is being exported.
        logindels (List[List[str]]): Accumulated list of indel log rows to be serialized.
        debug (bool): Flag indicating whether to include debug details in error handling.
    """
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


def _enrich_variants_worker(
    args: Tuple[str, str, str, bool, str, str, bool, bool, int, bool],
) -> None:
    """Execute variant enrichment for a single contig in an isolated worker.
    This function is intended to be run in parallel across contigs.

    The function loads the contig FASTA, parses the corresponding VCF to insert
    SNPs and optional indels, writes enriched sequences and auxiliary outputs to
    disk, and reports any parsing or I/O errors through the enrichment exception
    handler.

    Args:
        args (Tuple[str,str,str,bool,str,str,bool,bool,int,bool]): Tuple bundling the
            contig identifier, input FASTA and VCF paths, keep flag, SNP and indel
            output directories, analysis and dictionary flags, verbosity level and
            debug flag.
    """
    (
        contig,
        fasta,
        vcf,
        keep,
        snpsdir,
        indelsdir,
        indels_analysis,
        store_dictionary,
        verbosity,
        debug,
    ) = args
    chrom_snps_dict: Dict[str, str] = {}
    logindels: List[List[str]] = []  # initialize variants dictionaries
    print_verbosity(f"Enriching contig {contig}", verbosity, VERBOSITY_LVL[3])
    start = time()  # track enrichment running time
    reader = GenomeReader(fasta, debug)
    reader.read()  # read contig sequence
    try:
        with gzip.open(vcf, mode="rt") as fin:
            samples = _retrieve_samples(fin, vcf, debug)
            # enrich contig sequences with snps and indels
            chrom_snps_dict, indels_contig, logindels = _insert_variants(
                fin,
                reader,
                samples,
                contig,
                keep,
                indels_analysis,
                chrom_snps_dict,
                logindels,
                store_dictionary,
                debug,
            )
            # store enriched contig sequence
            _save_enriched_contig(reader, snpsdir, debug)
            if indels_analysis:
                _save_indels_fasta(indelsdir, contig, indels_contig, debug)
            if store_dictionary:  # dump snp dictionary in json
                _store_dictionary_json(
                    chrom_snps_dict, contig, snpsdir, verbosity, debug
                )
                if indels_analysis:
                    _store_indels_log(indelsdir, contig, logindels, debug)
    except (IOError, Exception) as e:
        exception_handler(
            CrispritzEnrichmentError,
            f"Failed parsing VCF: {vcf}",
            os.EX_IOERR,
            debug,
            e,
        )
    print_verbosity(
        f"Enrichment on contig  {contig} completed in {time() - start:.2f}s",
        verbosity,
        VERBOSITY_LVL[3],
    )


def _enrich_variants(
    fasta_vcf_map: Dict[str, EnrichPair],
    contigs: List[str],
    keep: bool,
    indels_analysis: bool,
    snpsdir: str,
    indelsdir: str,
    store_dictionary: bool,
    threads: int,
    verbosity: int,
    debug: bool,
) -> None:
    """Coordinate per-contig variant enrichment across a set of contigs.
    This function prepares worker tasks and runs them either serially or in
    parallel to generate enriched contig FASTA files and auxiliary outputs.

    The function builds a list of enrichment jobs from the FASTA/VCF mapping,
    chooses between a simple loop and a multiprocessing pool based on the
    requested thread count, and tracks progress using the configured verbosity
    settings.

    Args:
        fasta_vcf_map (Dict[str,EnrichPair]): Mapping from contig name to paired
            ASTA/VCF inputs.
        contigs (List[str]): List of contig identifiers to be enriched.
        keep (bool): Flag indicating whether to keep all variants, including those
            failing FILTER.
        indels_analysis (bool): Flag indicating whether indel variants should be
            processed in addition to SNPs.
        snpsdir (str): Output directory where SNP-enriched contig FASTA files are written.
        indelsdir (str): Output directory where indel-related outputs are written.
        store_dictionary (bool): Flag indicating whether SNP/indel metadata dictionaries
            and logs should be generated.
        threads (int): Maximum number of worker processes to use for enrichment.
        verbosity (int): Verbosity level controlling progress messages.
        debug (bool): Flag indicating whether to use debug-aware error handling.
    """
    # initialize tasks for genome enrichment
    tasks = [
        (
            contig,
            fasta_vcf_map[contig].fasta,
            fasta_vcf_map[contig].vcf,
            keep,
            snpsdir,
            indelsdir,
            indels_analysis,
            store_dictionary,
            verbosity,
            debug,
        )
        for contig in contigs
    ]
    if threads == 1 or len(tasks) == 1:  # only one task to process
        for t in progress_bar(tasks, "Enriched contigs", verbosity):
            _enrich_variants_worker(t)
    else:  # use multiprocessing
        with Pool(processes=set_processes(len(tasks), threads)) as pool:
            with progress_bar_parallel(
                len(tasks), "Enriched contigs", verbosity
            ) as pbar:
                for _ in pool.imap_unordered(_enrich_variants_worker, tasks):
                    pbar.update(1)


def _run_enrich_genome(
    fasta_vcf_map: Dict[str, EnrichPair],
    keep: bool,
    indels_analysis: bool,
    outdir: str,
    store_dictionary: bool,
    threads: int,
    verbosity: int,
    debug: bool,
) -> None:
    """Run the internal genome enrichment workflow for all contigs in the map.
    This orchestrates directory setup, handling of contigs with and without
    variants, and parallel execution of per-contig enrichment.

    The function splits contigs based on VCF availability, prepares output
    folders, copies unchanged FASTA sequences for contigs lacking variants, and
    delegates variant-aware enrichment to `_enrich_variants` with the requested
    threading and logging configuration.

    Args:
        fasta_vcf_map (Dict[str,EnrichPair]): Mapping from contig name to paired
            FASTA/VCF inputs.
        keep (bool): Flag indicating whether to keep all variants, including those
            failing FILTER.
        indels_analysis (bool): Flag indicating whether indel variants should be
            processed in addition to SNPs.
        outdir (str): Base output directory where enrichment results will be created.
        store_dictionary (bool): Flag indicating whether SNP/indel metadata dictionaries
            and logs should be generated.
        threads (int): Maximum number of worker processes to use for enrichment.
        verbosity (int): Verbosity level controlling progress messages.
        debug (bool): Flag indicating whether to use debug-aware error handling.
    """
    # retrieve contig to enrich with variants and those without variants associated
    contigs_vcf, contigs_wo_vcf = _split_contigs(fasta_vcf_map, verbosity)
    snpsdir, indelsdir = _prepare_output_dir(outdir)  # prepare enrichment output folder
    # copy content of original fasta for contig without variants
    _enrich_no_variants(fasta_vcf_map, contigs_wo_vcf, snpsdir, verbosity, debug)
    # enrich contig fasta with vcf variants
    _enrich_variants(
        fasta_vcf_map,
        contigs_vcf,
        keep,
        indels_analysis,
        snpsdir,
        indelsdir,
        store_dictionary,
        threads,
        verbosity,
        debug,
    )


def enrich_genome(
    fastas: List[str],
    vcfs: List[str],
    keep: bool,
    process_indels: bool,
    store_dictionary: bool,
    outdir: str,
    threads: int,
    verbosity: int,
    debug: bool,
) -> None:
    """Top-level entry point to enrich an entire genome with variants from VCF files.
    This function wires together input discovery, per-contig enrichment and high-level
    progress reporting.

    The function builds a FASTA/VCF mapping from the provided file lists, logs the
    start and end of the genome-wide workflow, and delegates detailed processing to
    `_run_enrich_genome` using the requested keep, indel, dictionary and threading
    options.

    Args:
        fastas (List[str]): List of input contig FASTA file paths.
        vcfs (List[str]): List of input VCF file paths containing variants.
        keep (bool): Flag indicating whether to keep all variants, including those
            failing FILTER.
        process_indels (bool): Flag indicating whether indel variants should be
            processed in addition to SNPs.
        store_dictionary (bool): Flag indicating whether SNP/indel metadata dictionaries
            and logs should be generated.
        outdir (str): Base output directory where enrichment results will be created.
        threads (int): Maximum number of worker processes to use for enrichment.
        verbosity (int): Verbosity level controlling progress messages.
        debug (bool): Flag indicating whether to use debug-aware error handling.
    """
    # construct a fasta-vcf files map
    fasta_vcf_map = _construct_fasta_vcf_map(fastas, vcfs, verbosity, debug)
    print_verbosity("Enriching genome with input variants", verbosity, VERBOSITY_LVL[1])
    start = time()  # genome enrichment start point
    _run_enrich_genome(
        fasta_vcf_map,
        keep,
        process_indels,
        outdir,
        store_dictionary,
        threads,
        verbosity,
        debug,
    )  # genome enrichment
    print_verbosity(
        f"Genome enrichment on {len(fasta_vcf_map)} contigs completed in "
        f"{time() - start:.2f}s",
        verbosity,
        VERBOSITY_LVL[2],
    )
