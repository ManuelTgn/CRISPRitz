""" """

from .crispritz_error import CrispritzEnrichmentError, EnrichmentPairError
from .exception_handlers import exception_handler
from .genome_io import GenomeReader, GenomeWriter, INDELOFFSET
from .dna_alphabet import IUPAC_ENCODER, IUPACTABLE
from .crispritz_argparse import CrispritzEnrichmentInputArgs
from .utils import print_verbosity, create_folder, find_tabix_index, VERBOSITYLVL

from typing import List, Dict, Set, Optional, Tuple, Union
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


@dataclass
class Snps:
    snps: List[str]
    gtidxs: List[str]


@dataclass
class Indels:
    indels: List[str]
    gtidxs: List[str]


class IndelsSet:

    def __init__(self, debug: bool) -> None:
        self._debug = debug  # store debug flag
        self._start = 0  # indel start index in fake fasta
        self._indel_seqs: List[List[str]] = []  # indel sequences
        self._i = 1  # indel index

    def push(self, indel_seq: List[str]) -> Tuple[int, int, int]:
        self._indel_seqs.append(indel_seq)  # push indel sequence in list
        start_i, stop_i = self._start, self._start + len(
            indel_seq
        )  # compute indel stop position in fake fasta
        self._start = stop_i + 1  # update start position
        self._i += 1  # update indel id
        return self._i - 1, start_i, stop_i

    @property
    def start_i(self) -> int:
        return self._start

    @property
    def sequences(self) -> List[List[str]]:
        return self._indel_seqs


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


def _tabix_index(vcf_fname: str, verbosity: int, debug: bool) -> None:
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


def _extract_af(info: str, vcf_fname: str, debug: bool) -> int:
    for i, e in enumerate(info.split(";")):  # look in INFO field
        if e[:2] == "AF":
            return i
    exception_handler(
        CrispritzEnrichmentError,
        f"Failed retrieving AF index from VCF: {vcf_fname}",
        os.EX_IOERR,
        debug,
    )


def _skip_variant(variant_filter: str) -> bool:
    return variant_filter != "PASS"


def _retrieve_af(info: str, idx: int, gtidx: int) -> str:
    return info.split(";")[idx][3:].split(",")[gtidx - 1]


def _compute_vid(chrom: str, pos: Union[int, str], ref: str, alt: str) -> str:
    chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
    return f"{chrom}-{pos}-{ref}/{alt}"


def _retrieve_carriers(genotypes: List[str], samples: List[str], gtidx: str) -> str:
    carriers = [
        f"{samples[i]}:{g}"
        for i, gt in enumerate(genotypes)
        if gtidx in (g := gt.split(":")[0])
    ]
    return ",".join(sorted(carriers))


def _create_snp_dict_entry(carriers: str, alleles: str, vid: str, af: str) -> str:
    if carriers:
        return f"{carriers};{alleles};{vid};{af}"
    return f";{alleles};{vid};{af}"


def insert_snp_in_dict(
    chrom_snps_dict: Dict,
    contig: str,
    variant: List[str],
    samples: List[str],
    afidx: int,
    snps: Snps,
) -> Dict[str, str]:
    snpkey = f"{contig},{variant[1]}"  # retrieve snp key
    # compute dictionary entry for each snp (multiallelic sites)
    entries = []
    for snp, gtidx in zip(snps.snps, snps.gtidxs):
        carriers = _retrieve_carriers(variant[9:], samples, gtidx)
        af = _retrieve_af(variant[7], afidx, int(gtidx))
        alleles = f"{variant[3]},{snp}"
        vid = _compute_vid(contig, variant[1], variant[3], snp)
        entries.append(_create_snp_dict_entry(carriers, alleles, vid, af))
    chrom_snps_dict[snpkey] = "$".join(entries) if len(entries) > 1 else entries[0]
    return chrom_snps_dict


def _process_snp(
    variant: List[str],
    snps: Snps,
    reader: GenomeReader,
    chrom_snps_dict: Dict[str, str],
    samples: List[str],
    afidx: int,
    store_dictionary: bool,
    debug: bool,
) -> Dict[str, str]:
    # retrieve variant contig and position
    contig = variant[0] if variant[0].startswith("chr") else f"chr{variant[0]}"
    pos = int(variant[1]) - 1
    ref, ref_nt = (
        variant[3],
        reader.sequence[pos],
    )  # retrieve ref allele from contig sequence
    if ref not in IUPACTABLE[ref_nt]:  # mismatch between VCF and contig FASTA data
        variantid = _compute_vid(contig, pos, ref, ",".join(snps.snps))
        exception_handler(
            CrispritzEnrichmentError,
            f"Mismatching REF alleles in VCF and FASTA: {ref} - {ref_nt} (variant: {variantid})",
            os.EX_DATAERR,
            debug,
        )
    # enrich contig sequence with iupac character
    reader.insert_snp(IUPAC_ENCODER["".join(snps.snps + [ref])], pos)
    if store_dictionary:  # insert snp in dictionary
        insert_snp_in_dict(chrom_snps_dict, contig, variant, samples, afidx, snps)
    return chrom_snps_dict


def _split_snps_indels(ref: str, alts: str) -> Tuple[Snps, Indels]:
    # retrieve reference, snps and indels for current variant
    snps, indels, snps_gtidxs, indels_gtidxs = [], [], [], []
    for i, alt in enumerate(alts.strip().split(",")):
        if len(alt) == len(ref) == 1:  # snp found
            snps.append(alt)
            snps_gtidxs.append(str(i + 1))
        else:  # indel found
            indels.append(alt)
            indels_gtidxs.append(str(i + 1))
    return Snps(snps=snps, gtidxs=snps_gtidxs), Indels(
        indels=indels, gtidxs=indels_gtidxs
    )


def _initialize_samples_dict_indels(
    indels: Indels, genotypes: List[str], samples: List[str]
) -> Dict[str, str]:
    samples_dict: Dict[str, str] = {indel: "" for indel in indels.indels}
    for i, indel in enumerate(indels.indels):
        samples_dict[indel] = _retrieve_carriers(genotypes, samples, indels.gtidxs[i])
    return samples_dict


def _compute_indel_coordinates(ref: str, pos: int) -> Tuple[int, int]:
    # compute start/stop coordinates for indel
    start = pos - INDELOFFSET
    stop = pos + INDELOFFSET + len(ref)
    return start, stop


def insert_indel_in_dict(
    logindels: List[List[str]],
    variant: List[str],
    indel: str,
    contig: str,
    gtidx: str,
    start_i: int,
    stop_i: int,
    indel_idx: int,
    samples: str,
    vid: str,
    afidx: int,
    refseq: str,
) -> List[List[str]]:
    # compute indel info and extended id
    indel_info = f"{contig}_{variant[1]}_{variant[3]}_{indel}"
    indel_start, indel_stop = _compute_indel_coordinates(variant[3], int(variant[1]))
    indel_id_ext = f"{contig}_{indel_start}-{indel_stop}_{indel_idx}"
    # retrieve allele-specific af
    af = _retrieve_af(variant[7], afidx, int(gtidx))
    fakepos = f"{start_i},{stop_i}"  # position key
    # fill indels log file with current indel data
    logindels.append([indel_id_ext, samples, vid, af, indel_info, fakepos, refseq])
    return logindels


def _process_indel(
    variant: List[str],
    indels: Indels,
    reader: GenomeReader,
    indels_contig: IndelsSet,
    logindels: List[List[str]],
    samples: List[str],
    afidx: int,
    store_dictionary: bool,
    debug: bool,
) -> Tuple[IndelsSet, List[List[str]]]:
    # retrieve variant contig and position
    contig = variant[0] if variant[0].startswith("chr") else f"chr{variant[0]}"
    pos = int(variant[1]) - 1
    # retrieve ref allele from contig sequence
    ref, ref_nt = variant[3], "".join(reader.sequence[pos : pos + len(variant[3])])
    if ref not in IUPACTABLE[ref_nt]:  # mismatch between VCF and contig FASTA data
        variantid = _compute_vid(variant[0], pos, ref, ",".join(indels.indels))
        exception_handler(
            CrispritzEnrichmentError,
            f"Mismatching REF alleles in VCF and FASTA: {ref} - {ref_nt} (variant: {variantid})",
            os.EX_DATAERR,
            debug,
        )
    if "<" in ref:
        return indels_contig, logindels  # TODO: add return
    # initialize samples dictionary for indels
    samples_dict = _initialize_samples_dict_indels(indels, variant[9:], samples)
    for i, indel in enumerate(indels.indels):
        if samples_dict[indel]:  # carriers found for indel
            # reconstruct indel sequence
            indel_pair = reader.insert_indel(ref, list(indel), pos)
            indel_i, start_i, stop_i = indels_contig.push(indel_pair.indelseq)
            # compute variant id
            vid = _compute_vid(variant[0], pos, ref, indel)
            if store_dictionary:
                logindels = insert_indel_in_dict(
                    logindels,
                    variant,
                    indel,
                    contig,
                    indels.gtidxs[i],
                    start_i,
                    stop_i,
                    indel_i,
                    samples_dict[indel],
                    vid,
                    afidx,
                    "".join(indel_pair.refseq),
                )
    return indels_contig, logindels


def insert_variants(
    vcfin: TextIOWrapper,
    reader: GenomeReader,
    vcf_fname: str,
    samples: List[str],
    indels_analysis: bool,
    chrom_snps_dict: Dict[str, str],
    logindels: List[List[str]],
    store_dictionary: bool,
    debug: bool,
) -> Tuple[Dict[str, str], IndelsSet, List[List[str]]]:
    # initialize indel-specific variables (updated and only used for indels)
    indels_contig = IndelsSet(debug)
    # allele frequency position in info field
    afidx = -1
    for line in vcfin:  # iterate over variants
        variant = line.strip().split("\t")  # split variant in its fields
        if _skip_variant(variant[6]):  # filter != PASS
            continue
        if afidx == -1:  # retrieve AF position in info field
            afidx = _extract_af(variant[7], vcf_fname, debug)
        assert afidx >= -1
        # retrieve ref, snps, and indel alleles
        snps, indels = _split_snps_indels(variant[3], variant[4])
        if snps.snps:  # insert snp in contig sequence
            chrom_snps_dict = _process_snp(
                variant,
                snps,
                reader,
                chrom_snps_dict,
                samples,
                afidx,
                store_dictionary,
                debug,
            )
        if indels_analysis and indels.indels:
            indels_contig, logindels = _process_indel(
                variant,
                indels,
                reader,
                indels_contig,
                logindels,
                samples,
                afidx,
                store_dictionary,
                debug,
            )
    return chrom_snps_dict, indels_contig, logindels


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
        print_verbosity(
            f"Enrichment on contig  {contig} completed in {time() - start:.2f}s",
            verbosity,
            VERBOSITYLVL[3],
        )


def enrich_genome(
    fasta_vcf_map: Dict[str, EnrichPair],
    indels_analysis: bool,
    outdir: str,
    store_dictionary: bool,
    verbosity: int,
    debug: bool,
):
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


def add_variants(args: CrispritzEnrichmentInputArgs) -> None:
    # construct a fasta-vcf files map
    fasta_vcf_map = construct_fasta_vcf_map(
        args.fastas, args.vcfs, args.verbosity, args.debug
    )
    start = time()  # genome enrichment start point
    print_verbosity(
        "Enriching genome with input variants", args.verbosity, VERBOSITYLVL[1]
    )
    enrich_genome(
        fasta_vcf_map, args.indels, args.outdir, True, args.verbosity, args.debug
    )  # genome enrichment
    print_verbosity(
        f"Genome enrichment on {len(fasta_vcf_map)} contigs completed in "
        f"{time() - start:.2f}s",
        args.verbosity,
        VERBOSITYLVL[2],
    )
