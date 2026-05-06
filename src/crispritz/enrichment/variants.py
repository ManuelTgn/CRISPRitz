""" """

from typing import List, Iterator, Optional
from dataclasses import dataclass


class Snp:
    """Represent a single SNP allele from a VCF record. This stores its genomic
    position, reference and alternate bases, and the genotype index pointing to
    the allele.

    The object is immutable after creation and is used as a lightweight carrier
    when building SNP collections for enrichment and annotation.
    """

    def __init__(self, pos: int, ref: str, alt: str, gtidx: int) -> None:
        # NOTE: pos are 0-based
        self._pos = pos  # snp position
        self._ref = ref  # ref allele
        self._alt = alt  # alt allele
        self._gtidx = gtidx  # genotype index for variant

    def __repr__(self) -> str:
        return (
            f"<Snp object; pos={self._pos}, ref={self._ref!r}, "
            f"alt={self._alt!r}, gtidx={self._gtidx}>"
        )

    @property
    def pos(self) -> int:
        return self._pos

    @property
    def ref(self) -> str:
        return self._ref

    @property
    def alt(self) -> str:
        return self._alt

    @property
    def gtidx(self) -> int:
        return self._gtidx


class Snps:
    """Collect multiple SNP alleles that originate from the same VCF record.
    This groups individual `Snp` objects and exposes convenience accessors over
    the set.

    The container behaves like a lightweight sequence with iteration, truthiness,
    and helper methods to retrieve shared position, reference base and allele
    metadata across all SNPs.
    """

    def __init__(self, items: Optional[List[Snp]] = None) -> None:
        self._items: List[Snp] = items if items is not None else []

    def __repr__(self) -> str:
        return f"Snps(n={len(self._items)}, items={self._items!r})"

    def __len__(self) -> int:
        return len(self._items)

    def __iter__(self) -> Iterator[Snp]:
        return iter(self._items)

    def __bool__(self) -> bool:
        return bool(self._items)

    def add(self, snp: Snp) -> None:
        if not isinstance(snp, Snp):
            raise TypeError(f"Snps.add expects a Snp, got {type(snp).__name__}")
        self._items.append(snp)

    @property
    def items(self) -> List[Snp]:
        # expose a copy to avoid accidental external mutation
        return list(self._items)

    @property
    def alts(self) -> List[str]:
        return [s.alt for s in self._items]

    @property
    def gtidxs(self) -> List[int]:
        return [s.gtidx for s in self._items]

    @property
    def pos(self) -> int:
        # all SNPs from one record share the same position
        assert self._items
        return self._items[0].pos

    @property
    def ref(self) -> str:
        assert self._items
        return self._items[0].ref


class Indel:
    """Represent a single indel allele from a VCF record. This stores its genomic
    position, reference and alternate sequences, and the genotype index for the
    allele.

    The object also exposes a convenience check for symbolic alleles so that
    callers can consistently skip non-literal representations.

    Args:
        pos: Zero-based genomic position of the indel.
        ref: Reference allele sequence from the VCF record.
        alt: Alternate allele sequence from the VCF record.
        gtidx: Genotype index referring to this allele in sample genotypes.
    """

    def __init__(self, pos: int, ref: str, alt: str, gtidx: int) -> None:
        # NOTE: pos are 0-based
        self._pos = pos
        self._ref = ref
        self._alt = alt
        self._gtidx = gtidx

    def __repr__(self) -> str:
        return f"<Indel object; pos={self._pos}, ref={self._ref!r}, alt={self._alt!r}, gtidx={self._gtidx}>"

    @property
    def pos(self) -> int:
        return self._pos

    @property
    def ref(self) -> str:
        return self._ref

    @property
    def alt(self) -> str:
        return self._alt

    @property
    def gtidx(self) -> int:
        return self._gtidx

    def is_symbolic(self) -> bool:
        # skipped <...> alleles; this helps keep that rule centralized
        return "<" in self._ref or "<" in self._alt


class Indels:
    """Collect multiple indel alleles that originate from the same VCF record.
    This groups individual `Indel` objects and provides helpers to access shared
    coordinates and filter alleles.

    The container behaves like a lightweight sequence with iteration, truthiness
    and convenience methods for retrieving positions, reference bases and
    non-symbolic subsets of the stored indels.
    """

    def __init__(self, items: Optional[List[Indel]] = None) -> None:
        self._items: List[Indel] = items if items is not None else []

    def __repr__(self) -> str:
        return f"Indels(n={len(self._items)}, items={self._items!r})"

    def __len__(self) -> int:
        return len(self._items)

    def __iter__(self) -> Iterator[Indel]:
        return iter(self._items)

    def __bool__(self) -> bool:
        return bool(self._items)

    def add(self, indel: Indel) -> None:
        if not isinstance(indel, Indel):
            raise TypeError(f"Indels.add expects an Indel, got {type(indel).__name__}")
        self._items.append(indel)

    @property
    def items(self) -> List[Indel]:
        return list(self._items)

    def alts(self) -> List[str]:
        return [i.alt for i in self._items]

    def gtidxs(self) -> List[int]:
        return [i.gtidx for i in self._items]

    def pos(self) -> int:
        assert self._items
        return self._items[0].pos

    def ref(self) -> str:
        assert self._items
        return self._items[0].ref

    def non_symbolic(self) -> "Indels":
        return Indels([i for i in self._items if not i.is_symbolic()])


@dataclass
class IndelInfo:
    idx: int
    start: int
    stop: int


@dataclass
class IndelPair:
    refseq: List[str]
    indelseq: List[str]


class IndelsSet:

    def __init__(self, debug: bool) -> None:
        self._debug = debug  # store debug flag
        self._start = 0  # indel start index in fake fasta
        self._indel_seqs: List[List[str]] = []  # indel sequences
        self._i = 1  # indel index

    def push(self, indel_seq: List[str]) -> IndelInfo:
        self._indel_seqs.append(indel_seq)  # push indel sequence in list
        start_i, stop_i = self._start, self._start + len(
            indel_seq
        )  # compute indel stop position in fake fasta
        self._start = stop_i + 1  # update start position
        self._i += 1  # update indel id
        return IndelInfo(idx=self._i, start=start_i, stop=stop_i)

    @property
    def start_i(self) -> int:
        return self._start

    @property
    def sequences(self) -> List[List[str]]:
        return self._indel_seqs
