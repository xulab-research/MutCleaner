# mutcleaner/core/codon_table.py
from __future__ import annotations


from collections.abc import Collection
from typing import TYPE_CHECKING

from .constants import (
    STANDARD_GENETIC_CODE_DNA,
    STANDARD_START_CODONS_DNA,
    STANDARD_GENETIC_CODE_RNA,
    STANDARD_START_CODONS_RNA,
)

if TYPE_CHECKING:
    from typing import Dict, List, Literal, Optional

__all__ = ["CodonTable"]


def __dir__() -> List[str]:
    return __all__


class CodonTable:
    """codon table used to translate codons to amino acids"""

    def __init__(
        self,
        name: str,
        codon_map: Dict[str, str],
        start_codons: Optional[Collection[str]] = None,
        stop_codons: Optional[Collection[str]] = None,
    ):
        self.name = name
        self.codon_map = {k.upper(): v for k, v in codon_map.items()}

        # auto detect stop codons
        if stop_codons is None:
            self.stop_codons = [codon for codon, aa in self.codon_map.items() if aa == "*"]
        else:
            self.stop_codons = set([c.upper() for c in stop_codons])

        # set start codons
        if start_codons is None:
            self.start_codons = STANDARD_START_CODONS_DNA
        else:
            self.start_codons = set([c.upper() for c in start_codons])

    def translate_codon(self, codon: str) -> str:
        """translate single codon to corresponding amino acid"""
        return self.codon_map.get(codon.upper(), "X")

    def is_stop_codon(self, codon: str) -> bool:
        """check if codon is a stop codon"""
        return codon.upper() in self.stop_codons

    def is_start_codon(self, codon: str) -> bool:
        """check if codon is a start codon"""
        return codon.upper() in self.start_codons

    @classmethod
    def get_standard_table(cls, seq_type: Literal["DNA", "RNA"] = "DNA") -> "CodonTable":
        """get standard codon table (NCBI standard)"""
        if seq_type == "DNA":
            return cls("Standard", STANDARD_GENETIC_CODE_DNA, STANDARD_START_CODONS_DNA)
        elif seq_type == "RNA":
            return cls("Standard", STANDARD_GENETIC_CODE_RNA, STANDARD_START_CODONS_RNA)
        else:
            raise ValueError("Invalid sequence type")

    @classmethod
    def get_table_by_name(cls, name: str, seq_type: Literal["DNA", "RNA"] = "DNA") -> "CodonTable":
        """get codon table by name"""
        return {
            "Standard": cls.get_standard_table(seq_type),
        }[name]
