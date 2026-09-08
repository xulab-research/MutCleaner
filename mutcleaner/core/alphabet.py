# mutcleaner/core/alphabet.py
from __future__ import annotations

from abc import ABC
from typing import TYPE_CHECKING

from .constants import (
    STANDARD_DNA_BASES,
    AMBIGUOUS_DNA_BASES,
    STANDARD_RNA_BASES,
    AMBIGUOUS_RNA_BASES,
    STANDARD_AMINO_ACIDS,
    AMBIGUOUS_AMINO_ACIDS,
    AA1_TO_3,
    AA3_TO_1,
)

if TYPE_CHECKING:
    from typing import List, Set

__all__ = ["BaseAlphabet", "DNAAlphabet", "RNAAlphabet", "ProteinAlphabet"]


def __dir__() -> List[str]:
    return __all__


class BaseAlphabet(ABC):
    """Base class for biological alphabets"""

    def __init__(self, letters: Set[str], name: str):
        self.letters = set(letter.upper() for letter in letters)
        self.name = name

    def is_valid_char(self, char: str) -> bool:
        """Check if character is valid in this alphabet"""
        return char.upper() in self.letters

    def is_valid_sequence(self, sequence: str) -> bool:
        """Check if entire sequence is valid"""
        return all(self.is_valid_char(char) for char in sequence)

    def get_invalid_chars(self, sequence: str) -> Set[str]:
        """Get set of invalid characters in sequence"""
        return set(char.upper() for char in sequence) - self.letters

    def validate_sequence(self, sequence: str) -> str:
        """Validate sequence and raise error if invalid"""
        invalid = self.get_invalid_chars(sequence)
        if invalid:
            raise ValueError(f"Invalid characters in {self.name} sequence: {invalid}")
        return sequence.upper()

    def __contains__(self, char: str) -> bool:
        return self.is_valid_char(char)

    def __str__(self) -> str:
        return f"{self.name}Alphabet: {''.join(sorted(self.letters))}"


class DNAAlphabet(BaseAlphabet):
    """DNA alphabet (A, T, C, G)"""

    def __init__(self, include_ambiguous: bool = False):
        standard = STANDARD_DNA_BASES
        if include_ambiguous:
            # IUPAC ambiguous nucleotide codes
            ambiguous = AMBIGUOUS_DNA_BASES
            letters = standard | ambiguous
        else:
            letters = standard

        super().__init__(letters, "DNA")
        self.include_ambiguous = include_ambiguous


class RNAAlphabet(BaseAlphabet):
    """RNA alphabet (A, U, C, G)"""

    def __init__(self, include_ambiguous: bool = False):
        standard = STANDARD_RNA_BASES
        if include_ambiguous:
            ambiguous = AMBIGUOUS_RNA_BASES
            letters = standard | ambiguous
        else:
            letters = standard

        super().__init__(letters, "RNA")
        self.include_ambiguous = include_ambiguous


class ProteinAlphabet(BaseAlphabet):
    """Protein alphabet (20 standard amino acids + stop codon)"""

    def __init__(self, include_stop: bool = True, include_ambiguous: bool = False):
        # 20 standard amino acids
        standard = STANDARD_AMINO_ACIDS

        letters = standard.copy()

        if include_stop:
            letters.add("*")

        if include_ambiguous:
            # Ambiguous amino acids
            letters.update(AMBIGUOUS_AMINO_ACIDS)

        super().__init__(letters, "Protein")
        self.include_stop = include_stop
        self.include_ambiguous = include_ambiguous

    def get_three_letter_code(self, one_letter: str, strict: bool = True) -> str:
        """Convert one-letter to three-letter amino acid code"""
        if strict:
            if one_letter not in AA1_TO_3.keys():
                raise KeyError(f"Invalid character: {one_letter}")
        return AA1_TO_3.get(one_letter.upper(), "Unk")

    def get_one_letter_code(self, three_letter: str, strict: bool = True) -> str:
        """Convert three-letter to one-letter amino acid code"""
        if strict:
            if three_letter not in AA3_TO_1.keys():
                raise KeyError(f"Invalid amino acid code: {three_letter}")
        return AA3_TO_1.get(three_letter.upper(), "X")
