# mutcleaner/core/mutation.py
from __future__ import annotations

import re
from abc import ABC, abstractmethod
from typing import Generic, TYPE_CHECKING

from .alphabet import ProteinAlphabet, DNAAlphabet, RNAAlphabet
from .codon import CodonTable
from .types import MutationType

if TYPE_CHECKING:
    from typing import Any, Dict, List, Literal, Optional, Sequence, Set, Type

    from .alphabet import BaseAlphabet

__all__ = [
    "AminoAcidMutation",
    "AminoAcidMutationSet",
    "BaseMutation",
    "CodonMutation",
    "CodonMutationSet",
    "MutationSet",
]


def __dir__() -> List[str]:
    return __all__


class BaseMutation(ABC):
    """Base class for all mutations"""

    def __init__(
        self,
        wild_type: str,
        mutant_type: str,
        position: int,
        alphabet: Optional[BaseAlphabet] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ):
        if position < 0:
            raise ValueError("Position must be non-negative (0-indexed)")
        self.wild_type = wild_type
        self.mutant_type = mutant_type
        self.position = position
        self.alphabet = alphabet
        self.metadata = metadata or {}

    @property
    @abstractmethod
    def type(self) -> str:
        """Get the type of mutation"""
        pass

    @abstractmethod
    def __str__(self) -> str:
        """Return string representation of the mutation"""
        pass

    @abstractmethod
    def _is_valid(self) -> bool:
        """Check if the mutation is valid"""
        pass

    @abstractmethod
    def get_mutation_category(self) -> str:
        """Get mutation category"""
        pass

    @classmethod
    @abstractmethod
    def from_string(
        cls,
        mutation_string: str,
        is_zero_based: bool,
        alphabet: Optional[BaseAlphabet] = None,
    ) -> "BaseMutation":
        """Parse mutation from string format like 'A123V' or 'Ala123Val'"""
        pass

    def __eq__(self, other) -> bool:
        """Check if two mutations are equal"""
        if not isinstance(other, self.__class__):
            return False
        return self.position == other.position and str(self) == str(other) and self.type == other.type

    def __hash__(self) -> int:
        """Enable use in sets and as dict keys"""
        return hash((self.__class__, self.position, str(self), self.type))


class AminoAcidMutation(BaseMutation):
    """Represents an amino acid mutation (e.g., A123V)"""

    def __init__(
        self,
        wild_type: str,
        position: int,
        mutant_type: str,
        alphabet: Optional[ProteinAlphabet] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ):
        super().__init__(wild_type, mutant_type, position, alphabet, metadata)

        self.alphabet = alphabet or ProteinAlphabet(include_stop=True)

        self.wild_amino_acid = wild_type.upper()
        self.mutant_amino_acid = mutant_type.upper()

        if not self._is_valid():
            raise ValueError(f"Invalid amino acid mutation: {self}")

    @property
    def type(self) -> str:
        """Get the type of mutation"""
        return "amino_acid"

    @property
    def effect_type(self) -> Literal["synonymous", "nonsense", "missense"]:
        """Get the effect type of the mutation (synonymous, nonsense, or missense)"""
        return self.get_mutation_category()

    def __str__(self) -> str:
        return f"{self.wild_amino_acid}{self.position}{self.mutant_amino_acid}"

    def _is_valid(self) -> bool:
        """Check if mutation uses valid amino acid codes"""
        return self.wild_amino_acid in self.alphabet and self.mutant_amino_acid in self.alphabet and isinstance(self.position, int) and self.position >= 0

    def is_synonymous(self) -> bool:
        """Check if mutation is synonymous (no change)"""
        return self.wild_amino_acid == self.mutant_amino_acid

    def is_nonsense(self) -> bool:
        """Check if mutation introduces stop codon"""
        return self.mutant_amino_acid == "*"

    def is_missense(self) -> bool:
        """Check if mutation is missense (changes amino acid)"""
        return not self.is_synonymous() and not self.is_nonsense()

    def get_mutation_category(self) -> Literal["synonymous", "nonsense", "missense"]:
        """Get mutation classification"""
        if self.is_synonymous():
            return "synonymous"
        elif self.is_nonsense():
            return "nonsense"
        else:
            return "missense"

    @classmethod
    def from_string(
        cls,
        mutation_str: str,
        is_zero_based: bool = False,
        alphabet: Optional[ProteinAlphabet] = None,
    ) -> "AminoAcidMutation":
        """Parse mutation from string format like 'A123V' or 'Ala123Val'"""
        mutation_str = mutation_str.strip()

        # Handle three-letter codes first
        three_letter_pattern = r"^([A-Za-z]{3})(\d+)([A-Za-z]{3})$"
        match = re.match(three_letter_pattern, mutation_str)

        if match:
            alphabet = alphabet or ProteinAlphabet(include_stop=True)
            wild_three, position, mutant_three = match.groups()

            try:
                wild_amino_acid = alphabet.get_one_letter_code(wild_three)
                mutant_amino_acid = alphabet.get_one_letter_code(mutant_three)
            except KeyError as e:
                raise ValueError(f"Unknown three-letter amino acid code: {e}")
            if is_zero_based:
                return cls(wild_amino_acid, int(position), mutant_amino_acid, alphabet)
            else:
                return cls(wild_amino_acid, int(position) - 1, mutant_amino_acid, alphabet)

        # Handle one-letter codes
        one_letter_pattern = r"^([A-Z\*])(\d+)([A-Z\*])$"
        match = re.match(one_letter_pattern, mutation_str)

        if not match:
            raise ValueError(f"Invalid mutation format: {mutation_str}. " f"Expected formats: 'A123V' or 'Ala123Val'")

        wild_amino_acid, position, mutant_amino_acid = match.groups()
        if is_zero_based:
            return cls(wild_amino_acid, int(position), mutant_amino_acid, alphabet)
        else:
            return cls(wild_amino_acid, int(position) - 1, mutant_amino_acid, alphabet)


class CodonMutation(BaseMutation):
    """Represent a codon substitution.

    The ``position`` attribute is always stored as a zero-based codon
    index relative to the coding sequence. The corresponding nucleotide
    interval is ``position * 3 : position * 3 + 3``.
    """

    def __init__(
        self,
        wild_type: str,
        position: int,
        mutant_type: str,
        alphabet: Optional[BaseAlphabet] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ):
        super().__init__(wild_type, mutant_type, position, alphabet, metadata)

        self.wild_codon = wild_type.upper()
        self.mutant_codon = mutant_type.upper()

        # Auto-detect sequence type based on presence of T/U
        self.seq_type: Literal["DNA", "RNA", "Both"] = self._detect_seq_type(self.wild_codon, self.mutant_codon)

        # Use appropriate alphabet based on detected sequence type
        self.alphabet = alphabet if alphabet is not None else (RNAAlphabet() if self.seq_type == "RNA" else DNAAlphabet())
        # print("Alphabet:", self.alphabet, "seq_type:", self.seq_type)

        if not self._is_valid():
            raise ValueError(f"Invalid codon mutation: {self}")

    @property
    def nucleotide_start(self) -> int:
        """Return the zero-based nucleotide start of this codon."""
        return self.position * 3

    @property
    def nucleotide_end(self) -> int:
        """Return the exclusive nucleotide end of this codon."""
        return self.nucleotide_start + 3

    @staticmethod
    def _detect_seq_type(wild_codon: str, mutant_codon: str) -> Literal["DNA", "RNA", "Both"]:
        """Auto-detect sequence type based on T/U presence"""
        combined_sequence = (wild_codon + mutant_codon).upper()

        has_t = "T" in combined_sequence
        has_u = "U" in combined_sequence

        if has_t and has_u:
            raise ValueError("Codons cannot contain both T and U")
        elif has_t:
            return "DNA"
        elif has_u:
            return "RNA"
        else:
            return "Both"

    @property
    def type(self) -> str:
        """Get the type of mutation"""
        return f"codon_{self.seq_type.lower()}"

    def __str__(self) -> str:
        return f"{self.wild_codon}{self.position}{self.mutant_codon}"

    def _is_valid(self) -> bool:
        """Check if codons are valid"""
        return len(self.wild_codon) == 3 and len(self.mutant_codon) == 3 and self.alphabet.is_valid_sequence(self.wild_codon) and self.alphabet.is_valid_sequence(self.mutant_codon) and isinstance(self.position, int) and self.position >= 0

    def get_mutation_category(self) -> str:
        return f"codon_{self.seq_type.lower()}"

    def to_amino_acid_mutation(self, codon_table: Optional[CodonTable] = None) -> AminoAcidMutation:
        """Convert codon mutation to amino acid mutation"""
        if codon_table is None:
            # Use appropriate codon table based on detected sequence type
            if self.seq_type == "Both":
                # Default to DNA when ambiguous
                codon_table = CodonTable.get_standard_table("DNA")
            else:
                codon_table = CodonTable.get_standard_table(self.seq_type)

        wild_aa = codon_table.translate_codon(self.wild_codon)
        mutant_aa = codon_table.translate_codon(self.mutant_codon)

        return AminoAcidMutation(wild_aa, self.position, mutant_aa, metadata=self.metadata.copy())

    @classmethod
    def from_string(
        cls,
        mutation_str: str,
        is_zero_based: bool = False,
        alphabet: Optional[BaseAlphabet] = None,
    ) -> "CodonMutation":
        """Parse mutation from string format like 'ATG123TAA' or 'AUG123UAA'"""
        mutation_str = mutation_str.strip()

        # Handle codon format first
        codon_pattern = r"^([ATUCG]{3})(\d+)([ATUCG]{3})$"
        match = re.match(codon_pattern, mutation_str)

        if match:
            wild_codon, position, mutant_codon = match.groups()
            alphabet = alphabet if alphabet is not None else (RNAAlphabet() if cls._detect_seq_type(wild_codon, mutant_codon) == "RNA" else DNAAlphabet())
            if is_zero_based:
                return cls(wild_codon, int(position), mutant_codon, alphabet)
            else:
                return cls(wild_codon, int(position) - 1, mutant_codon, alphabet)

        raise ValueError(f"Invalid codon mutation format: {mutation_str}")


class MutationSet(Generic[MutationType]):
    """Represents a set of mutations of the same type"""

    def __init__(
        self,
        mutations: Sequence[MutationType],
        mutation_type: Optional[Type[MutationType]],
        name: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ):
        if not mutations:
            raise ValueError("MutationSet must contain at least one mutation")

        if mutation_type is None:
            # guess the mutation type based on the first mutation
            mutation_type = type(mutations[0])

        # Validate that all mutations are of the same type and mutation type
        self._validate_mutation_types(mutations, mutation_type)

        # Validate that all mutations have the same type property
        self._validate_mutation_type_consistency(mutations)

        # Check for duplicate positions
        self._validate_unique_positions(mutations)

        self.mutations: List = sorted(list(mutations), key=lambda m: m.position)
        self.mutation_type = mutation_type
        self.name = name
        self.metadata = metadata or {}

    def __str__(self) -> str:
        return ",".join(str(m) for m in self.mutations)

    def __repr__(self) -> str:
        return f"MutationSet(mutations={self.mutations!r}, mutation_type={self.mutation_type}, name={self.name!r})"

    def __len__(self) -> int:
        """Return number of mutations"""
        return len(self.mutations)

    def __iter__(self):
        """Make the mutation set iterable"""
        return iter(self.mutations)

    @property
    def mutation_subtype(self) -> str:
        """Get the specific mutation subtype (e.g., 'amino_acid', 'codon_dna', 'codon_rna', 'codon_both')"""
        if self.mutations:
            return self.mutations[0].type
        return "unknown"

    def _validate_mutation_types(self, mutations: Sequence[MutationType], expected_type: Type[MutationType]) -> None:
        """Validate that all mutations are of the expected class type"""
        if not issubclass(expected_type, BaseMutation):
            raise TypeError(f"Expected mutations to be of type {expected_type.__name__}")

        invalid_mutations = [m for m in mutations if not isinstance(m, expected_type)]
        if invalid_mutations:
            raise ValueError(f"All mutations must be of type {expected_type.__name__}. " f"Found {len(invalid_mutations)} mutations of different types.")

    def _validate_mutation_type_consistency(self, mutations: Sequence[MutationType]) -> None:
        """Validate that all mutations have the same type property"""
        if not mutations:
            return

        seq_types = {
            mutation.seq_type
            for mutation in mutations
            if isinstance(mutation, CodonMutation)
            and mutation.seq_type != "Both"
        }

        if seq_types <= {"DNA"}:
            return

        if seq_types <= {"RNA"}:
            return
    
        if not seq_types:
            return

        raise ValueError(
            "Codon mutations contain incompatible DNA and RNA types"
        )

    def _validate_unique_positions(self, mutations: Sequence[MutationType]) -> None:
        """Validate that mutations have unique positions"""
        positions = [m.position for m in mutations]
        if len(positions) != len(set(positions)):
            duplicates = [pos for pos in positions if positions.count(pos) > 1]
            raise ValueError(f"Duplicate mutations at positions: {set(duplicates)}")

    def add_mutation(self, mutation: MutationType) -> None:
        """Add a mutation to this set"""
        # Validate mutation class type
        if not isinstance(mutation, self.mutation_type):
            raise ValueError(f"Mutation must be of type {self.mutation_type.__name__}, " f"got {type(mutation).__name__}")

        # Validate mutation type property consistency
        if self.mutations and mutation.type != self.mutation_subtype:
            raise ValueError(f"Mutation type property must match existing mutations. " f"Expected '{self.mutation_subtype}', got '{mutation.type}'")

        # Check for position conflict
        if mutation.position in self.get_positions():
            raise ValueError(f"Mutation already exists at position {mutation.position}")

        self.mutations.append(mutation)
        self.mutations.sort(key=lambda m: m.position)

    def remove_mutation(self, position: int) -> bool:
        """Remove mutation at specified position, return True if removed"""
        original_length = len(self.mutations)
        self.mutations = [m for m in self.mutations if m.position != position]
        return len(self.mutations) < original_length

    def get_mutation_at(self, position: int) -> Optional[MutationType]:
        """Get mutation at specified position"""
        for mutation in self.mutations:
            if mutation.position == position:
                return mutation
        return None

    def has_mutation_at(self, position: int) -> bool:
        """Check if there is a mutation at specified position"""
        return position in self.get_positions()

    def is_single_mutation(self) -> bool:
        """Check if this is a single mutation"""
        return len(self.mutations) == 1

    def is_multiple_mutations(self) -> bool:
        """Check if this contains multiple mutations"""
        return len(self.mutations) > 1

    def get_mutation_count(self) -> int:
        """Get number of mutations"""
        return len(self.mutations)

    def validate_all(self) -> bool:
        """Validate all mutations"""
        return all(mutation._is_valid() for mutation in self.mutations)

    def get_positions(self) -> List[int]:
        """Get all mutation positions"""
        return [mutation.position for mutation in self.mutations]

    def get_positions_set(self) -> Set[int]:
        """Get all mutation positions as a set"""
        return set(self.get_positions())

    def get_mutation_categories(self) -> Dict[str, int]:
        """Get mutation category statistics"""
        categories = {}
        for mutation in self.mutations:
            category = mutation.get_mutation_category()
            categories[category] = categories.get(category, 0) + 1
        return categories

    def filter_by_category(self, category: str) -> List[MutationType]:
        """Filter mutations by category"""
        return [m for m in self.mutations if m.get_mutation_category() == category]

    def sort_by_position(self) -> None:
        """Sort mutations by position in ascending order"""
        self.mutations.sort(key=lambda m: m.position)

    def get_sorted_by_position(self) -> List[MutationType]:
        """Get mutations sorted by position without modifying the original list"""
        return sorted(self.mutations, key=lambda m: m.position)

    @classmethod
    def from_string(
        cls,
        string: str,
        sep: Optional[str] = None,
        is_zero_based: bool = False,
        reject_redundant: bool = True,
        reject_nonsense: bool = True,
        mutation_type: Optional[Type[MutationType]] = None,
        alphabet: Optional[BaseAlphabet] = None,
        name: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> "MutationSet":
        """
        Create a mutation set from a string

        Parameters
        ----------
        string : str
            String containing mutations separated by delimiter
        sep : Optional[str], default=None
            Separator to use. If None, will try to guess
        is_zero_based : bool, default=False
            Whether origin mutation positions are zero-based
        reject_redundant : bool, default=True
            Whether to reject synonymous/no-op mutations
        reject_nonsense : bool, default=True
            Whether to reject nonsense mutations.
        mutation_type : Type[MutationType], default=None
            The type of mutations to create. If None, will infer from first mutation
        alphabet : Optional[BaseAlphabet], default=None
            Alphabet to use for mutation parsing (if applicable)
        name : str, default=None
            Optional name for the mutation set
        metadata : Optional[Dict[str, Any]], default=None
            Optional metadata for the mutation set

        Returns
        -------
        MutationSet
            A MutationSet created from the input string.
            Returns AminoAcidMutationSet for amino acid mutations,
            CodonMutationSet for codon mutations, or generic MutationSet for others.

        Raises
        ------
        ValueError
            If string is empty, no valid mutations found, or mutations are inconsistent
        TypeError
            If mutation_type is not a subclass of BaseMutation
        """
        if not string or not string.strip():
            raise ValueError("Input string cannot be empty")

        string = string.strip()

        # Guess separator if not provided
        if sep is None:
            sep = cls._guess_sep(string)
            if sep is None:
                # Assume single mutation if no separator found
                mutation_parts = [string]
            else:
                mutation_parts = string.split(sep)
        else:
            mutation_parts = string.split(sep)

        # Clean up mutation strings
        mutation_parts = [part.strip() for part in mutation_parts if part.strip()]

        if not mutation_parts:
            raise ValueError("No valid mutation strings found after splitting")

        mutations = []
        errors = []

        for i, mutation_str in enumerate(mutation_parts):
            try:
                if mutation_type is None:
                    # Try to infer mutation type from the first mutation string
                    mutation = cls._infer_and_create_mutation(mutation_str, is_zero_based, alphabet)
                    mutation_type = type(mutation)
                else:
                    # Use specified mutation type
                    if hasattr(mutation_type, "from_string"):
                        mutation = mutation_type.from_string(mutation_str, is_zero_based, alphabet)
                    else:
                        raise NotImplementedError(f"Mutation type {mutation_type.__name__} does not have from_string method")

                mutations.append(mutation)

            except Exception as e:
                errors.append(f"Error parsing mutation '{mutation_str}' at position {i}: {str(e)}")

        if not mutations:
            error_msg = "No valid mutations could be parsed"
            if errors:
                error_msg += f". Errors encountered: {'; '.join(errors)}"
            raise ValueError(error_msg)

        if errors:
            # Raise an exception if we have some valid mutations but also errors
            raise ValueError(f"Some mutations could not be parsed: {'; '.join(errors)}")

        if reject_redundant:
            redundant = [str(m) for m in mutations if hasattr(m, "is_synonymous") and callable(getattr(m, "is_synonymous")) and m.is_synonymous()]
            if redundant:
                raise ValueError(f"Redundant/no-op mutations detected:{', '.join(redundant)}")
        if reject_nonsense:
            nonsense = [str(m) for m in mutations if hasattr(m, "is_nonsense") and callable(getattr(m, "is_nonsense")) and m.is_nonsense()]
            if nonsense:
                raise ValueError(f"Nonsense/stop-gain mutations detected: {', '.join(nonsense)}")
        # Return appropriate mutation set type based on detected mutation type
        if mutation_type == AminoAcidMutation:
            return AminoAcidMutationSet(
                mutations=mutations,
                name=name,
                metadata=metadata,
            )
        elif mutation_type == CodonMutation:
            return CodonMutationSet(
                mutations=mutations,
                name=name,
                metadata=metadata,
            )
        else:
            # For other mutation types or when called on subclasses, use the generic approach
            return cls(
                mutations=mutations,
                mutation_type=mutation_type,
                name=name,
                metadata=metadata,
            )

    @classmethod
    def _infer_and_create_mutation(
        cls,
        mutation_str: str,
        is_zero_based: bool = False,
        alphabet: Optional[BaseAlphabet] = None,
    ) -> MutationType:
        """
        Infer mutation type and create mutation from string

        This method tries different mutation types to see which one can parse the string.

        Parameters
        ----------
        mutation_str : str
            Mutation string to parse
        is_zero_based : bool, default=False
            Whether origin mutation positions are zero-based
        alphabet : Optional[BaseAlphabet]
            Optional alphabet to use for parsing

        Returns
        -------
        MutationType
            The inferred mutation type that can parse the string
        """
        # List of mutation types to try (order matters - most common first)
        mutation_types = [
            AminoAcidMutation,
            CodonMutation,
        ]

        last_error = None

        for mutation_type in mutation_types:
            try:
                return mutation_type.from_string(mutation_str, is_zero_based, alphabet)
            except Exception as e:
                last_error = e
                continue

        # none of the mutation types could parse the string
        raise ValueError(f"Could not parse mutation string '{mutation_str}' with any known mutation type. " f"Last error: {last_error}")

    @classmethod
    def _create_mutation(
        cls,
        mutation_str: str,
        mutation_type: Type[BaseMutation],
        is_zero_based: bool = False,
        alphabet: Optional[BaseAlphabet] = None,
    ) -> BaseMutation:
        return mutation_type.from_string(mutation_str, is_zero_based, alphabet)

    @staticmethod
    def _guess_sep(string: str) -> Optional[str]:
        """Guess the separator for a string of mutations"""
        if not string:
            return None

        candidate_separators = [";", ",", "|", ":", "/", "\\", "\t"]

        # Count occurrences of each separator
        separator_counts = {sep: string.count(sep) for sep in candidate_separators}

        # Filter out separators that don't appear
        valid_separators = {sep: count for sep, count in separator_counts.items() if count > 0}

        if not valid_separators:
            return None

        # Choose the separator with the highest count
        # In case of tie, prefer the order in candidate_separators
        best_sep = None
        best_count = 0

        for sep in candidate_separators:
            if sep in valid_separators and valid_separators[sep] > best_count:
                best_sep = sep
                best_count = valid_separators[sep]

        return best_sep


class AminoAcidMutationSet(MutationSet[AminoAcidMutation]):
    """Represents a set of amino acid mutations"""

    def __init__(
        self,
        mutations: Sequence[AminoAcidMutation],
        name: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ):
        super().__init__(mutations, AminoAcidMutation, name, metadata)

    def get_synonymous_mutations(self) -> List[AminoAcidMutation]:
        """Get all synonymous mutations"""
        return [m for m in self.mutations if m.is_synonymous()]

    def get_missense_mutations(self) -> List[AminoAcidMutation]:
        """Get all missense mutations"""
        return [m for m in self.mutations if m.is_missense()]

    def get_nonsense_mutations(self) -> List[AminoAcidMutation]:
        """Get all nonsense mutations"""
        return [m for m in self.mutations if m.is_nonsense()]

    def has_stop_codon_mutations(self) -> bool:
        """Check if any mutations introduce stop codons"""
        return any(m.is_nonsense() for m in self.mutations)

    def count_by_effect_type(self) -> Dict[str, int]:
        """Count mutations by effect type"""
        return {
            "synonymous": len(self.get_synonymous_mutations()),
            "missense": len(self.get_missense_mutations()),
            "nonsense": len(self.get_nonsense_mutations()),
        }


class CodonMutationSet(MutationSet[CodonMutation]):
    """Represents a set of codon mutations"""

    def __init__(
        self,
        mutations: Sequence[CodonMutation],
        name: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ):
        super().__init__(mutations, CodonMutation, name, metadata)

    @property
    def seq_type(self) -> Literal["DNA", "RNA", "Both"]:
        """Get the sequence type (DNA, RNA, or Both) of the codon mutations"""
        if self.mutations:
            return self.mutations[0].seq_type
        return "Both"

    def to_amino_acid_mutation_set(self, codon_table: Optional[CodonTable] = None) -> AminoAcidMutationSet:
        """Convert all codon mutations to amino acid mutations"""
        if codon_table is None:
            # Use appropriate codon table based on detected sequence type
            if self.seq_type == "Both":
                # Default to DNA when ambiguous
                codon_table = CodonTable.get_standard_table("DNA")
            else:
                codon_table = CodonTable.get_standard_table(self.seq_type)

        aa_mutations = [mutation.to_amino_acid_mutation(codon_table) for mutation in self.mutations]

        return AminoAcidMutationSet(
            aa_mutations,
            name=f"{self.name}_aa" if self.name else None,
            metadata=self.metadata.copy(),
        )
