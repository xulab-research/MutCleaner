# mutcleaner/core/sequence.py
from __future__ import annotations

import warnings
from abc import ABC
from pathlib import Path
from typing import TYPE_CHECKING

from .alphabet import DNAAlphabet, RNAAlphabet, ProteinAlphabet
from .codon import CodonTable
from .constants import DNA_BASE_COMPLEMENTS, RNA_BASE_COMPLEMENTS
from .mutation import (
    BaseMutation,
    CodonMutation,
    AminoAcidMutation,
    MutationSet,
    CodonMutationSet,
    AminoAcidMutationSet,
)

if TYPE_CHECKING:
    from typing import Callable, Dict, List, Literal, Optional, Type, Union

    from .alphabet import BaseAlphabet
    from .types import SequenceType


__all__ = ["BaseSequence", "DNASequence", "ProteinSequence", "RNASequence"]


def __dir__() -> List[str]:
    return __all__


class BaseSequence(ABC):
    """Base class for biological sequences"""

    def __init__(
        self,
        sequence: str,
        alphabet: Optional[BaseAlphabet] = None,
        name: Optional[str] = None,
        metadata: Optional[Dict] = None,
    ):
        if not sequence:
            raise ValueError("Sequence cannot be empty")
        if alphabet is None:
            # If not provided, request the subclass's default
            alphabet = type(self).default_alphabet()
        if alphabet is None:
            # No default provided by the subclass → raise in the base class
            raise TypeError(
                f"{type(self).__name__} requires 'alphabet' (no default provided)"
            )
        self.alphabet = alphabet
        self.sequence = self.alphabet.validate_sequence(sequence)
        self.name = name
        self.metadata = metadata or {}

    def __len__(self) -> int:
        return len(self.sequence)

    def __str__(self) -> str:
        return self.sequence

    def __getitem__(self: SequenceType, key) -> SequenceType:
        if isinstance(key, slice):
            return type(self)(
                self.sequence[key],
                self.alphabet,
                (
                    f"{self.name}_{key.start}_{key.stop if key.stop is not None else len(self.sequence)}"
                    if self.name
                    else None
                ),
                self.metadata.copy(),
            )
        elif isinstance(key, int):
            return type(self)(
                self.sequence[key],
                self.alphabet,
                f"{self.name}_{key}_{key+1}" if self.name else None,
                self.metadata.copy(),
            )
        else:
            raise TypeError(f"Invalid argument type: {type(key).__name__}")

    def __eq__(self, other) -> bool:
        if isinstance(other, str):
            return self.sequence == other
        if not isinstance(other, BaseSequence):
            return False
        return self.sequence == other.sequence and type(self.alphabet) == type(
            other.alphabet
        )

    @classmethod
    def default_alphabet(cls) -> Optional[BaseAlphabet]:
        """Subclasses may override this method to provide a default alphabet.

        By default, it returns None, indicating no default is provided
        and callers must pass `alphabet` explicitly."""
        raise NotImplementedError("Subclasses must implement this method.")

    def get_subsequence(
        self: SequenceType, start: int, end: Optional[int] = None
    ) -> SequenceType:
        """get subsequence (0-indexed, inclusive)"""
        if start < 0:
            raise IndexError("Start position must be greater than or equal to 0")
        if end is not None and end < start:
            raise ValueError("End position must be greater than or equal to start.")

        subseq = self[start:end]
        return subseq

    def apply_mutation(
        self: SequenceType,
        mutation: Union[BaseMutation, MutationSet, MutationSet[BaseMutation]],
    ) -> SequenceType:
        """
        Apply a mutation or set of mutations to the sequence and return a new sequence.

        Parameters
        ----------
        mutation : Union[BaseMutation, CodonMutationSet, AminoAcidMutationSet]
            A single mutation or a set of mutations to apply

        Returns
        -------
        SequenceType
            A new sequence with the mutation(s) applied

        Raises
        ------
        ValueError
            If mutation position is invalid or mutation is incompatible
        TypeError
            If mutation type is not supported
        """
        # Handle mutation sets (multiple mutations)
        if isinstance(mutation, (AminoAcidMutationSet, CodonMutationSet)):
            # Apply mutations in reverse order of position to avoid index shifting
            mutations = sorted(
                mutation.mutations, key=lambda m: m.position, reverse=True
            )
            result_sequence = self

            for single_mutation in mutations:
                result_sequence = result_sequence.apply_mutation(single_mutation)

            return result_sequence

        # Handle single mutations
        elif isinstance(mutation, BaseMutation):
            # Validate mutation position and get mutation details
            if isinstance(mutation, CodonMutation):
                codon_position = mutation.position
                start_pos = mutation.nucleotide_start
                end_pos = mutation.nucleotide_end

                # For codon mutations, check 3 bases starting at position
                if start_pos < 0 or end_pos > len(self.sequence):
                    raise ValueError(
                        f"Codon position {codon_position} is out of bounds "
                        f"for a sequence of length {len(self.sequence)} "
                        f"({len(self.sequence) // 3} complete codons)"
                    )

                # Check mutation subtypes (DNA or RNA)
                valid_combinations = {
                    "DNA": DNAAlphabet,
                    "Both": DNAAlphabet,
                    "RNA": RNAAlphabet,
                }
                expected_alphabet = valid_combinations.get(mutation.seq_type)

                if expected_alphabet is None or not isinstance(
                    self.alphabet, expected_alphabet
                ):
                    raise TypeError(
                        f"Unmatching mutation subtype: {mutation.seq_type} with {mutation.seq_type} sequence"
                    )

                # Validate original codon matches expected
                actual_codon = self.sequence[start_pos:end_pos]
                if actual_codon != mutation.wild_codon:
                    raise ValueError(
                        f"Expected codon '{mutation.wild_codon}' at position {mutation.position}, "
                        f"but found '{actual_codon}'"
                    )

                # Apply codon mutation (replace 3 bases)
                new_sequence = (
                    self.sequence[:start_pos]
                    + mutation.mutant_codon
                    + self.sequence[end_pos:]
                )

            elif isinstance(mutation, AminoAcidMutation):
                # For amino acid mutations, check single position
                if mutation.position < 0 or mutation.position >= len(self.sequence):
                    raise ValueError(
                        f"Amino acid mutation position {mutation.position} is out of bounds for sequence of length {len(self.sequence)}"
                    )

                # Validate original amino acid matches expected
                actual_aa = self.sequence[mutation.position]
                if actual_aa != mutation.wild_amino_acid:
                    raise ValueError(
                        f"Expected amino acid '{mutation.wild_amino_acid}' at position {mutation.position}, "
                        f"but found '{actual_aa}'"
                    )

                # Apply amino acid mutation (replace single position)
                new_sequence = (
                    self.sequence[: mutation.position]
                    + mutation.mutant_amino_acid
                    + self.sequence[mutation.position + 1 :]
                )

            else:
                # Handle other BaseMutation subclasses generically
                if mutation.position < 0 or mutation.position >= len(self.sequence):
                    raise ValueError(
                        f"Mutation position {mutation.position} is out of bounds for sequence of length {len(self.sequence)}"
                    )

                # For generic mutations, we can't validate original or determine replacement length
                # This is a fallback for custom mutation types
                raise TypeError(
                    f"Unsupported mutation subtype: {type(mutation).__name__}. "
                    f"Only CodonMutation and AminoAcidMutation are supported."
                )

            # Update metadata to track mutation
            new_metadata = self.metadata.copy()
            if "mutations_applied" not in new_metadata:
                new_metadata["mutations_applied"] = []

            mutation_record = {
                "type": type(mutation).__name__,
                "mutation_type": mutation.type,
                "position": mutation.position,
            }

            # Add type-specific information
            if isinstance(mutation, CodonMutation):
                mutation_record.update(
                    {
                        "wild_codon": mutation.wild_codon,
                        "mutant_codon": mutation.mutant_codon,
                        "seq_type": mutation.seq_type,
                    }
                )
            elif isinstance(mutation, AminoAcidMutation):
                mutation_record.update(
                    {
                        "wild_amino_acid": mutation.wild_amino_acid,
                        "mutant_amino_acid": mutation.mutant_amino_acid,
                        "effect_type": mutation.effect_type,
                    }
                )

            new_metadata["mutations_applied"].append(mutation_record)

            # Create new sequence instance
            return type(self)(new_sequence, self.alphabet, self.name, new_metadata)

        else:
            raise TypeError(f"Unsupported mutation type: {type(mutation).__name__}")

    def infer_mutation(self: SequenceType, other: SequenceType) -> MutationSet:
        """Infer a mutation that leads to a specific sequence"""
        if type(self) is not type(other):
            raise TypeError(
                "Sequences must be of the same type, "
                f"got {type(self).__name__} and {type(other).__name__}"
            )

        if len(self.sequence) != len(other.sequence):
            raise ValueError(
                "Sequences must have the same length, "
                f"got {len(self.sequence)} and {len(other.sequence)}"
            )

        if isinstance(self, ProteinSequence):
            mutations = [
                AminoAcidMutation(
                    wild_type=wild_amino_acid,
                    position=position,
                    mutant_type=mutant_amino_acid,
                    alphabet=self.alphabet,
                )
                for position, (wild_amino_acid, mutant_amino_acid) in enumerate(
                    zip(self.sequence, other.sequence)
                )
                if wild_amino_acid != mutant_amino_acid
            ]

            if not mutations:
                raise ValueError("No amino-acid mutations were inferred")

            return AminoAcidMutationSet(mutations)

        if isinstance(self, (DNASequence, RNASequence)):
            if len(self.sequence) % 3 != 0:
                raise ValueError(
                    f"{type(self).__name__} mutation inference requires sequence "
                    f"lengths divisible by 3, got {len(self.sequence)}"
                )

            mutations = []

            for codon_position in range(len(self.sequence) // 3):
                start = codon_position * 3
                end = start + 3

                wild_codon = self.sequence[start:end]
                mutant_codon = other.sequence[start:end]

                if wild_codon != mutant_codon:
                    mutations.append(
                        CodonMutation(
                            wild_type=wild_codon,
                            position=codon_position,
                            mutant_type=mutant_codon,
                            alphabet=self.alphabet,
                        )
                    )

            if not mutations:
                raise ValueError("No codon mutations were inferred")

            return CodonMutationSet(mutations)

        raise TypeError(
            f"Mutation inference is not supported for "
            f"{type(self).__name__}"
        )


class ProteinSequence(BaseSequence):
    """Protein sequence with amino acid validation"""

    def __init__(
        self,
        sequence: str,
        alphabet: Optional[ProteinAlphabet] = None,
        name: Optional[str] = None,
        metadata: Optional[Dict] = None,
    ):
        super().__init__(sequence, alphabet, name, metadata)

    @classmethod
    def default_alphabet(cls) -> Optional[BaseAlphabet]:
        return ProteinAlphabet(include_stop=True)

    def get_residue(self, position: int) -> str:
        """Get amino acid at specific position (0-indexed)"""
        if position < 0 or position >= len(self.sequence):
            raise IndexError(
                f"Position {position} out of range (0-{len(self.sequence)})"
            )
        return self.sequence[position]

    def find_motif(self, motif: str) -> List[int]:
        """Find all positions where motif occurs (0-indexed)"""
        positions = []
        motif = motif.upper()
        start = 0
        while True:
            pos = self.sequence.find(motif, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1
        return positions


class RNASequence(BaseSequence):
    """RNA sequence with nucleotide validation"""

    def __init__(
        self,
        sequence: str,
        alphabet: Optional[RNAAlphabet] = None,
        name: Optional[str] = None,
        metadata: Optional[Dict] = None,
    ):
        super().__init__(sequence, alphabet, name, metadata)

    @classmethod
    def default_alphabet(cls) -> Optional[BaseAlphabet]:
        return RNAAlphabet()

    def reverse_complement(self) -> "RNASequence":
        """Get reverse complement of RNA sequence"""
        try:
            rev_comp = "".join(
                RNA_BASE_COMPLEMENTS[base] for base in self.sequence[::-1]
            )
        except KeyError as e:
            raise ValueError(f"Invalid RNA base found: {e}")
        return RNASequence(
            sequence=rev_comp,
            name=f"{self.name}_rc" if self.name else None,
            metadata=self.metadata,
        )

    def back_transcribe(self) -> "DNASequence":
        """Back-transcribe RNA sequence into DNA sequence"""
        dna_seq = self.sequence.replace("U", "T")
        return DNASequence(
            sequence=dna_seq,
            name=f"{self.name}_back_transcribe" if self.name else None,
            metadata=self.metadata,
        )

    def translate(
        self,
        codon_table: Optional[CodonTable] = None,
        start_at_first_met: bool = False,
        stop_at_stop_codon: bool = False,
        require_mod3: bool = True,
        start: Optional[int] = None,
        end: Optional[int] = None,
    ) -> ProteinSequence:
        """
        Translate RNA sequence into amino acid sequence using this codon table.

        Parameters
        ----------
        codon_table : Optional[CodonTable], default=None
            Codon table to use for translation. If None, uses standard genetic code.
        start_at_first_met : bool, default=False
            Start translation at the first start codon if found.
        stop_at_stop_codon : bool, default=False
            Stop translation when a stop codon is encountered.
        require_mod3 : bool, default=True
            Whether the sequence must be a multiple of 3 in length.
        start : Option[int], default=None
            Custom 0-based start position. Overrides `start_at_first_met`.
        end : Option[int], default=None
            Custom 0-based end position. Overrides `stop_at_stop_codon`.

        Returns
        -------
        ProteinSequence
            Translated amino acid sequence.
        """
        aa_seq = translate(
            sequence=self.sequence,
            seq_type="RNA",
            codon_table=codon_table,
            start_at_first_met=start_at_first_met,
            stop_at_stop_codon=stop_at_stop_codon,
            require_mod3=require_mod3,
            start=start,
            end=end,
        )
        return ProteinSequence(
            sequence=aa_seq,
            name=f"{self.name}_translation" if self.name else None,
            metadata=self.metadata,
        )


class DNASequence(BaseSequence):
    """DNA sequence with nucleotide validation"""

    def __init__(
        self,
        sequence: str,
        alphabet: Optional[DNAAlphabet] = None,
        name: Optional[str] = None,
        metadata: Optional[Dict] = None,
    ):
        super().__init__(sequence, alphabet, name, metadata)

    @classmethod
    def default_alphabet(cls) -> Optional[BaseAlphabet]:
        return DNAAlphabet()

    def reverse_complement(self) -> "DNASequence":
        """Get reverse complement of DNA sequence"""
        try:
            rev_comp = "".join(
                DNA_BASE_COMPLEMENTS[base] for base in self.sequence[::-1]
            )
        except KeyError as e:
            raise ValueError(f"Invalid DNA base found: {e}")
        return DNASequence(
            sequence=rev_comp,
            name=f"{self.name}_rc" if self.name else None,
            metadata=self.metadata,
        )

    def translate(
        self,
        codon_table: Optional[CodonTable] = None,
        start_at_first_met: bool = False,
        stop_at_stop_codon: bool = False,
        require_mod3: bool = True,
        start: Optional[int] = None,
        end: Optional[int] = None,
    ) -> ProteinSequence:
        """
        Translate DNA sequence into amino acid sequence using this codon table.

        Parameters
        ----------
        codon_table : Optional[CodonTable], default=None
            Codon table to use for translation. If None, uses standard genetic code.
        start_at_first_met : bool, default=False
            Start translation at the first start codon if found.
        stop_at_stop_codon : bool, default=False
            Stop translation when a stop codon is encountered.
        require_mod3 : bool, default=True
            Whether the sequence must be a multiple of 3 in length.
        start : Optional[int], default=None
            Custom 0-based start position. Overrides `start_at_first_met`.
        end : Optional[int], default=None
            Custom 0-based end position. Overrides `stop_at_stop_codon`.

        Returns
        -------
        ProteinSequence
            Translated amino acid sequence.
        """
        aa_seq = translate(
            sequence=self.sequence,
            seq_type="DNA",
            codon_table=codon_table,
            start_at_first_met=start_at_first_met,
            stop_at_stop_codon=stop_at_stop_codon,
            require_mod3=require_mod3,
            start=start,
            end=end,
        )
        return ProteinSequence(
            sequence=aa_seq,
            name=f"{self.name}_translation" if self.name else None,
            metadata=self.metadata,
        )

    def transcribe(self) -> "RNASequence":
        """Transcribe DNA sequence into RNA sequence"""
        rna_seq = self.sequence.replace("T", "U")
        return RNASequence(
            sequence=rna_seq,
            name=f"{self.name}_transcribed" if self.name else None,
            metadata=self.metadata,
        )


def translate(
    sequence: str,
    seq_type: Literal["DNA", "RNA"] = "DNA",
    codon_table: Optional[CodonTable] = None,
    start_at_first_met: bool = False,
    stop_at_stop_codon: bool = False,
    require_mod3: bool = True,
    start: Optional[int] = None,
    end: Optional[int] = None,
) -> str:
    """
    Translate DNA or RNA sequence into amino acid sequence using this codon table.
    this function should not be called directly. Instead, use the DNASequence or RNASequence classes.

    Parameters
    ----------
    sequence : str
        DNA or RNA sequence to translate.
    seq_type : Literal["DNA", "RNA"], default="DNA"
    codon_table : Optional[CodonTable], default=None
        Codon table to use for translation. If None, uses standard genetic code.
    start_at_first_met : bool, default=False
        Start translation at the first start codon if found.
    stop_at_stop_codon : bool, default=False
        Stop translation when a stop codon is encountered.
    require_mod3 : bool, default=True
        Whether the sequence must be a multiple of 3 in length.
    start : Option[int], default=None
        Custom 0-based start position. Overrides `start_at_first_met`.
    end : Option[int], default=None
        Custom 0-based end position. Overrides `stop_at_stop_codon`.

    Returns
    -------
    str
        Translated amino acid sequence.
    """
    n = len(sequence)

    # Use standard table if none provided
    if codon_table is None:
        codon_table = CodonTable.get_standard_table(seq_type=seq_type)

    # Auto detect start if `start` not provided
    if start is None:
        if start_at_first_met:
            for i in range(0, n - 2, 3):
                codon = sequence[i : i + 3]
                if codon_table.is_start_codon(codon):
                    start = i
                    break
            else:
                return ""  # No start codon found
        else:
            start = 0

    # Auto detect end if `end` not provided
    if end is None:
        if stop_at_stop_codon:
            for i in range(start, n - 2, 3):
                codon = sequence[i : i + 3]
                if codon_table.is_stop_codon(codon):
                    end = i + 3
                    break
            else:
                end = n
        else:
            end = n

    sub_seq = sequence[start:end]
    if len(sub_seq) % 3 != 0:
        remainder = len(sub_seq) % 3
        if require_mod3:
            raise ValueError(
                f"Sequence length from start={start} to end={end} is not divisible by 3 "
                f"(remainder = {remainder})."
            )
        else:
            warnings.warn(
                f"Sequence length from start={start} to end={end} is not divisible by 3. "
                f"Discarding {remainder} trailing nucleotide(s): {sub_seq[-remainder:]}"
            )
            sub_seq = sub_seq[: len(sub_seq) - remainder]

    # Translate using this codon table
    codons = [sub_seq[i : i + 3] for i in range(0, len(sub_seq), 3)]
    aa_seq = "".join([codon_table.translate_codon(codon) for codon in codons])

    return aa_seq


def load_sequences_from_fasta(
    fasta_path: Union[str, Path],
    sequence_class: Type[BaseSequence],
    alphabet: Optional[BaseAlphabet] = None,
    header_func: Optional[Callable[[str], tuple[str, str]]] = None,
    allow_duplicates: bool = False,
) -> Dict[str, BaseSequence]:
    """
    Load sequences from a FASTA file into a dictionary of BaseSequence-derived objects.

    Parameters
    ----------
    fasta_path : Union[str, Path]
        Path to the FASTA file
    sequence_class : Type[BaseSequence]
        A subclass of BaseSequence to instantiate (e.g. DNASequence)
    alphabet : Optional[BaseAlphabet], default=None
        Optional alphabet to validate the sequence. If None, uses default for sequence class.
    header_func : Optional[Callable[[str], tuple[str, str]]], default=None
        Function to process header line and extract (sequence_id, description).
        If None, uses default parsing (first word as ID, rest as description).
    allow_duplicates : bool, default=False
        If False, raises error on duplicate sequence IDs. If True, overwrites.

    Returns
    -------
    Dict[str, BaseSequence]
        Dictionary of {sequence_id: sequence_object}

    Raises
    ------
    FileNotFoundError
        If FASTA file doesn't exist
    ValueError
        If duplicate sequence IDs found and allow_duplicates=False
    TypeError
        If sequence_class is not a subclass of BaseSequence
    """
    # Validate inputs
    if not issubclass(sequence_class, BaseSequence):
        raise TypeError(
            f"sequence_class must be a subclass of BaseSequence, got {sequence_class}"
        )

    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    # Default header processing function
    def default_header_func(header_line: str) -> tuple[str, str]:
        """Extract ID and description from header line (without '>')"""
        parts = header_line.split(maxsplit=1)
        seq_id = parts[0]
        description = parts[1] if len(parts) > 1 else ""
        return seq_id, description

    # Use provided header function or default
    process_header = header_func or default_header_func

    # Get default alphabet for the sequence class if not provided
    def get_default_alphabet():
        """Get default alphabet for the sequence class"""
        if sequence_class == DNASequence:
            return DNAAlphabet()
        elif sequence_class == RNASequence:
            return RNAAlphabet()
        elif sequence_class == ProteinSequence:
            return ProteinAlphabet(include_stop=True)
        else:
            # TODO: For custom sequence classes, try to instantiate with None
            # and let the class handle default alphabet
            # but for easy, raise an error
            raise TypeError(
                f"sequence_class {sequence_class} does not have a default alphabet, "
                f"please provide it with the alphabet parameter."
            )

    sequences: Dict[str, BaseSequence] = {}
    current_id = None
    current_description = ""
    current_seq_lines = []
    line_number = 0

    try:
        with fasta_path.open("r", encoding="utf-8") as f:
            for line in f:
                line_number += 1
                line = line.strip()

                # Skip empty lines
                if not line:
                    continue

                if line.startswith(">"):
                    # Save previous entry if exists
                    if current_id is not None:
                        _save_sequence(
                            sequences,
                            current_id,
                            current_description,
                            current_seq_lines,
                            sequence_class,
                            alphabet or get_default_alphabet(),
                            allow_duplicates,
                            line_number,
                        )

                    # Process new header
                    try:
                        current_id, current_description = process_header(line[1:])
                    except Exception as e:
                        raise ValueError(
                            f"Error processing header at line {line_number}: '{line}'. {str(e)}"
                        )

                    current_seq_lines = []

                else:
                    # Accumulate sequence lines
                    # Remove any whitespace and validate characters
                    clean_line = "".join(line.split())  # Remove all whitespace
                    if clean_line:  # Only add non-empty lines
                        current_seq_lines.append(clean_line)

            # Save last entry
            if current_id is not None:
                _save_sequence(
                    sequences,
                    current_id,
                    current_description,
                    current_seq_lines,
                    sequence_class,
                    alphabet or get_default_alphabet(),
                    allow_duplicates,
                    line_number,
                )

    except UnicodeDecodeError:
        raise ValueError(
            f"Unable to decode file {fasta_path}. Please ensure it's a valid text file."
        )
    except Exception as e:
        raise ValueError(f"Error reading FASTA file at line {line_number}: {str(e)}")

    if not sequences:
        raise ValueError(f"No valid sequences found in {fasta_path}")

    return sequences


def _save_sequence(
    sequences: Dict[str, BaseSequence],
    seq_id: str,
    description: str,
    seq_lines: List[str],
    sequence_class: Type[BaseSequence],
    alphabet: BaseAlphabet,
    allow_duplicates: bool,
    line_number: int,
) -> None:
    """Helper function to save a sequence to the dictionary"""
    # Check for duplicate IDs
    if seq_id in sequences and not allow_duplicates:
        raise ValueError(
            f"Duplicate sequence ID '{seq_id}' found near line {line_number}"
        )

    # Join sequence lines and convert to uppercase
    full_seq = "".join(seq_lines).upper()

    # Skip empty sequences
    if not full_seq:
        warnings.warn(f"Empty sequence found for ID '{seq_id}', skipping")
        return

    try:
        # Create sequence object
        sequences[seq_id] = sequence_class(
            sequence=full_seq,
            alphabet=alphabet,
            name=seq_id,
            metadata={"description": description, "source": "fasta"},
        )
    except Exception as e:
        raise ValueError(f"Error creating sequence object for ID '{seq_id}': {str(e)}")
