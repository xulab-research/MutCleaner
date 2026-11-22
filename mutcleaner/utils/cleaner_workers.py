# mutcleaner/utils/cleaner_workers.py
from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd
from tqdm import tqdm

from ..core.alphabet import BaseAlphabet
from ..core.mutation import AminoAcidMutation, BaseMutation, CodonMutation, MutationSet
from ..core.sequence import ProteinSequence
from ..core.types import SequenceType
from .mutation_converter import invert_mutation_set

if TYPE_CHECKING:
    from typing import Any, Dict, List, Literal, Optional, Set, Tuple, Type, Union

    from pandas import Index

    from ..core.sequence import DNASequence, RNASequence

__all__ = [
    "valid_single_mutation",
    "apply_single_mutation",
    "infer_wt_sequence_grouped",
    "infer_single_mutationset",
    "validate_single_mutation_and_sequence",
]


def __dir__() -> List[str]:
    return __all__


def valid_single_mutation(args: Tuple) -> Tuple[Optional[str], Optional[str]]:
    """
    Process a single mutation string.

    Parameters
    ----------
    args : Tuple
        (mut_info, format_mutations, mutation_sep, offset, mutation_type, alphabet,
        include_ambiguous, include_sec_pyl, include_stop, reject_redundant, shared_cache).
        Share a cache only across calls with the same parsing and formatting options.

    Returns
    -------
    Tuple[Optional[str], Optional[str]]
        (formatted_mutation, error_message) - one will be None
    """
    mut_info, format_mutations, mutation_sep, offset, mutation_type, alphabet, include_ambiguous, include_sec_pyl, include_stop, reject_redundant, shared_cache = args

    if pd.isna(mut_info):
        return None, "Missing mutation information"

    try:
        # Check cache first
        if shared_cache is not None and mut_info in shared_cache:
            return shared_cache[mut_info]

        mutation_set = MutationSet.from_string(mut_info, sep=mutation_sep, offset=offset, mutation_type=mutation_type, alphabet=alphabet, include_ambiguous=include_ambiguous, include_sec_pyl=include_sec_pyl, include_stop=include_stop, reject_redundant=reject_redundant)
        result = (str(mutation_set) if format_mutations else mut_info, None)
        if shared_cache is not None:
            shared_cache[mut_info] = result
        return result

    except (ValueError, TypeError) as e:
        # MutationSet.from_string raises ValueError for invalid mutations
        # and TypeError for type-related issues
        error_msg = f"{type(e).__name__}: {str(e)}"
        if shared_cache is not None:
            shared_cache[mut_info] = (None, error_msg)
        return None, error_msg
    except Exception as e:
        error_msg = f"Unexpected error - {type(e).__name__}: {str(e)}"
        tqdm.write(f"Warning: Unexpected error processing mutation '{mut_info}': {e}")
        if shared_cache is not None:
            shared_cache[mut_info] = (None, error_msg)
        return None, error_msg


def apply_single_mutation(row_data: Tuple, dataset_columns: Index, sequence_column: str, name_column: str, mutation_column: str, position_columns: Optional[Dict[str, str]], mutation_sep: str, offset: int, mutation_type: Optional[Type[BaseMutation]], alphabet: Optional[BaseAlphabet], sequence_class: Type[Union[ProteinSequence, DNASequence, RNASequence]], *, include_ambiguous: bool = False, include_sec_pyl: bool = False, include_stop: bool = False, reject_redundant: bool = True) -> Tuple[Optional[str], Optional[str]]:
    """
    Apply mutations to a single sequence.

    Parameters
    ----------
    row_data : Tuple
        Row data from the dataset
    dataset_columns : Index
        Column names of the dataset
    sequence_column : str
        Column name containing sequences
    name_column : str
        Column name containing protein identifiers
    mutation_column : str
        Column name containing mutation information
    position_columns : Optional[Dict[str, str]]
        Position column mapping for sequence extraction
    mutation_sep : str
        Separator for splitting multiple mutations
    offset : int
        Input position of the first residue or codon.
    mutation_type : type of BaseMutation, optional
        Mutation class to use when parsing the annotation.
    alphabet : BaseAlphabet, optional
        Alphabet for parsing mutations. For proteins, it also validates the sequence.
    sequence_class : Type[Union[ProteinSequence, DNASequence, RNASequence]]
        Sequence class to use for mutation application
    include_ambiguous : bool, default=False
        Include B, Z, J, and X in the default protein alphabet.
    include_sec_pyl : bool, default=False
        Include U and O in the default protein alphabet.
    include_stop : bool, default=False
        Include the stop symbol in the default protein alphabet.
    reject_redundant : bool, default=True
        Reject an annotation containing an unchanged amino acid.

    Returns
    -------
    Tuple[Optional[str], Optional[str]]
        (mutated_sequence, error_message) - either sequence or error, not both
    """
    try:
        row_dict = dict(zip(dataset_columns, row_data))

        name = row_dict.get(name_column)
        mut_info = row_dict.get(mutation_column)
        sequence_str = row_dict.get(sequence_column)

        if not name or not mut_info:
            return None, f"Missing name or mutation info"

        if not sequence_str:
            return None, f"Missing sequence for {name}"

        # Apply position-based sequence extraction if specified
        if position_columns:
            start_pos = row_dict.get(position_columns.get("start", "start_pos"))
            end_pos = row_dict.get(position_columns.get("end", "end_pos"))
            if start_pos is not None and end_pos is not None:
                sequence_str = sequence_str[int(start_pos) : int(end_pos)]

        if sequence_class is ProteinSequence:
            sequence = ProteinSequence(sequence_str, alphabet=alphabet, name=name, include_ambiguous=include_ambiguous, include_sec_pyl=include_sec_pyl, include_stop=include_stop)
            alphabet = sequence.alphabet
        else:
            sequence = sequence_class(sequence_str, name=name)
        mutation_set = MutationSet.from_string(mut_info, sep=mutation_sep, offset=offset, mutation_type=mutation_type, alphabet=alphabet, reject_redundant=reject_redundant)
        mutated_sequence = sequence.apply_mutation(mutation_set)

        return str(mutated_sequence), None

    except Exception as e:
        return None, f"{type(e).__name__}: {str(e)}"


def _build_wt_row(
    measurement_row: Dict[str, Any],
    wt_sequence: str,
    mutation_column: str,
    sequence_column: str,
    label_columns: List[str],
    wt_label: float,
) -> Dict[str, Any]:
    """Build one wild-type row from the first measurement of its reference group."""
    wt_row = measurement_row.copy()
    wt_row[mutation_column] = "WT"
    wt_row[sequence_column] = wt_sequence
    for label_col in label_columns:
        if label_col in wt_row:
            wt_row[label_col] = wt_label
    return wt_row


def _rewrite_mutation_annotation(
    reference_sequence: str,
    mutated_sequence: str,
    mutation_type: Type[BaseMutation],
    alphabet: BaseAlphabet,
    mutation_sep: str,
    offset: int,
) -> str:
    """Rewrite a mutated sequence as a mutation annotation relative to a reference sequence."""
    if len(reference_sequence) != len(mutated_sequence):
        raise ValueError(f"sequence length differs from the reference: {len(mutated_sequence)} != {len(reference_sequence)}")

    if issubclass(mutation_type, AminoAcidMutation):
        mutations = [AminoAcidMutation(wild, position, mutant, alphabet) for position, (wild, mutant) in enumerate(zip(reference_sequence, mutated_sequence)) if wild != mutant]
        if not mutations:
            raise ValueError("mutated sequence equals the reference sequence")
        return mutation_sep.join(f"{mutation.wild_amino_acid}{mutation.position + offset}{mutation.mutant_amino_acid}" for mutation in mutations)

    if issubclass(mutation_type, CodonMutation):
        if len(reference_sequence) % 3 != 0:
            raise ValueError(f"sequence length {len(reference_sequence)} is not a whole number of codons")
        mutations = []
        for position in range(len(reference_sequence) // 3):
            wild_codon = reference_sequence[position * 3 : position * 3 + 3]
            mutant_codon = mutated_sequence[position * 3 : position * 3 + 3]
            if wild_codon != mutant_codon:
                mutations.append(CodonMutation(wild_codon, position, mutant_codon, alphabet))
        if not mutations:
            raise ValueError("mutated sequence equals the reference sequence")
        return mutation_sep.join(f"{mutation.wild_codon}{mutation.position + offset}{mutation.mutant_codon}" for mutation in mutations)

    raise TypeError(f"unsupported mutation type: {mutation_type.__name__}")


def infer_wt_sequence_grouped(
    group_data: Tuple[Any, pd.DataFrame],
    name_column: str,
    mutation_column: str,
    sequence_column: str,
    label_columns: List[str],
    wt_label: float,
    mutation_sep: str,
    offset: int,
    handle_multiple_wt: Literal["error", "separate", "first"],
    reserved_names: Set[str],
    sequence_class: Type[Union[ProteinSequence, DNASequence, RNASequence]],
    *,
    alphabet: Optional[BaseAlphabet],
    include_ambiguous: bool,
    include_sec_pyl: bool,
    include_stop: bool,
    reject_redundant: bool,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Infer wild-type sequences by reversing each row's substitutions.

    ``handle_multiple_wt`` decides what happens when a group yields several wild
    types: ``error`` fails the group, ``first`` keeps every measurement by
    rewriting annotations relative to the first inferred reference, and
    ``separate`` splits the group per reference and names the additional
    references ``{name}_variant{index}`` (skipping names in ``reserved_names``).

    Returns the rows kept for the dataset (plus appended wild-type rows) and the
    rows that could not be mapped to a reference, each carrying ``error_message``.
    """
    protein_name = group_data[0]
    group = group_data[1]

    try:
        measurements: List[Tuple[Dict[str, Any], str, Type[BaseMutation], Any]] = []
        for _, row in group.iterrows():
            mutated_sequence = str(row[sequence_column]).strip()
            if sequence_class is ProteinSequence:
                sequence = ProteinSequence(mutated_sequence, alphabet=alphabet, name=str(protein_name), include_ambiguous=include_ambiguous, include_sec_pyl=include_sec_pyl, include_stop=include_stop)
            else:
                sequence = sequence_class(mutated_sequence, alphabet=alphabet, name=str(protein_name))
            mutation_set = MutationSet.from_string(row[mutation_column], sep=mutation_sep, offset=offset, alphabet=sequence.alphabet, reject_redundant=reject_redundant)
            wt_sequence = str(sequence.apply_mutation(invert_mutation_set(mutation_set)))
            measurements.append((row.to_dict(), wt_sequence, mutation_set.mutation_type, sequence))

        reference_sequences = list(dict.fromkeys(wt_sequence for _, wt_sequence, _, _ in measurements))

        if len(reference_sequences) == 1:
            success_rows = [measurement_row for measurement_row, _, _, _ in measurements]
            success_rows.append(_build_wt_row(measurements[0][0], reference_sequences[0], mutation_column, sequence_column, label_columns, wt_label))
            return success_rows, []

        if handle_multiple_wt == "first":
            selected_reference = reference_sequences[0]
            success_rows = []
            failed_rows = []
            for measurement_row, wt_sequence, mutation_type, sequence in measurements:
                if wt_sequence == selected_reference:
                    success_rows.append(measurement_row)
                    continue

                try:
                    measurement_row[mutation_column] = _rewrite_mutation_annotation(selected_reference, str(sequence), mutation_type, sequence.alphabet, mutation_sep, offset)
                except (TypeError, ValueError) as error:
                    failed_row = measurement_row.copy()
                    failed_row["error_message"] = f"Cannot rewrite the measurement as a mutation of the first reference selected by handle_multiple_wt='first' for {protein_name}: {type(error).__name__}: {error}"
                    failed_rows.append(failed_row)
                    continue

                success_rows.append(measurement_row)

            success_rows.append(_build_wt_row(measurements[0][0], selected_reference, mutation_column, sequence_column, label_columns, wt_label))
            return success_rows, failed_rows

        if handle_multiple_wt == "separate":
            success_rows = []
            for index, reference_sequence in enumerate(reference_sequences):
                reference_rows = []
                for measurement_row, wt_sequence, _, _ in measurements:
                    if wt_sequence == reference_sequence:
                        reference_rows.append(measurement_row.copy())

                if index == 0:
                    reference_name = reference_rows[0][name_column]
                else:
                    variant_index = 1
                    reference_name = f"{protein_name}_variant{variant_index}"
                    while reference_name in reserved_names:
                        variant_index += 1
                        reference_name = f"{protein_name}_variant{variant_index}"
                    reserved_names.add(reference_name)

                for reference_row in reference_rows:
                    reference_row[name_column] = reference_name
                wt_row = _build_wt_row(reference_rows[0], reference_sequence, mutation_column, sequence_column, label_columns, wt_label)
                wt_row[name_column] = reference_name
                success_rows.extend(reference_rows)
                success_rows.append(wt_row)
            return success_rows, []

        error_row = group.iloc[0].to_dict()
        error_row["error_message"] = f"Multiple wildtype sequences inferred for {protein_name}: {len(reference_sequences)}"
        return [], [error_row]

    except Exception as e:
        error_row = group.iloc[0].to_dict() if len(group) > 0 else {name_column: str(protein_name)}
        error_row["error_message"] = f"{type(e).__name__}: {str(e)}"
        return [], [error_row]


def infer_single_mutationset(row_data: Tuple, dataset_columns: Index, wt_sequence_column: str, mut_sequence_column: str, mutation_sep: str, sequence_class: Type[SequenceType], *, alphabet: Optional[BaseAlphabet] = None, include_ambiguous: bool = False, include_sec_pyl: bool = False, include_stop: bool = False) -> Tuple[Optional[str], Optional[str]]:
    """
    Process a single row to infer mutations between WT and mutated sequences.

    This function is designed to be used with parallel processing.

    Parameters
    ----------
    row_data : Tuple
        Single row data containing sequence information
    dataset_columns : Index
        Column names of the dataset
    wt_sequence_column : str
        Column key for wild-type sequence
    mut_sequence_column : str
        Column key for mutated sequence
    mutation_sep : str
        Separator for joining multiple mutations
    sequence_class : Type[SequenceType]
        Sequence class to use for mutation inference
    alphabet : BaseAlphabet, optional
        Alphabet shared by the wild-type and mutated sequences.
    include_ambiguous : bool, default=False
        Allow ambiguous symbols in the default protein alphabet.
    include_sec_pyl : bool, default=False
        Allow U and O in the default protein alphabet.
    include_stop : bool, default=False
        Allow the stop symbol in the default protein alphabet.

    Returns
    -------
    tuple
        (result, error_message) where error_message is None on success
    """
    try:
        row_dict = dict(zip(dataset_columns, row_data))

        # Extract sequences
        wt_sequence = row_dict.get(wt_sequence_column, "")
        if sequence_class is ProteinSequence:
            wt_seq = ProteinSequence(wt_sequence, alphabet=alphabet, include_ambiguous=include_ambiguous, include_sec_pyl=include_sec_pyl, include_stop=include_stop)
        else:
            wt_seq = sequence_class(wt_sequence, alphabet=alphabet)
        mut_seq = sequence_class(row_dict.get(mut_sequence_column, ""), alphabet=wt_seq.alphabet)

        # Check sequence lengths are equal
        if len(wt_seq) != len(mut_seq):
            raise ValueError(f"Sequence length mismatch: WT={len(wt_seq)}, MUT={len(mut_seq)}")

        mutations = wt_seq.infer_mutation(mut_seq)
        inferred_mutations = mutation_sep.join(tuple(map(str, mutations.mutations)))
        return inferred_mutations, None
    except Exception as e:
        return None, str(e)


def validate_single_mutation_and_sequence(
    row_data: Tuple,
    dataset_columns: Index,
    wt_sequence_column: str,
    name_column: str,
    mutation_column: str,
    mut_sequence_column: str,
    mutation_sep: str,
    offset: int,
    sequence_class: Type[Union[ProteinSequence, DNASequence, RNASequence]],
) -> Tuple[Optional[str], Optional[str]]:

    try:
        row_dict = dict(zip(dataset_columns, row_data))

        name = row_dict.get(name_column)
        mut_info = row_dict.get(mutation_column)
        wt_sequence_str = row_dict.get(wt_sequence_column)
        mut_sequence_str = row_dict.get(mut_sequence_column)

        if not name or not mut_info:
            return None, f"Missing name or mutation info"

        if not wt_sequence_str or not mut_sequence_str:
            return None, f"Missing wildtype sequence or mutant sequence for {name}"

        wt_sequence = sequence_class(wt_sequence_str, name=name)
        mutation_set = MutationSet.from_string(mut_info, sep=mutation_sep, offset=offset)
        mut_sequence = sequence_class(mut_sequence_str, name=name)
        inferred_mutated_sequence = wt_sequence.apply_mutation(mutation_set)
        if str(inferred_mutated_sequence) != str(mut_sequence):
            return None, f"Mutated sequence does not match applied mutations for {name}"

        return str(inferred_mutated_sequence), None

    except Exception as e:
        return None, f"{type(e).__name__}: {str(e)}"
