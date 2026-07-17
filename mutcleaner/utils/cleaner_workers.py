# mutcleaner/utils/cleaner_workers.py
from __future__ import annotations

import pandas as pd

from tqdm import tqdm
from typing import TYPE_CHECKING

from .mutation_converter import invert_mutation_set
from ..core.mutation import BaseMutation, MutationSet
from ..core.alphabet import BaseAlphabet
from ..core.types import SequenceType

if TYPE_CHECKING:
    from pandas import Index
    from typing import Any, Dict, List, Literal, Optional, Tuple, Type, Union

    from ..core.alphabet import ProteinAlphabet, DNAAlphabet, RNAAlphabet
    from ..core.sequence import ProteinSequence, DNASequence, RNASequence

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
        (mut_info, format_mutations, mutation_sep, is_zero_based, cache)

    Returns
    -------
    Tuple[Optional[str], Optional[str]]
        (formatted_mutation, error_message) - one will be None
    """
    mut_info, format_mutations, mutation_sep, is_zero_based, mutation_type, alphabet, shared_cache = args

    if pd.isna(mut_info):
        return None, "Missing mutation information"

    try:
        # Check cache first
        if shared_cache is not None and mut_info in shared_cache:
            return shared_cache[mut_info]

        if format_mutations:
            # Parse and format mutation
            mutation_set = MutationSet.from_string(
                mut_info,
                sep=mutation_sep,
                is_zero_based=is_zero_based,
                mutation_type=mutation_type,
                alphabet=alphabet,
            )
            formatted_mut = str(mutation_set)

            if shared_cache is not None:
                shared_cache[mut_info] = (formatted_mut, None)
            return formatted_mut, None
        else:
            # Only validate, don't format
            # Try to create MutationSet to validate - if it succeeds, mutation is valid
            MutationSet.from_string(
                mut_info, sep=mutation_sep, is_zero_based=is_zero_based, mutation_type=mutation_type, alphabet=alphabet
            )
            # If no exception was raised, the mutation is valid
            if shared_cache is not None:
                shared_cache[mut_info] = (mut_info, None)
            return mut_info, None

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


def apply_single_mutation(
    row_data: Tuple,
    dataset_columns: Index,
    sequence_column: str,
    name_column: str,
    mutation_column: str,
    position_columns: Optional[Dict[str, str]],
    mutation_sep: str,
    is_zero_based: bool,
    mutation_type: Optional[Type[BaseMutation]],
    alphabet: Optional[BaseAlphabet],
    sequence_class: Type[Union[ProteinSequence, DNASequence, RNASequence]],
) -> Tuple[Optional[str], Optional[str]]:
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
    is_zero_based : bool
        Whether the mutation position is zero-based.
    sequence_class : Type[Union[ProteinSequence, DNASequence, RNASequence]]
        Sequence class to use for mutation application

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

        sequence = sequence_class(sequence_str, name=name)
        mutation_set = MutationSet.from_string(
            mut_info, sep=mutation_sep, is_zero_based=is_zero_based, mutation_type=mutation_type, alphabet=alphabet
        )
        mutated_sequence = sequence.apply_mutation(mutation_set)

        return str(mutated_sequence), None

    except Exception as e:
        return None, f"{type(e).__name__}: {str(e)}"


def infer_wt_sequence_grouped(
    group_data: Tuple[Any, pd.DataFrame],
    name_column: str,
    mutation_column: str,
    sequence_column: str,
    label_columns: List[str],
    wt_label: float,
    mutation_sep: str,
    is_zero_based: bool,
    handle_multiple_wt: Literal["error", "separate", "first"],
    sequence_class: Type[Union[ProteinSequence, DNASequence, RNASequence]],
    alphabet_class: Type[Union[ProteinAlphabet, DNAAlphabet, RNAAlphabet]],
) -> Tuple[List[Dict[str, Any]], str]:
    """
    Process a single protein group and return list of rows (including WT).

    This is a module-level function that processes protein groups independently.
    """
    protein_name = group_data[0]
    group = group_data[1]

    try:
        # Infer wild-type sequences
        inferred_wt_seqs = set()
        result_rows = []

        # First, add all original rows to result
        for _, row in group.iterrows():
            result_rows.append(row.to_dict())

        # Then, infer WT sequences
        for _, row in group.iterrows():
            mut_info = row[mutation_column]
            mut_seq = row[sequence_column]

            # Parse mutation and create sequence
            mutation_set = MutationSet.from_string(
                mut_info,
                sep=mutation_sep,
                is_zero_based=is_zero_based,
                alphabet=alphabet_class(),
            )
            sequence = sequence_class(str(mut_seq).strip(), name=str(protein_name))

            # Infer wild-type sequence by applying inverted mutations
            inverted_mutation_set = invert_mutation_set(mutation_set)
            wt_seq = sequence.apply_mutation(inverted_mutation_set)
            inferred_wt_seqs.add(str(wt_seq))

        # Handle wild-type sequences
        if len(inferred_wt_seqs) == 1:
            # Single wild-type sequence - add one WT row
            wt_seq_str = inferred_wt_seqs.pop()

            # Create WT row based on first row of the group
            first_row = group.iloc[0].to_dict()
            wt_row = first_row.copy()

            # Update WT-specific fields
            wt_row[mutation_column] = "WT"  # or empty string if preferred
            wt_row[sequence_column] = wt_seq_str

            # Set labels to `wt_label` for WT
            for label_col in label_columns:
                if label_col in wt_row:
                    wt_row[label_col] = wt_label

            result_rows.append(wt_row)
            return result_rows, "success"

        elif len(inferred_wt_seqs) > 1:
            # Multiple wild-type sequences
            if handle_multiple_wt == "first":
                # Take the first WT sequence
                wt_seq_str = list(inferred_wt_seqs)[0]

                first_row = group.iloc[0].to_dict()
                wt_row = first_row.copy()
                wt_row[mutation_column] = "WT"
                wt_row[sequence_column] = wt_seq_str

                for label_col in label_columns:
                    if label_col in wt_row:
                        wt_row[label_col] = wt_label

                result_rows.append(wt_row)
                return result_rows, "success"

            elif handle_multiple_wt == "separate":
                # Add multiple WT rows
                for i, wt_seq_str in enumerate(inferred_wt_seqs):
                    first_row = group.iloc[0].to_dict()
                    wt_row = first_row.copy()

                    if i == 0:
                        wt_row[mutation_column] = "WT"
                    else:
                        wt_row[mutation_column] = f"WT_variant{i}"

                    wt_row[sequence_column] = wt_seq_str

                    for label_col in label_columns:
                        if label_col in wt_row:
                            wt_row[label_col] = wt_label

                    result_rows.append(wt_row)

                return result_rows, "failed"

            else:  # handle_multiple_wt == "error"
                # Add error information to the first row
                error_row = group.iloc[0].to_dict()
                error_row["error_message"] = (
                    f"Multiple wildtype sequences inferred for {protein_name}: {len(inferred_wt_seqs)}"
                )
                return [error_row], "failed"

        else:
            # No wild-type sequences inferred
            error_row = group.iloc[0].to_dict()
            error_row["error_message"] = (
                f"No wildtype sequences could be inferred for {protein_name}"
            )
            return [error_row], "failed"

    except Exception as e:
        # Save error information in first row
        error_row = (
            group.iloc[0].to_dict()
            if len(group) > 0
            else {name_column: str(protein_name)}
        )
        error_row["error_message"] = f"{type(e).__name__}: {str(e)}"
        return [error_row], "failed"


def infer_single_mutationset(
    row_data: Tuple,
    dataset_columns: Index,
    wt_sequence_column: str,
    mut_sequence_column: str,
    mutation_sep: str,
    sequence_class: Type[SequenceType],
) -> Tuple[Optional[str], Optional[str]]:
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

    Returns
    -------
    tuple
        (result, error_message) where error_message is None on success
    """
    try:
        row_dict = dict(zip(dataset_columns, row_data))

        # Extract sequences
        wt_seq = sequence_class(row_dict.get(wt_sequence_column, ""))
        mut_seq = sequence_class(row_dict.get(mut_sequence_column, ""))

        # Check sequence lengths are equal
        if len(wt_seq) != len(mut_seq):
            raise ValueError(
                f"Sequence length mismatch: WT={len(wt_seq)}, MUT={len(mut_seq)}"
            )

        mutations = wt_seq.infer_mutation(mut_seq)  # type: ignore[arg-type]
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
    is_zero_based: bool,
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
        mutation_set = MutationSet.from_string(
            mut_info, sep=mutation_sep, is_zero_based=is_zero_based
        )
        mut_sequence = sequence_class(mut_sequence_str, name=name)
        inferred_mutated_sequence = wt_sequence.apply_mutation(mutation_set)
        if str(inferred_mutated_sequence) != str(mut_sequence):
            return None, f"Mutated sequence does not match applied mutations for {name}"

        return str(inferred_mutated_sequence), None

    except Exception as e:
        return None, f"{type(e).__name__}: {str(e)}"