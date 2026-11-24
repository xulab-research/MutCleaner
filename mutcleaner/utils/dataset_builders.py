# mutcleaner/utils/dataset_builders.py
from __future__ import annotations
from ..core.mutation import (
    AminoAcidMutation,
    BaseMutation,
    CodonMutation,
    MutationSet,
)
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any, Dict, List, Optional, Type, Tuple, Union

    from ..core.alphabet import BaseAlphabet
    from ..core.sequence import (
        BaseSequence,
        DNASequence,
        ProteinSequence,
        RNASequence,
    )
import pandas as pd
from tqdm import tqdm

"""
Functions are used in mutcleaner.cleaners.basic_cleaners.convert_to_mutation_dataset_format()

format 1:

>>> pd.DataFrame({
...     'name': ['prot1', 'prot1', 'prot1', 'prot2', 'prot2'],
...     'mut_info': ['A0S,Q1D', 'C2D', 'WT', 'E0F', 'WT'],
...     'mut_seq': ['SDCDEF', 'AQDDEF', 'AQCDEF', 'FGHIGHK', 'EGHIGHK'],
...     'score': [1.5, 2.0, 0.0, 3.0, 0.0]
... })

format 2:

>>> df2 = pd.DataFrame({
...     'name': ['prot1', 'prot1', 'prot2'],
...     'sequence': ['AKCDEF', 'AKCDEF', 'FEGHIS'],
...     'mut_info': ['A0K,C2D', 'Q1P', 'E1F'],
...     'score': [1.5, 2.0, 3.0],
...     'mut_seq': ['KKDDEF', 'APCDEF', 'FFGHIS']
... })
"""

__all__ = ["convert_format_1", "convert_format_2"]


def __dir__() -> List[str]:
    return __all__


def convert_format_1(
    df: pd.DataFrame,
    name_column: str,
    mutation_column: str,
    mutated_sequence_column: str,
    score_column: str,
    include_wild_type: bool,
    mutation_set_prefix: str,
    offset: int,
    additional_metadata: Optional[Dict[str, Any]],
    sequence_class: Type[Union[ProteinSequence, DNASequence, RNASequence]],
    mutation_type: Type[BaseMutation],
    alphabet: BaseAlphabet,
) -> Tuple[pd.DataFrame, Dict[str, BaseSequence]]:
    """Convert Format 1 (with WT rows) to mutation dataset format."""

    input_df = df.copy()

    # Extract reference sequences from WT rows
    wt_rows = input_df[input_df[mutation_column] == "WT"]
    if wt_rows.empty:
        raise ValueError("No wild-type (WT) entries found in the dataset")

    reference_sequences = {}
    for _, row in wt_rows.iterrows():
        name = row[name_column]
        sequence = row[mutated_sequence_column]
        reference_sequences[name] = sequence_class(sequence, alphabet=alphabet)

    # Filter out wild-type entries if requested
    if not include_wild_type:
        input_df = input_df[input_df[mutation_column] != "WT"].copy()

    if input_df.empty:
        raise ValueError("No mutation data remaining after filtering")

    # Process mutations (now supporting multi-mutations)
    output_rows = []
    total_rows = len(input_df)
    for idx, (mut_info, name, score) in tqdm(enumerate(zip(input_df[mutation_column], input_df[name_column], input_df[score_column])), total=total_rows):
        # Skip wild-type if it somehow made it through filtering
        if mut_info == "WT":
            continue

        # Parse mutations (single or multiple)
        try:
            mutations = _parse_mutations_string(
                mut_info,
                offset=offset,
                mutation_type=mutation_type,
                alphabet=alphabet,
            )
        except ValueError as error:
            raise ValueError(f"Cannot parse mutation {mut_info!r} in row {idx}: {error}") from error
        # Create one output row per individual mutation within the set
        mutation_set_id = f"{mutation_set_prefix}_{idx + 1}"
        mutation_set_name = f"{name}_{mut_info}"

        for mutation in mutations:
            output_row = _create_output_row_from_mutation(
                mutation_set_id,
                mutation_set_name,
                mut_info,
                name,
                score,
                mutation,
                additional_metadata,
            )
            output_rows.append(output_row)

    _check_mutation_offset(output_rows, reference_sequences, offset)
    output_df = pd.DataFrame(output_rows)
    return output_df, reference_sequences


def convert_format_2(
    df: pd.DataFrame,
    name_column: str,
    mutation_column: str,
    sequence_column: str,
    score_column: str,
    mutation_set_prefix: str,
    offset: int,
    additional_metadata: Optional[Dict[str, Any]],
    sequence_class: Type[Union[ProteinSequence, DNASequence, RNASequence]],
    mutation_type: Type[BaseMutation],
    alphabet: BaseAlphabet,
) -> Tuple[pd.DataFrame, Dict[str, BaseSequence]]:
    """Convert Format 2 (with sequence column) to mutation dataset format."""

    input_df = df.copy()

    # Extract reference sequences from sequence column
    reference_sequences = {}
    for name, group in tqdm(input_df.groupby(name_column)):
        sequences = group[sequence_column].unique()
        if len(sequences) > 1:
            raise ValueError(f"Multiple different sequences found for protein '{name}': {sequences}")
        reference_sequences[name] = sequence_class(sequences[0], alphabet=alphabet)

    # Process mutations (now supporting multi-mutations)
    output_rows = []
    total_rows = len(input_df)
    for idx, (mut_info, name, score) in tqdm(enumerate(zip(input_df[mutation_column], input_df[name_column], input_df[score_column])), total=total_rows):
        # Parse mutations (single or multiple)
        try:
            mutations = _parse_mutations_string(
                mut_info,
                offset=offset,
                mutation_type=mutation_type,
                alphabet=alphabet,
            )
        except ValueError as e:
            raise ValueError(f"Cannot parse mutation '{mut_info}' in row {idx}: {e}")

        # Create one output row per individual mutation within the set
        mutation_set_id = f"{mutation_set_prefix}_{idx + 1}"
        mutation_set_name = f"{name}_{mut_info}"

        for mutation in mutations:
            output_row = _create_output_row_from_mutation(
                mutation_set_id,
                mutation_set_name,
                mut_info,
                name,
                score,
                mutation,
                additional_metadata,
            )
            output_rows.append(output_row)

    _check_mutation_offset(output_rows, reference_sequences, offset)
    output_df = pd.DataFrame(output_rows)
    return output_df, reference_sequences


def _check_mutation_offset(output_rows: List[Dict[str, Any]], reference_sequences: Dict[str, BaseSequence], offset: int) -> None:
    """Flag a one-position shift when it yields a perfect reference match."""
    candidate_offsets = (offset - 1, offset, offset + 1)
    reference_text = {name: str(sequence) for name, sequence in reference_sequences.items()}
    matches_by_reference = {name: dict.fromkeys(candidate_offsets, 0) for name in reference_sequences}
    totals_by_reference = dict.fromkeys(reference_sequences, 0)

    for row in output_rows:
        name = row["reference_id"]
        if name not in reference_text:
            continue
        sequence = reference_text[name]
        wild_type = row["wild_codon"] if "wild_codon" in row else row["wild_amino_acid"]
        width = len(wild_type)
        totals_by_reference[name] += 1
        for candidate in candidate_offsets:
            start = (row["position"] + offset - candidate) * width
            matches_by_reference[name][candidate] += start >= 0 and sequence[start : start + width] == wild_type

    for name, matches in matches_by_reference.items():
        total = totals_by_reference[name]
        if matches[offset] < total and total in matches.values():
            evidence = ", ".join(f"offset {candidate}: {matched}/{total}" for candidate, matched in matches.items())
            raise ValueError(f"Possible mutation numbering mismatch for reference {name!r}. Wild-type matches by offset: {evidence}. Verify the input numbering before converting. If the annotations were already formatted, check the offset used to format them.")


def _create_output_row_from_mutation(
    mutation_set_id: str,
    mutation_set_name: str,
    original_mutation_string: str,
    name: str,
    score: float,
    mutation: BaseMutation,
    additional_metadata: Optional[Dict[str, Any]],
) -> Dict[str, Any]:
    """Create one flattened MutationDataset row from a mutation object."""
    output_row = {
        "mutation_set_id": mutation_set_id,
        "reference_id": name,
        "mutation_string": str(mutation),
        "position": mutation.position,
        "mutation_type": mutation.type,
        "mutation_set_name": mutation_set_name,
        "label": score,
        "set_original_mutation_string": original_mutation_string,
    }

    if isinstance(mutation, AminoAcidMutation):
        output_row.update(
            {
                "wild_amino_acid": mutation.wild_amino_acid,
                "mutant_amino_acid": mutation.mutant_amino_acid,
            }
        )

    elif isinstance(mutation, CodonMutation):
        output_row.update(
            {
                "wild_codon": mutation.wild_codon,
                "mutant_codon": mutation.mutant_codon,
                "sequence_type": mutation.seq_type,
                "position_unit": "codon",
            }
        )

    else:
        raise TypeError(f"Unsupported mutation type: {type(mutation).__name__}")

    if additional_metadata:
        output_row.update({f"set_{key}": value for key, value in additional_metadata.items()})

    return output_row


def _parse_mutations_string(
    mutation_string: str,
    offset: int,
    mutation_type: Type[BaseMutation],
    alphabet: BaseAlphabet,
) -> List[BaseMutation]:
    """Parse one or more mutations from a mutation string.

    Parameters
    ----------
    mutation_string : str
        String containing one or more mutations.
    offset : int
        Input position of the first residue or codon.
    mutation_type : Type[BaseMutation]
        Mutation class used for parsing.
    alphabet : BaseAlphabet
        Alphabet used to validate mutation symbols.

    Returns
    -------
    List[BaseMutation]
        Parsed mutation objects.
    """
    if not isinstance(mutation_string, str) or not mutation_string.strip():
        raise ValueError("Mutation string cannot be empty")

    mutation_set = MutationSet.from_string(
        mutation_string.strip(),
        offset=offset,
        mutation_type=mutation_type,
        alphabet=alphabet,
    )

    return list(mutation_set.mutations)
