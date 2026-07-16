# cleaners/cdna_proteolysis_custom_cleaners.py
from __future__ import annotations

import pandas as pd
from functools import partial
from joblib import Parallel, delayed
from tqdm import tqdm
from typing import TYPE_CHECKING

from ..core.mutation import MutationSet
from ..core.pipeline import multiout_step
from ..core.sequence import ProteinSequence
from ..utils.mutation_converter import invert_mutation_set

if TYPE_CHECKING:
    from typing import Any, Dict, List, Sequence, Tuple, Union


__all__ = [
    "validate_wt_sequence",
]


def __dir__() -> List[str]:
    return __all__


@multiout_step(main="success", failed="failed")
def validate_wt_sequence(
    dataset: pd.DataFrame,
    name_column: str,
    mutation_column: str,
    sequence_column: str,
    wt_identifier: str = "wt",
    num_workers: int = 4,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Validate wild-type (WT) sequences per protein/name group by checking
    consistency between the *explicit* WT row (if present) and the WT sequence
    *inferred* from each mutant row. If the explicit WT row is missing, a new
    one is synthesized using the inferred sequence.

    This step groups the dataset by ``name_column`` and, for each group:
      1) Reconstructs the WT sequence by inverting mutation sets from all mutant rows.
      2) If an explicit WT row is present (i.e., ``mutation_column == wt_identifier``),
         ensures its sequence matches the consensus inferred WT sequence.
      3) If the explicit WT row is missing, synthesizes a new WT row with the
         inferred sequence and null experimental values.
      4) Fails the group if (a) multiple distinct inferred WT sequences are found,
         or (b) the inferred WT sequence does not match the explicit WT sequence.
         Otherwise, the group is marked as success.

    Parallel execution is performed with ``joblib.Parallel(backend="loky")`` and
    a progress bar via ``tqdm``. If parallel execution fails, a sequential
    fallback is used. Diagnostic warnings are printed with ``tqdm.write``.

    Parameters
    ----------
    dataset : pd.DataFrame
        Input dataframe containing at least the name, mutation, and sequence columns.
    name_column : str
        Column that identifies a protein/name group for validation.
    mutation_column : str
        Column containing mutation annotations. Its format must be parsable by
        ``MutationSet.from_string(sep=",", is_zero_based=True)``.
    sequence_column : str
        Column containing the amino-acid sequence for that row (treated as the
        mutated sequence for mutants and the WT sequence for the WT row).
    wt_identifier : str, default "wt"
        The exact token used in ``mutation_column`` to mark the explicit WT row.
    num_workers : int, default 4
        Number of worker processes for parallel execution.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        - successful_df:
            Concatenation of original and synthesized rows for validated groups.
            Both mutant rows and the explicit (or synthesized) WT rows from successful
            groups are returned. May be empty if no group passes.
        - failed_df:
            Rows summarizing groups that failed validation. Each failed group is
            represented by at least one row with an additional column ``"error_message"``
            describing the reason, e.g., multiple inferred WT sequences, or mismatch
            between inferred and explicit WT sequences. May be empty.

    Examples
    --------
    >>> # Minimal schema example showing both explicit and synthesized WT cases:
    >>> df = pd.DataFrame({
    ...     "protein": ["P1", "P1", "P1", "P2", "P2"],
    ...     "mut":     ["wt", "A0G", "L4F", "K3R", "M0A"],
    ...     "seq":     ["ACDELG", "GCDELG", "ACDEFG", "MNPRQ", "ANPKQ"],
    ... })
    >>> successful, failed = validate_wt_sequence(
    ...     df, name_column="protein", mutation_column="mut", sequence_column="seq",
    ...     wt_identifier="wt", num_workers=2
    ... )
    >>> successful
      protein  mut     seq
    0      P1  A0G  GCDELG
    1      P1  L4F  ACDEFG
    2      P1   wt  ACDELG
    3      P2  K3R  MNPRQ
    4      P2  M0A  ANPKQ
    5      P2   wt  MNPKQ

    See Also
    --------
    validate_wt_sequence_grouped : The per-group worker that validates a single name/protein group.
    """
    _process_protein_group = partial(
        validate_wt_sequence_grouped,
        name_column=name_column,
        mutation_column=mutation_column,
        sequence_column=sequence_column,
        wt_identifier=wt_identifier,
    )

    # Group by protein and process in parallel
    grouped = list(dataset.groupby(name_column, sort=False))

    try:
        results = Parallel(n_jobs=num_workers, backend="loky")(
            delayed(_process_protein_group)(group_data) 
            for group_data in tqdm(grouped, desc="Processing proteins")
        )
    except Exception as e:
        tqdm.write(
            f"Warning: Parallel processing failed, falling back to sequential: {e}"
        )
        # Fallback to sequential processing
        results = []
        for group_data in tqdm(grouped, desc="Processing proteins (sequential)"):
            try:
                result = _process_protein_group(group_data)
                results.append(result)
            except Exception as group_e:
                # Create error entry for this specific group
                protein_name = group_data[0]
                error_row = {
                    name_column: str(protein_name),
                    "error_message": f"Sequential processing error: {type(group_e).__name__}: {str(group_e)}",
                }
                results.append(([error_row], "failed"))

    # Filter out None results and validate structure
    valid_results = []
    invalid_count = 0

    for i, result in enumerate(results):
        if result is None:
            invalid_count += 1
            tqdm.write(f"Warning: Result {i} is None, skipping")
            continue

        if not isinstance(result, tuple) or len(result) != 2:
            invalid_count += 1
            tqdm.write(f"Warning: Result {i} has invalid format, skipping: {result}")
            continue

        rows_list, category = result
        if category not in ("success", "failed"):
            invalid_count += 1
            tqdm.write(
                f"Warning: Result {i} has invalid category '{category}', skipping"
            )
            continue

        if not isinstance(rows_list, list):
            invalid_count += 1
            tqdm.write(f"Warning: Result {i} has invalid rows format, skipping")
            continue

        valid_results.append(result)

    if invalid_count > 0:
        tqdm.write(f"Warning: {invalid_count} invalid results were skipped")

    # Collect all rows
    successful_rows = []
    failed_rows = []

    for rows_list, category in valid_results:
        if category == "success":
            successful_rows.extend(rows_list)
        else:
            failed_rows.extend(rows_list)

    # Convert to DataFrame format
    successful_df = pd.DataFrame(successful_rows) if successful_rows else pd.DataFrame()
    failed_df = pd.DataFrame(failed_rows) if failed_rows else pd.DataFrame()

    tqdm.write(f"Success: {len(successful_df)}, Failed: {len(failed_df)}")

    return successful_df, failed_df


def validate_wt_sequence_grouped(
    group_data: Tuple[Any, pd.DataFrame],
    name_column: str,
    mutation_column: str,
    sequence_column: str,
    wt_identifier: str = "wt",
) -> Tuple[List[Dict[str, Any]], str]:
    """
    Validate one protein/name group by comparing the explicit WT row with
    the WT sequence inferred from each mutant row.
    If explicit WT row is missing, synthesize it using the inferred sequence.
    """
    protein_name, original_group = group_data

    try:
        # Check WT sequence is valid
        wt_rows = original_group[original_group[mutation_column] == wt_identifier]

        # Fail immediately if multiple explicit WT rows exist
        if len(wt_rows) > 1:
            error_row = wt_rows.iloc[0].to_dict()
            error_row["error_message"] = "Multiple explicit WT rows"
            return [error_row], "failed"

        has_wt = not wt_rows.empty

        # Extract original WT sequence
        explicit_wt_seq = (
            str(wt_rows.iloc[0][sequence_column]).strip() if has_wt else None
        )
        wt_row_dict = wt_rows.iloc[0].to_dict() if has_wt else None

        # Extract mutant rows
        mutants = original_group[original_group[mutation_column] != wt_identifier]

        # Handle edge case where no mutants exist
        if mutants.empty and has_wt:
            return [wt_row_dict], "success"
        elif mutants.empty and not has_wt:
            return [], "failed"

        # Infer wild-type sequences
        inferred_wt_seqs = set()
        result_rows = []

        # First, add all original rows to result
        for _, row in mutants.iterrows():
            result_rows.append(row.to_dict())

        # Then, infer WT sequences
        for _, row in mutants.iterrows():
            mut_info = row[mutation_column]
            mut_seq = row[sequence_column]

            # Parse mutation and create sequence
            mutation_set = MutationSet.from_string(mut_info, sep=",", is_zero_based=True)
            sequence = ProteinSequence(str(mut_seq).strip())

            # Infer wild-type sequence by applying inverted mutations
            inverted_mutation_set = invert_mutation_set(mutation_set)
            wt_seq = sequence.apply_mutation(inverted_mutation_set)
            inferred_wt_seqs.add(str(wt_seq))

        # Check consistency among inferred WT sequences
        if len(inferred_wt_seqs) > 1:
            error_row = mutants.iloc[0].to_dict()
            error_row["error_message"] = (
                f"Multiple wildtype sequences inferred for {protein_name}: "
                f"{len(inferred_wt_seqs)}"
            )
            return [error_row], "failed"

        inferred_wt = inferred_wt_seqs.pop()

        # Perform validation checks
        if has_wt:
            # Compare inferred and explicit WT sequences
            if explicit_wt_seq != inferred_wt:
                error_row = mutants.iloc[0].to_dict()
                error_row["error_message"] = (
                    "Inferred WT sequence differs from original WT sequence: "
                    f"{explicit_wt_seq} vs. {inferred_wt}"
                    )
                return [error_row], "failed"
            else:
                result_rows.append(wt_row_dict)
                return result_rows, "success"
        else:
            # Synthesize WT row if missing
            synthetic_wt_row = mutants.iloc[0].to_dict()
            synthetic_wt_row[mutation_column] = wt_identifier
            synthetic_wt_row[sequence_column] = inferred_wt

            # Clear experimental values for synthetic row
            for col in synthetic_wt_row.keys():
                if col not in [name_column, mutation_column, sequence_column]:
                    synthetic_wt_row[col] = None

            result_rows.append(synthetic_wt_row)
            return result_rows, "success"

    except Exception as e:
        # Save error information in first row
        error_row = mutants.iloc[0].to_dict() if "mutants" in locals() and len(mutants) > 0 else {name_column: str(protein_name)}
        error_row["error_message"] = f"{type(e).__name__}: {str(e)}"
        return [error_row], "failed"
