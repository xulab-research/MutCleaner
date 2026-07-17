# cleaners/human_domainome_custom_cleaners.py
from __future__ import annotations

import pandas as pd
from functools import partial
from joblib import Parallel, delayed
from tqdm import tqdm
from typing import TYPE_CHECKING

from .basic_cleaners import split_columns, merge_columns
from ..core.pipeline import multiout_step, pipeline_step

if TYPE_CHECKING:
    from typing import Dict, List, Optional, Tuple

__all__ = [
    "generate_mutation_strings",
    "process_domain_positions",
    "extract_domain_sequences",
]


def __dir__() -> List[str]:
    return __all__


@pipeline_step
def generate_mutation_strings(
    dataset: pd.DataFrame,
    name_column: str = "domain_ID",
    wt_aa_column: str = "wt_aa",
    mut_aa_column: str = "mut_aa",
    aa_pos_column: str = "pos",
) -> pd.DataFrame:
    """
    Generate per-row mutation notations from wild-type amino acids, positions, and a
    domain identifier that encodes a sequence offset.

    This step expects the identifier column (``name_column``) to be formatted as
    ``"<Uniprot_ID>_<Pfam_ID>_<sequence_offset>"``. It splits that column, computes a
    relative residue position as ``pos - sequence_offset``, builds a simple mutation
    string ``<wt_aa><relative_pos><mut_aa>`` (e.g., ``A15K``), stores it in
    ``mut_info``. Helper columns introduced during the transformation are dropped
    before returning.

    Parameters
    ----------
    dataset : pandas.DataFrame
        Input table containing at least the identifier, wild-type AA, mutant AA,
        and absolute position columns.
    name_column : str, default "domain_ID"
        Column holding the domain identifier in the form
        ``Uniprot_ID_Pfam_ID_sequence_offset``.
    wt_aa_column : str, default "wt_aa"
        Column with the wild-type amino acid (single-letter code).
    mut_aa_column : str, default "mut_aa"
        Column with the mutant amino acid (single-letter code).
    aa_pos_column : str, default "pos"
        Column with the absolute residue position (integer).

    Returns
    -------
    pd.DataFrame
        The input dataframe with:
        * A new column ``mut_info`` containing strings like ``A15K``;
        * All other columns are preserved.

    Notes
    -----
    - The relative position is computed as ``pos - sequence_offset``; interpret it
    as 0-based or 1-based according to how ``sequence_offset`` is defined in your
    data.
    - The function relies on utility helpers ``split_columns`` and ``merge_columns``.
    It also treats ``sequence_offset`` and ``pos`` as integers.

    Raises
    ------
    KeyError
        If any of the required columns are missing.
    ValueError
        If position or offset values cannot be interpreted as integers.

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ...     "domain_ID": ["P12345_PF00001_10"],
    ...     "wt_aa": ["A"],
    ...     "mut_aa": ["K"],
    ...     "pos": [25],
    ... })
    >>> generate_mutation_strings(df)[["domain_ID", "mut_info"]]
    Splitting column 'domain_ID' into ['Uniprot_ID', 'Pfam_ID', 'sequence_offset']...
    Splitting using separator: '_'
    Successfully created split columns: ['Uniprot_ID', 'Pfam_ID', 'sequence_offset']
    Merging columns ['wt_aa', 'aa_pos', 'mut_aa'] into 'mut_info'...
    Successfully created merged column 'mut_info'
    Merging columns ['Uniprot_ID', 'Pfam_ID', 'aa_pos'] into 'domain_ID'...
    Successfully created merged column 'domain_ID'
               domain_ID mut_info
    0  P12345_PF00001_10     A15K
    """
    # Separate 'name' column to get sequence offset
    dataset = split_columns(
        dataset,
        column_to_split=name_column,
        new_column_names=["Uniprot_ID", "Pfam_ID", "sequence_offset"],
        separator="_",
    )
    # Convert columns to appropriate data types
    dataset.loc[:, ["sequence_offset", aa_pos_column]] = (
        dataset.loc[:, ["sequence_offset", aa_pos_column]]
        .apply(pd.to_numeric, errors="coerce")
        .astype("Int64")
    )
    # Generate mutation strings
    dataset["aa_pos"] = (dataset[aa_pos_column] - dataset["sequence_offset"]).astype(
        str
    )
    dataset = merge_columns(
        dataset,
        columns_to_merge=[wt_aa_column, "aa_pos", mut_aa_column],
        new_column_name="mut_info",
        separator="",
    )

    dataset = dataset.drop(
        columns=["Uniprot_ID", "Pfam_ID", "sequence_offset", "aa_pos"]
    )
    return dataset


# ==================================================================
# ======== Pipeline Functions Used in Supplementary Table 4 ========
# ==================================================================
@multiout_step(main="success", failed="failed")
def process_domain_positions(
    dataset: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Process domain position information from PFAM entries

    This function extracts position information from PFAM entries and calculates
    relative mutation positions. It handles parsing errors by separating failed records.

    Parameters
    ----------
    dataset : pd.DataFrame
        Dataset with PFAM_entry column containing position information

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        (successful_dataset, failed_dataset) - datasets with and without parsing errors

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ...     'PFAM_entry': ['PF00001/10-100', 'PF00002/20-200', 'invalid_entry'],
    ...     'pos': [15, 25, 30],
    ...     'wt_aa': ['A', 'C', 'D'],
    ...     'mut_aa': ['K', 'Y', 'E']
    ... })
    >>> successful, failed = process_domain_positions(df)
    >>> print(len(successful))  # Should be 2
    2
    >>> print(len(failed))  # Should be 1
    1
    """
    tqdm.write("Extracting peptide position information...")

    # Create a copy to avoid modifying the original
    result_dataset = dataset.copy()
    result_dataset["error_message"] = None

    # Extract start and end positions from PFAM_entry (format: "PFAM_ENTRY/start-end")
    try:
        position_info = result_dataset["PFAM_entry"].str.extract(r"/(\d+)-(\d+)")

        # Track which entries failed to parse
        parse_failed = position_info.isnull().any(axis=1)
        result_dataset.loc[parse_failed, "error_message"] = (
            "Failed to parse PFAM_entry position information"
        )

        # Process successful entries
        success_mask = ~parse_failed

        if success_mask.any():
            result_dataset.loc[success_mask, "start_pos"] = (
                position_info.loc[success_mask, 0].astype(int) - 1
            )  # Convert to 0-based
            result_dataset.loc[success_mask, "end_pos"] = position_info.loc[
                success_mask, 1
            ].astype(int)

            # Convert absolute position to 0-based
            result_dataset.loc[success_mask, "pos"] = (
                result_dataset.loc[success_mask, "pos"] - 1
            )

            # Calculate relative position within the domain
            result_dataset.loc[success_mask, "mut_rel_pos"] = (
                result_dataset.loc[success_mask, "pos"]
                - result_dataset.loc[success_mask, "start_pos"]
            )

            # Generate mutation info using relative position
            result_dataset.loc[success_mask, "mut_info"] = (
                result_dataset.loc[success_mask, "wt_aa"]
                + result_dataset.loc[success_mask, "mut_rel_pos"]
                .astype(int)
                .astype(str)
                + result_dataset.loc[success_mask, "mut_aa"]
            )

    except Exception as e:
        # If something goes wrong, mark all as failed
        result_dataset["error_message"] = f"Error processing positions: {str(e)}"
        success_mask = pd.Series([False] * len(result_dataset))

    # Separate successful and failed datasets
    successful_dataset = result_dataset[success_mask].drop(columns=["error_message"])
    failed_dataset = result_dataset[~success_mask]

    successful_dataset[["start_pos", "end_pos"]] = successful_dataset[
        ["start_pos", "end_pos"]
    ].astype(int)

    tqdm.write(
        f"Position processing: {len(successful_dataset)} successful, {len(failed_dataset)} failed"
    )

    return successful_dataset, failed_dataset


@multiout_step(main="success", failed="failed")
def extract_domain_sequences(
    dataset: pd.DataFrame,
    sequence_column: str = "sequence",
    start_pos_column: str = "start_pos",
    end_pos_column: str = "end_pos",
    num_workers: int = 4,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Extract domain sequences from full sequences using position information

    This function extracts domain subsequences based on start and end positions.
    Records with invalid positions or missing sequences are separated into the failed dataset.

    Parameters
    ----------
    dataset : pd.DataFrame
        Dataset with full sequences and position information
    sequence_column : str, default='sequence'
        Column containing full wild-type sequences
    start_pos_column : str, default='start_pos'
        Column containing domain start positions (0-based)
    end_pos_column : str, default='end_pos'
        Column containing domain end positions
    num_workers : int, default=4
        Number of parallel workers for processing, set to -1 for all available CPUs

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        (successful_dataset, failed_dataset) - datasets with and without extraction errors

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ...     'name': ['prot1', 'prot2', 'prot3'],
    ...     'sequence': ['ABCDEFGHIJ', 'KLMNOPQRST', None],
    ...     'start_pos': [2, 0, 5],
    ...     'end_pos': [7, 4, 10]
    ... })
    >>> successful, failed = extract_domain_sequences(df)
    >>> print(successful['sequence'].tolist())
    ['CDEFG', 'KLMN']
    >>> print(len(failed))  # Should be 1 (the None sequence)
    1
    """
    tqdm.write("Extracting domain sequences from full sequences...")

    # Validate required columns exist
    required_columns = [sequence_column, start_pos_column, end_pos_column]
    missing_columns = [col for col in required_columns if col not in dataset.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    # Create a copy to avoid modifying the original
    result_dataset = dataset.copy()

    # Prepare partial function with fixed parameters
    _extract_domain = partial(
        _extract_single_domain,
        dataset_columns=dataset.columns,
        sequence_column=sequence_column,
        start_pos_column=start_pos_column,
        end_pos_column=end_pos_column,
    )

    # Parallel processing
    rows = dataset.itertuples(index=False, name=None)
    results = Parallel(n_jobs=num_workers, backend="loky")(
        delayed(_extract_domain)(row)
        for row in tqdm(rows, total=len(dataset), desc="Extracting domains")
    )

    # Separate domain sequences and error messages
    domain_sequences, error_messages = map(list, zip(*results))

    result_dataset["domain_sequence"] = domain_sequences
    result_dataset["error_message"] = error_messages

    # Success mask is where we have domain sequences
    success_mask = pd.notnull(result_dataset["domain_sequence"])

    # For successful extractions, replace the original sequence with the domain sequence
    result_dataset.loc[success_mask, sequence_column] = result_dataset.loc[
        success_mask, "domain_sequence"
    ]

    # Separate successful and failed datasets
    successful_dataset = result_dataset[success_mask].drop(
        columns=["error_message", "domain_sequence", start_pos_column, end_pos_column]
    )
    failed_dataset = result_dataset[~success_mask].drop(columns=["domain_sequence"])

    tqdm.write(
        f"Domain extraction: {len(successful_dataset)} successful, {len(failed_dataset)} failed"
    )

    return successful_dataset, failed_dataset


def _extract_single_domain(
    row_data: Tuple,
    dataset_columns: pd.Index,
    sequence_column: str,
    start_pos_column: str,
    end_pos_column: str,
) -> Tuple[Optional[str], Optional[str]]:
    """Extract domain sequence for a single row

    Returns
    -------
    Tuple[Optional[str], Optional[str]]
        (domain_sequence, error_message)
    """
    # Convert tuple to dict for easier access
    row = dict(zip(dataset_columns, row_data))

    try:
        full_seq = row.get(sequence_column)
        start_pos = row.get(start_pos_column)
        end_pos = row.get(end_pos_column)

        if pd.isnull(full_seq) or full_seq is None:
            return None, "Missing sequence"

        if pd.isnull(start_pos) or pd.isnull(end_pos):
            return None, "Missing position information"

        start_pos = int(start_pos)
        end_pos = int(end_pos)

        if start_pos < 0:
            return None, f"Invalid start position: {start_pos} < 0"

        if end_pos > len(full_seq):
            return (
                None,
                f"End position {end_pos} exceeds sequence length {len(full_seq)}",
            )

        if start_pos >= end_pos:
            return None, f"Invalid position range: start={start_pos} >= end={end_pos}"

        if start_pos >= len(full_seq):
            return (
                None,
                f"Start position {start_pos} exceeds sequence length {len(full_seq)}",
            )

        return full_seq[start_pos:end_pos], None

    except Exception as e:
        return None, f"Extraction error: {str(e)}"
