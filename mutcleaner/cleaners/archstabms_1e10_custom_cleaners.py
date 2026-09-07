from __future__ import annotations

from typing import TYPE_CHECKING

import tqdm
import numpy as np
import pandas as pd
from tqdm import tqdm
from ..core.mutation import MutationSet
from ..core.pipeline import multiout_step, pipeline_step

if TYPE_CHECKING:
    from typing import Dict, List, Optional, Tuple

__all__ = [
    "compute_mutations",
    "add_wild_type_sequences_by_library",
    "convert_pairwise_couplings_to_ddg",
]


def __dir__() -> List[str]:
    return __all__


@pipeline_step
def compute_mutations(
    dataset: pd.DataFrame,
    name_column: str = "name",
    WT_column: str = "WT",
    mut_seq: str = "mut_seq",
) -> pd.DataFrame:
    """compute the mutations by the wt_seq and mut_seq and generate mutation column

    Parameters
    ----------
    dataset : pandas.DataFrame
    name_column : str
        Grouping key column.
    WT_column : str
        column used to check whether wt or mut
    mut_seq : str
        Column containing the amino-acid sequence for that row (treated as the
        mutated sequence for mutants and the WT sequence for the WT row).
    """

    def get_mut_info(group):
        # get wt sequence
        name = group.name

        wt_rows = group[group[WT_column] == True]
        if len(wt_rows) == 0:
            return pd.Series([""] * len(group), index=group.index)

        wt_seq = wt_rows[mut_seq].values[0]  # means wt_seq
        wt_array = np.array(list(wt_seq))

        # convert sequence to character matrix
        aa_array = np.array([list(seq) for seq in group[mut_seq]])
        diff_mask = aa_array != wt_array

        # generate mutation info for each sequence
        desc = f"Processing name_column {name}"
        mut_info_list = []
        for i, row in enumerate(tqdm(aa_array, desc=desc)):
            positions = np.where(diff_mask[i])[0]
            if len(positions) == 0:  # WT
                mut_info_list.append("WT")
                continue
            muts = [f"{wt_array[pos]}{pos}{row[pos]}" for pos in positions]
            mut_str = ",".join(muts)
            mut_info_list.append(
                str(MutationSet.from_string(mut_str, is_zero_based=True))
                if mut_str
                else ""
            )
        return pd.Series(mut_info_list, index=group.index)

    dataset["mut_info"] = dataset.groupby(name_column, group_keys=False).apply(
        get_mut_info
    )

    dataset = dataset.drop(columns=WT_column)
    return dataset


@pipeline_step
def add_wild_type_sequences_by_library(
    dataset: pd.DataFrame,
    library_sequences: Dict[int, str],
    library_column: str = "library",
    sequence_column: str = "wt_sequence",
) -> pd.DataFrame:
    """Add wild-type sequences based on library identifiers.

    Parameters
    ----------
    dataset : pd.DataFrame
        Input dataset.
    library_sequences : Dict[int, str]
        Mapping from library identifiers to wild-type sequences.
    library_column : str, default="library"
        Column containing library identifiers.
    sequence_column : str, default="wt_sequence"
        Output column for wild-type sequences.

    Returns
    -------
    pd.DataFrame
        Dataset with the wild-type sequence column added.

    Raises
    ------
    ValueError
        If any library identifier has no corresponding wild-type sequence.
    """
    result = dataset.copy()
    result[sequence_column] = result[library_column].map(library_sequences)

    missing_libraries = (
        result.loc[result[sequence_column].isna(), library_column]
        .drop_duplicates()
        .tolist()
    )
    if missing_libraries:
        raise ValueError(
            f"No wild-type sequence configured for libraries: {missing_libraries}"
        )

    return result


@multiout_step(main="success", failed="failed")
def convert_pairwise_couplings_to_ddg(
    dataset: pd.DataFrame,
    group_columns: Optional[List[str]] = None,
    mutation_column: str = "mut_info",
    label_column: str = "label",
    mutation_separator: str = ",",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Convert pairwise energetic couplings to double-mutant free-energy changes.

    For a double mutant containing mutations A and B, the total free-energy
    change relative to the wild type is calculated as:

    ΔΔG(A,B) = ΔΔG(A) + ΔΔG(B) + ΔΔΔG(A,B)

    Single-mutant labels are left unchanged. Double mutants missing either
    constituent single-mutant label are moved to the failed dataset.

    Parameters
    ----------
    dataset : pd.DataFrame
        Input dataset containing single- and double-mutant energy terms.
    group_columns : Optional[List[str]], default=None
        Columns identifying the protein and biophysical trait within which
        single-mutant labels are matched. Defaults to ``["name"]``.
    mutation_column : str, default="mut_info"
        Column containing mutation descriptions.
    label_column : str, default="label"
        Column containing inferred free-energy terms.
    mutation_separator : str, default=","
        Separator between individual mutations.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        A tuple containing:

        - successful_dataset : pd.DataFrame
            Single mutants and double mutants successfully converted to ΔΔG.
        - failed_dataset : pd.DataFrame
            Double mutants missing at least one constituent single-mutant
            label, with an additional ``error_message`` column.

    Raises
    ------
    ValueError
        If required columns are missing or a single mutation has multiple
        label measurements within the same group.

    Examples
    --------
    >>> df = pd.DataFrame({
    ...     "name": ["proteinA"] * 5,
    ...     "mut_info": [
    ...         "A0G",
    ...         "B1V",
    ...         "A0G,B1V",
    ...         "C2D",
    ...         "A0G,D3E",
    ...     ],
    ...     "label": [0.4, -0.1, 0.2, 0.3, 0.5],
    ... })
    >>> successful, failed = convert_pairwise_couplings_to_ddg(df)
    >>> successful["label"].tolist()
    [0.4, -0.1, 0.5, 0.3]
    >>> failed["mut_info"].tolist()
    ['A0G,D3E']
    """
    group_columns = group_columns or ["name"]

    required_columns = {*group_columns, mutation_column, label_column}
    missing_columns = required_columns - set(dataset.columns)
    if missing_columns:
        raise ValueError(f"Columns not found in dataset: {sorted(missing_columns)}")

    result = dataset.copy()
    failed_messages = {}
    grouper = group_columns[0] if len(group_columns) == 1 else group_columns

    for group_name, group in result.groupby(grouper, sort=False, dropna=False):
        mutation_parts = (
            group[mutation_column]
            .astype("string")
            .str.split(mutation_separator)
            .map(
                lambda mutations: (
                    [mutation.strip() for mutation in mutations]
                    if isinstance(mutations, list)
                    else mutations
                )
            )
        )
        mutation_orders = mutation_parts.str.len()

        single_mask = mutation_orders.eq(1) & group[mutation_column].ne("WT")
        single_mutations = mutation_parts.loc[single_mask].str[0]

        duplicated = single_mutations[single_mutations.duplicated()].unique()
        if len(duplicated):
            raise ValueError(
                f"Duplicated single-mutant labels in group {group_name!r}: "
                f"{duplicated.tolist()}"
            )

        single_labels = pd.Series(
            group.loc[single_mask, label_column].to_numpy(),
            index=single_mutations,
        )

        for index, mutations in mutation_parts.loc[mutation_orders.eq(2)].items():
            missing_mutations = [
                mutation
                for mutation in mutations
                if mutation not in single_labels.index
            ]

            if missing_mutations:
                failed_messages[index] = (
                    "Missing constituent single-mutant label(s): "
                    f"{', '.join(missing_mutations)}"
                )
                continue

            result.at[index, label_column] += single_labels.loc[mutations].sum()

    failed_dataset = result.loc[list(failed_messages)].copy()
    failed_dataset["error_message"] = pd.Series(failed_messages)

    successful_dataset = result.drop(index=failed_dataset.index).copy()

    return successful_dataset, failed_dataset
