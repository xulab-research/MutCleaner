from __future__ import annotations
from typing import TYPE_CHECKING

import pandas as pd

from ..core.pipeline import pipeline_step

if TYPE_CHECKING:
    from typing import Dict, List, Optional

__all__ = [
    "standardize_rbd_target_names",
    "mark_wild_type_in_mut_info",
    "add_reference_sequences_by_target",
]


def __dir__() -> List[str]:
    return __all__


@pipeline_step
def standardize_rbd_target_names(
    dataset: pd.DataFrame,
    target_name_aliases: Dict[str, str],
    name_column: str = "name",
) -> pd.DataFrame:
    """Canonicalize RBD target/reference names.

    Parameters
    ----------
    dataset : pd.DataFrame
        Input RBD dataset.
    target_name_aliases : Dict[str, str]
        Mapping from alias names to canonical target names.
    name_column : str, default="name"
        Column containing the target/reference name.

    Returns
    -------
    pd.DataFrame
        Dataset with canonicalized target/reference names.
    """

    result = dataset.copy()
    name_values = result[name_column].astype("string").str.strip()
    result[name_column] = name_values.map(target_name_aliases).fillna(name_values)
    return result.reset_index(drop=True)


@pipeline_step
def mark_wild_type_in_mut_info(
    dataset: pd.DataFrame,
    mutation_column: str = "mut_info",
    variant_class_column: str = "variant_class",
) -> pd.DataFrame:
    """Mark wild-type RBD records in the mutation column.

    Parameters
    ----------
    dataset : pd.DataFrame
        Input RBD dataset.
    mutation_column : str, default="mut_info"
        Column containing mutation descriptions.
    variant_class_column : str, default="variant_class"
        Column indicating whether a row is wild type.

    Returns
    -------
    pd.DataFrame
        Dataset where wild-type rows have ``mutation_column`` set to ``"WT"``.
    """

    result = dataset.copy()
    variant_class = result[variant_class_column].astype("string").str.strip().str.lower()
    wt_mask = variant_class.eq("wildtype").fillna(False)
    result.loc[wt_mask, mutation_column] = "WT"
    return result.reset_index(drop=True)


@pipeline_step
def add_reference_sequences_by_target(
    dataset: pd.DataFrame,
    reference_sequences: Dict[str, str],
    name_column: str = "name",
    sequence_column: str = "sequence",
    fallback_reference_sequence: Optional[str] = None,
) -> pd.DataFrame:
    """Attach reference sequences to standardized RBD rows.

    Parameters
    ----------
    dataset : pd.DataFrame
        Input RBD dataset.
    reference_sequences : Dict[str, str]
        Mapping from target/reference names to RBD reference sequences.
    name_column : str, default="name"
        Column containing target/reference names.
    sequence_column : str, default="sequence"
        Output column for reference sequences.
    fallback_reference_sequence : Optional[str], default=None
        Sequence used when a target/reference name is missing from
        ``reference_sequences``.

    Returns
    -------
    pd.DataFrame
        Dataset with the reference sequence column attached.
    """

    result = dataset.copy()
    result[sequence_column] = result[name_column].map(reference_sequences)
    if fallback_reference_sequence is not None:
        result[sequence_column] = result[sequence_column].fillna(str(fallback_reference_sequence).strip())
    return result
