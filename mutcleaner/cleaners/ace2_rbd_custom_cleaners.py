from __future__ import annotations

import re
from typing import TYPE_CHECKING

import pandas as pd

from .basic_cleaners import apply_mutations_to_sequences
from ..core.pipeline import multiout_step, pipeline_step

if TYPE_CHECKING:
    from typing import Dict, List, Tuple

__all__ = [
    "prepare_ace2_binding_records",
    "add_reference_sequences_by_target",
    "apply_mutations_preserving_wild_type",
]

def __dir__() -> List[str]:
    """Return exported names.

    Returns
    -------
    List[str]
        Exported names.
    """

    return __all__


def _normalize_target_name(
    raw_target_name: str,
    target_name_aliases: Dict[str, str],
    reference_sequences: Dict[str, str],
) -> str:
    """Resolve one raw target label to a canonical target name.

    Parameters
    ----------
    raw_target_name : str
        Target name observed in the source data.
    target_name_aliases : Dict[str, str]
        Alias-to-canonical target mapping.
    reference_sequences : Dict[str, str]
        Canonical target reference sequences.

    Returns
    -------
    str
        Canonical target name.

    Raises
    ------
    ValueError
        If the target cannot be resolved.
    """

    def normalize_key(value: str) -> str:
        return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())

    normalized_target_name = str(raw_target_name).strip()
    if normalized_target_name in reference_sequences:
        return normalized_target_name

    normalized_key = normalize_key(normalized_target_name)
    normalized_aliases = {
        normalize_key(alias): canonical_name
        for alias, canonical_name in target_name_aliases.items()
    }
    if normalized_key in normalized_aliases:
        return normalized_aliases[normalized_key]

    direct_matches = [
        canonical_name
        for canonical_name in reference_sequences
        if normalize_key(canonical_name) == normalized_key
    ]
    if len(direct_matches) == 1:
        return direct_matches[0]

    raise ValueError(
        f"Unknown RBD target '{raw_target_name}'. Add it to target_name_aliases "
        "or reference_sequences in the cleaner config."
    )


@pipeline_step
def prepare_ace2_binding_records(
    dataset: pd.DataFrame,
    reference_sequences: Dict[str, str],
    target_name_aliases: Dict[str, str],
    column_mapping: Dict[str, str],
    target_column: str = "target",
    mutation_column: str = "aa_substitutions",
    label_column: str = "log10Ka",
    variant_class_column: str = "variant_class",
    default_target_name: str = "Wuhan-Hu-1",
    output_name_column: str = "name",
    output_mutation_column: str = "mut_info",
    output_label_column: str = "label",
) -> pd.DataFrame:
    """Standardize one ACE2 table into cleaner-ready columns.

    Parameters
    ----------
    dataset : pd.DataFrame
        Raw ACE2 binding-score dataframe.
    reference_sequences : Dict[str, str]
        Canonical RBD target reference sequences.
    target_name_aliases : Dict[str, str]
        RBD target alias-to-canonical-name mapping.
    column_mapping : Dict[str, str]
        Raw-to-standard column mapping. Known columns are renamed when present
        while unrelated columns are preserved.
    target_column : str, default="target"
        Source target-name column.
    mutation_column : str, default="aa_substitutions"
        Source mutation column.
    label_column : str, default="log10Ka"
        Source label column.
    variant_class_column : str, default="variant_class"
        Source variant-class column.
    default_target_name : str, default="Wuhan-Hu-1"
        Fallback target name if the source table omits ``target_column``.
    output_name_column : str, default="name"
        Output canonical target-name column.
    output_mutation_column : str, default="mut_info"
        Output mutation column.
    output_label_column : str, default="label"
        Output label column.

    Returns
    -------
    pd.DataFrame
        Standardized dataframe ready for downstream cleaning.

    Raises
    ------
    ValueError
        If required source columns are missing.
    """

    def normalize_bool_series(series: pd.Series) -> pd.Series:
        if pd.api.types.is_bool_dtype(series):
            return series.fillna(False)

        mapping = {
            "true": True,
            "false": False,
            "1": True,
            "0": False,
            "yes": True,
            "no": False,
            "t": True,
            "f": False,
        }
        return series.astype(str).str.strip().str.lower().map(mapping).fillna(False)

    result = dataset.copy()
    available_column_mapping = {
        raw_name: standardized_name
        for raw_name, standardized_name in column_mapping.items()
        if raw_name in result.columns and raw_name != standardized_name
    }
    if available_column_mapping:
        result = result.rename(columns=available_column_mapping)

    if mutation_column not in result.columns:
        raise ValueError(
            f"{mutation_column} column is missing in ACE2 raw input data"
        )
    if label_column not in result.columns:
        raise ValueError(f"{label_column} column is missing in ACE2 raw input data")

    if target_column not in result.columns:
        result[target_column] = default_target_name
    else:
        result[target_column] = (
            result[target_column].fillna(default_target_name).astype(str)
        )
    if variant_class_column not in result.columns:
        result[variant_class_column] = ""
    if "pass_pre_count_filter" in result.columns:
        result = result.loc[
            normalize_bool_series(result["pass_pre_count_filter"])
        ].copy()
    if "pass_ACE2bind_expr_filter" in result.columns:
        result = result.loc[
            normalize_bool_series(result["pass_ACE2bind_expr_filter"])
        ].copy()
    result[mutation_column] = result[mutation_column].fillna("").astype(str).str.strip()
    result = result.loc[
        ~result[mutation_column].str.contains(r"\*", regex=True, na=False)
    ].copy()

    result[output_name_column] = result[target_column].map(
        lambda value: _normalize_target_name(
            raw_target_name=str(value),
            target_name_aliases=target_name_aliases,
            reference_sequences=reference_sequences,
        )
    )
    result[output_mutation_column] = result[mutation_column]
    result[output_label_column] = pd.to_numeric(result[label_column], errors="coerce")

    standardized = result[
        [
            output_name_column,
            output_mutation_column,
            variant_class_column,
            output_label_column,
        ]
    ].rename(columns={variant_class_column: "variant_class"})
    standardized = standardized.loc[standardized[output_label_column].notna()].copy()

    mutation_text = standardized[output_mutation_column].fillna("").astype(str).str.strip()
    variant_class = standardized["variant_class"].fillna("").astype(str).str.strip().str.lower()
    synonymous_mask = variant_class.eq("synonymous")
    stop_mask = variant_class.eq("stop")
    deletion_mask = mutation_text.str.contains("-", regex=False)
    empty_mutation_mask = mutation_text.eq("")

    wt_mask = variant_class.eq("wildtype")
    standardized.loc[wt_mask, output_mutation_column] = "WT"
    keep_mask = wt_mask | ~(
        synonymous_mask | stop_mask | deletion_mask | empty_mutation_mask
    )
    standardized = standardized.loc[keep_mask].copy()
    return standardized.reset_index(drop=True)


@pipeline_step
def add_reference_sequences_by_target(
    dataset: pd.DataFrame,
    reference_sequences: Dict[str, str],
    name_column: str = "name",
    sequence_column: str = "sequence",
) -> pd.DataFrame:
    """Attach reference sequences to standardized ACE2 rows.

    Parameters
    ----------
    dataset : pd.DataFrame
        Standardized ACE2 dataframe.
    reference_sequences : Dict[str, str]
        Canonical RBD target reference sequences.
    name_column : str, default="name"
        Column containing canonical target names.
    sequence_column : str, default="sequence"
        Output sequence column.

    Returns
    -------
    pd.DataFrame
        Dataframe with reference sequences attached.

    Raises
    ------
    ValueError
        If a target lacks a configured reference sequence.
    """

    result = dataset.copy()
    result[sequence_column] = result[name_column].map(reference_sequences)
    missing_targets = sorted(
        {
            str(target)
            for target in result.loc[result[sequence_column].isna(), name_column].dropna()
        }
    )
    if missing_targets:
        raise ValueError(
            "Missing reference sequences for targets: "
            + ", ".join(missing_targets)
        )
    return result


@multiout_step(main="success", failed="failed")
def apply_mutations_preserving_wild_type(
    dataset: pd.DataFrame,
    sequence_column: str = "sequence",
    name_column: str = "name",
    mutation_column: str = "mut_info",
    wt_identifier: str = "WT",
    mutation_sep: str = ",",
    is_zero_based: bool = True,
    sequence_type: str = "protein",
    num_workers: int = 4,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Materialize ACE2 mutant sequences while keeping WT rows intact.

    Parameters
    ----------
    dataset : pd.DataFrame
        Standardized ACE2 dataframe with attached reference sequences.
    sequence_column : str, default="sequence"
        Column containing reference sequences.
    name_column : str, default="name"
        Column containing canonical target names.
    mutation_column : str, default="mut_info"
        Column containing mutation strings.
    wt_identifier : str, default="WT"
        Wild-type marker used in ``mutation_column``.
    mutation_sep : str, default=","
        Mutation separator used by the downstream mutation applier.
    is_zero_based : bool, default=True
        Whether mutation positions are already zero-based.
    sequence_type : str, default="protein"
        Sequence type passed to the basic mutation applier.
    num_workers : int, default=4
        Number of worker processes.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        Successful rows with ``mut_seq`` and failed rows with ``error_message``.
    """

    mutation_text = dataset[mutation_column].fillna("").astype(str).str.upper()
    wt_mask = mutation_text.eq(str(wt_identifier).upper())

    wt_rows = dataset.loc[wt_mask].copy()
    if not wt_rows.empty:
        wt_rows["mut_seq"] = wt_rows[sequence_column].astype(str)

    non_wt_rows = dataset.loc[~wt_mask].copy()
    if non_wt_rows.empty:
        return wt_rows.reset_index(drop=True), pd.DataFrame(
            columns=[*dataset.columns, "error_message"]
        )

    mutation_result = apply_mutations_to_sequences(
        non_wt_rows,
        sequence_column=sequence_column,
        name_column=name_column,
        mutation_column=mutation_column,
        mutation_sep=mutation_sep,
        is_zero_based=is_zero_based,
        sequence_type=sequence_type,
        num_workers=num_workers,
    )
    successful_rows = pd.concat([wt_rows, mutation_result.main], axis=0).sort_index(
        kind="stable"
    )
    return (
        successful_rows.reset_index(drop=True),
        mutation_result.side.get("failed", pd.DataFrame()).reset_index(drop=True),
    )
