from __future__ import annotations

import logging
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from .base_config import BaseCleanerConfig
from .basic_cleaners import (
    apply_mutations_to_sequences,
    average_labels_by_name,
    convert_to_mutation_dataset_format,
    convert_data_types,
    extract_and_rename_columns,
    filter_and_clean_data,
    merge_columns,
    read_dataset,
    subtract_labels_by_wt,
    validate_mutations,
)
from .rbd_custom_cleaner import (
    add_reference_sequences_by_target,
    mark_wild_type_in_mut_info,
    standardize_rbd_target_names,
)
from ..core.dataset import MutationDataset
from ..core.pipeline import Pipeline, create_pipeline

if TYPE_CHECKING:
    from typing import Any, Dict, List, Optional, Tuple, Union

__all__ = [
    "RBDAntibodyCleanerConfig",
    "create_rbd_antibody_cleaner",
    "clean_rbd_antibody_dataset",
]

logger = logging.getLogger(__name__)

DEFAULT_RBD_REFERENCE_SEQUENCES = {
    "Wuhan-Hu-1": "NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKST",
}

DEFAULT_RBD_TARGET_NAME_ALIASES = {
    "SARS-CoV-2": "Wuhan-Hu-1",
    "SARS_CoV_2": "Wuhan-Hu-1",
    "Wuhan_Hu_1": "Wuhan-Hu-1",
}


def __dir__() -> List[str]:
    """Return exported names."""

    return __all__


@dataclass
class RBDAntibodyCleanerConfig(BaseCleanerConfig):
    """Configuration for the RBD antibody cleaner.

    Attributes
    ----------
    reference_sequences : Dict[str, str]
        Canonical RBD reference sequences keyed by target/reference name.
    target_name_aliases : Dict[str, str]
        RBD target alias-to-canonical-name mapping.
    column_mapping : Dict[str, str]
        Mapping from raw source column names to standardized column names
        consumed by the RBD antibody cleaner.
    filters : Dict[str, Any]
        Filter conditions applied before mutation processing.
    drop_na_columns : List[str]
        Columns that must be present after filtering.
    type_conversions : Dict[str, str]
        Data type conversion specifications for label columns.
    fallback_reference_sequence : Optional[str]
        Reference sequence used when a row has an unknown reference name.
    require_known_reference_sequence : bool
        Whether reference names must be found in ``reference_sequences``.
    validate_mut_workers : int
        Worker count for mutation validation.
    process_workers : int
        Worker count for sequence materialization.
    label_columns : List[str]
        Label columns retained through the pipeline.
    primary_label_column : str
        Label column written into the final ``MutationDataset``.
    pipeline_name : str
        Pipeline name.
    """

    reference_sequences: Dict[str, str] = field(default_factory=lambda: deepcopy(DEFAULT_RBD_REFERENCE_SEQUENCES))

    target_name_aliases: Dict[str, str] = field(default_factory=lambda: deepcopy(DEFAULT_RBD_TARGET_NAME_ALIASES))

    column_mapping: Dict[str, str] = field(
        default_factory=lambda: {
            "name": "antibody_name",
            "target": "reference_id",
            "score": "label",
            "aa_substitutions": "mut_info",
            "variant_class": "variant_class",
            "n_aa_substitutions": "n_aa_substitutions",
            "pass_pre_count_filter": "pass_pre_count_filter",
            "pass_ACE2bind_expr_filter": "pass_ACE2bind_expr_filter",
        }
    )

    filters: Dict[str, Any] = field(
        default_factory=lambda: {
            "pass_pre_count_filter": True,
            "pass_ACE2bind_expr_filter": True,
        }
    )

    drop_na_columns: List[str] = field(default_factory=lambda: ["reference_id", "label"])

    type_conversions: Dict[str, str] = field(default_factory=lambda: {"label": "float"})

    fallback_reference_sequence: Optional[str] = None

    require_known_reference_sequence: bool = True

    validate_mut_workers: int = 16

    process_workers: int = 16

    label_columns: List[str] = field(default_factory=lambda: ["label"])

    primary_label_column: str = "label"

    pipeline_name: str = "RBDAntibody pipeline"

    def validate(self) -> None:
        super().validate()

        if self.require_known_reference_sequence and not self.reference_sequences:
            raise ValueError("reference_sequences cannot be empty when " "require_known_reference_sequence=True")
        if not self.require_known_reference_sequence and not self.reference_sequences and self.fallback_reference_sequence is None:
            raise ValueError("Provide reference_sequences or fallback_reference_sequence")
        if not self.label_columns:
            raise ValueError("label_columns cannot be empty")
        if self.primary_label_column not in self.label_columns:
            raise ValueError(f"primary_label_column '{self.primary_label_column}' must be in label_columns {self.label_columns}")

        required_standard_columns = {
            "antibody_name",
            "reference_id",
            "label",
            "mut_info",
            "variant_class",
        }
        missing_standard_columns = required_standard_columns - set(self.column_mapping.values())
        if missing_standard_columns:
            raise ValueError("column_mapping must provide standardized columns " f"{sorted(required_standard_columns)}, missing {sorted(missing_standard_columns)}")

        for target_name, sequence in self.reference_sequences.items():
            if len(str(sequence).strip()) <= 0:
                raise ValueError(f"Reference sequence for target '{target_name}' cannot be empty")

        if self.fallback_reference_sequence is not None and len(str(self.fallback_reference_sequence).strip()) <= 0:
            raise ValueError("fallback_reference_sequence cannot be empty when provided")


def create_rbd_antibody_cleaner(
    dataset_or_path: Optional[Union[pd.DataFrame, str, Path]] = None,
    config: Optional[Union[RBDAntibodyCleanerConfig, Dict[str, Any], str, Path]] = None,
) -> Pipeline:
    """Create the RBD antibody cleaning pipeline."""

    default_config = RBDAntibodyCleanerConfig()
    if config is None:
        final_config = default_config
    elif isinstance(config, RBDAntibodyCleanerConfig):
        final_config = config
    elif isinstance(config, dict):
        final_config = default_config.merge(config)
    elif isinstance(config, (str, Path)):
        final_config = RBDAntibodyCleanerConfig.from_json(config)
    else:
        raise TypeError(f"config must be RBDAntibodyCleanerConfig, dict, str, Path or None, got {type(config)}")

    antibody_name_column = final_config.column_mapping.get("name", "name")
    reference_name_column = final_config.column_mapping.get("target", "target")
    mutation_column = final_config.column_mapping.get("aa_substitutions", "aa_substitutions")
    variant_class_column = final_config.column_mapping.get("variant_class", "variant_class")

    logger.info(
        "RBD antibody dataset will be cleaned with pipeline: %s",
        final_config.pipeline_name,
    )
    logger.debug("Configuration:\n%s", final_config.get_summary())

    pipeline = create_pipeline(dataset_or_path, final_config.pipeline_name)
    pipeline = (
        pipeline.delayed_then(
            extract_and_rename_columns,
            column_mapping=final_config.column_mapping,
        )
        .delayed_then(
            convert_data_types,
            type_conversions=final_config.type_conversions,
        )
        .delayed_then(
            filter_and_clean_data,
            filters=final_config.filters,
            drop_na_columns=final_config.drop_na_columns,
        )
        .delayed_then(
            standardize_rbd_target_names,
            target_name_aliases=final_config.target_name_aliases,
            name_column=reference_name_column,
        )
        .delayed_then(
            mark_wild_type_in_mut_info,
            mutation_column=mutation_column,
            variant_class_column=variant_class_column,
        )
        .delayed_then(
            validate_mutations,
            mutation_column=mutation_column,
            format_mutations=True,
            mutation_sep=" ",
            is_zero_based=False,
            exclude_patterns="WT",
            cache_results=False,
            num_workers=final_config.validate_mut_workers,
        )
        .delayed_then(
            average_labels_by_name,
            name_columns=[reference_name_column, antibody_name_column, mutation_column],
            label_columns=final_config.label_columns,
        )
        # `group_name` is solely an internal WT-subtraction key used for multiple references and antibodies.
        # For the current antibody tables one file is usually one reference,
        # so `antibody_name` would often be enough on its own,
        # but we keep the combined key to preserve the existing downstream grouping behavior.
        .delayed_then(
            merge_columns,
            columns_to_merge=[reference_name_column, antibody_name_column],
            new_column_name="group_name",
            separator="_",
        )
        .delayed_then(
            subtract_labels_by_wt,
            name_column="group_name",
            label_columns=final_config.label_columns,
            mutation_column=mutation_column,
            wt_identifier="WT",
            in_place=True,
            drop_wt_row=True,
        )
        .delayed_then(
            add_reference_sequences_by_target,
            reference_sequences=final_config.reference_sequences,
            name_column=reference_name_column,
            sequence_column="sequence",
            fallback_reference_sequence=final_config.fallback_reference_sequence,
        )
        .delayed_then(
            apply_mutations_to_sequences,
            sequence_column="sequence",
            name_column="group_name",
            mutation_column=mutation_column,
            mutation_sep=",",
            is_zero_based=True,
            sequence_type="protein",
            num_workers=final_config.process_workers,
        )
        .delayed_then(
            convert_to_mutation_dataset_format,
            name_column=antibody_name_column,
            mutation_column=mutation_column,
            sequence_column="sequence",
            mutated_sequence_column="mut_seq",
            label_column=final_config.primary_label_column,
            include_wild_type=False,
            is_zero_based=True,
        )
    )

    if isinstance(dataset_or_path, (str, Path)):
        pipeline.add_delayed_step(read_dataset, 0, file_format="csv")
    elif dataset_or_path is not None and not isinstance(dataset_or_path, pd.DataFrame):
        raise TypeError(f"dataset_or_path must be pd.DataFrame, str, Path, or None, got {type(dataset_or_path)}")

    return pipeline


def clean_rbd_antibody_dataset(
    pipeline: Pipeline,
) -> Tuple[Pipeline, MutationDataset]:
    """Execute the RBD antibody cleaning pipeline."""

    pipeline.execute()
    dataset_df, reference_sequences = pipeline.data
    dataset = MutationDataset.from_dataframe(dataset_df, reference_sequences)
    logger.info(
        "Successfully cleaned RBD antibody dataset: %s mutation sets from %s references",
        len(dataset.mutation_sets),
        len(dataset.reference_sequences),
    )
    return pipeline, dataset
