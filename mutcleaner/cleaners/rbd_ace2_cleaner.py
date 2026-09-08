from __future__ import annotations

import logging
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from .rbd_custom_cleaner import (
    add_reference_sequences_by_target,
    mark_wild_type_in_mut_info,
    standardize_rbd_target_names,
)
from .base_config import BaseCleanerConfig
from .basic_cleaners import (
    apply_mutations_to_sequences,
    average_labels_by_name,
    convert_to_mutation_dataset_format,
    convert_data_types,
    extract_and_rename_columns,
    filter_and_clean_data,
    read_dataset,
    subtract_labels_by_wt,
    validate_mutations,
)
from ..core.dataset import MutationDataset
from ..core.pipeline import Pipeline, create_pipeline

if TYPE_CHECKING:
    from typing import Any, Dict, List, Optional, Tuple, Union

__all__ = [
    "RBDACE2CleanerConfig",
    "create_rbd_ace2_cleaner",
    "clean_rbd_ace2_dataset",
]

logger = logging.getLogger(__name__)

DEFAULT_RBD_REFERENCE_SEQUENCES = {
    "Wuhan-Hu-1": "NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKST",
    "Alpha": "NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTYGVGYQPYRVVVLSFELLHAPATVCGPKKST",
    "Beta": "NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGNIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVKGFNCYFPLQSYGFQPTYGVGYQPYRVVVLSFELLHAPATVCGPKKST",
    "Eta": "NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVKGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKST",
    "Delta": "NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYRYRLFRKSNLKPFERDISTEIYQAGSKPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKST",
    "Omicron_BA1": "NITNLCPFDEVFNATRFASVYAWNRKRISNCVADYSVLYNLAPFFTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGNIADYNYKLPDDFTGCVIAWNSNKLDSKVSGNYNYLYRLFRKSNLKPFERDISTEIYQAGNKPCNGVAGFNCYFPLRSYSFRPTYGVGHQPYRVVVLSFELLHAPATVCGPKKST",
    "Omicron_BA2": "NITNLCPFDEVFNATRFASVYAWNRKRISNCVADYSVLYNFAPFFAFKCYGVSPTKLNDLCFTNVYADSFVIRGNEVSQIAPGQTGNIADYNYKLPDDFTGCVIAWNSNKLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGNKPCNGVAGFNCYFPLRSYGFRPTYGVGHQPYRVVVLSFELLHAPATVCGPKKST",
    "Omicron_BQ11": "NITNLCPFDEVFNATTFASVYAWNRKRISNCVADYSVLYNFAPFFAFKCYGVSPTKLNDLCFTNVYADSFVIRGNEVSQIAPGQTGNIADYNYKLPDDFTGCVIAWNSNKLDSTVGGNYNYRYRLFRKSKLKPFERDISTEIYQAGNKPCNGVAGVNCYFPLQSYGFRPTYGVGHQPYRVVVLSFELLHAPATVCGPKKST",
    "Omicron_EG5": "NITNLCPFHEVFNATTFASVYAWNRKRISNCVADYSVIYNFAPFFAFKCYGVSPTKLNDLCFTNVYADSFVIRGNEVSQIAPGQTGNIADYNYKLPDDFTGCVIAWNSNKLDSKPSGNYNYLYRLLRKSKLKPFERDISTEIYQAGNKPCNGVAGPNCYSPLQSYGFRPTYGVGHQPYRVVVLSFELLHAPATVCGPKKST",
    "Omicron_FLip": "NITNLCPFHEVFNATTFASVYAWNRKRISNCVADYSVIYNFAPFFAFKCYGVSPTKLNDLCFTNVYADSFVIRGNEVSQIAPGQTGNIADYNYKLPDDFTGCVIAWNSNKLDSKPSGNYNYLYRFLRKSKLKPFERDISTEIYQAGNKPCNGVAGPNCYSPLQSYGFRPTYGVGHQPYRVVVLSFELLHAPATVCGPKKST",
    "Omicron_XBB15": "NITNLCPFHEVFNATTFASVYAWNRKRISNCVADYSVIYNFAPFFAFKCYGVSPTKLNDLCFTNVYADSFVIRGNEVSQIAPGQTGNIADYNYKLPDDFTGCVIAWNSNKLDSKPSGNYNYLYRLFRKSKLKPFERDISTEIYQAGNKPCNGVAGPNCYSPLQSYGFRPTYGVGHQPYRVVVLSFELLHAPATVCGPKKST",
    "Omicron_BA286": "NVTNLCPFHEVFNATRFASVYAWNRTRISNCVADYSVLYNFAPFFAFKCYGVSPTKLNDLCFTNVYADSFVIKGNEVSQIAPGQTGNIADYNYKLPDDFTGCVIAWNSNKLDSKHSGNYDYWYRLFRKSKLKPFERDISTEIYQAGNKPCKGKGPNCYFPLQSYGFRPTYGVGHQPYRVVVLSFELLHAPATVCGPKKST",
}

DEFAULT_RBD_TARGET_NAME_ALIASES = {
    "Wuhan_Hu_1": "Wuhan-Hu-1",
    "N501Y": "Alpha",
    "B1351": "Beta",
    "E484K": "Eta",
    "BA1": "Omicron_BA1",
    "BA2": "Omicron_BA2",
    "BQ11": "Omicron_BQ11",
    "EG5": "Omicron_EG5",
    "FLip": "Omicron_FLip",
    "XBB15": "Omicron_XBB15",
    "BA286": "Omicron_BA286",
}


def __dir__() -> List[str]:
    """Return exported names.

    Returns
    -------
    List[str]
        Exported names.
    """

    return __all__


@dataclass
class RBDACE2CleanerConfig(BaseCleanerConfig):
    """Configuration for the RBD ACE2 cleaner.

    Attributes
    ----------
    reference_sequences : Dict[str, str]
        Canonical RBD target reference sequences.
    target_name_aliases : Dict[str, str]
        RBD target alias-to-canonical-name mapping.
    column_mapping : Dict[str, str]
        Mapping from raw source column names to the standardized column names
        consumed by the RBD ACE2 cleaner.
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

    # Target/reference sequence configuration
    reference_sequences: Dict[str, str] = field(default_factory=lambda: deepcopy(DEFAULT_RBD_REFERENCE_SEQUENCES))

    target_name_aliases: Dict[str, str] = field(default_factory=lambda: deepcopy(DEFAULT_RBD_TARGET_NAME_ALIASES))

    # Column preparation configuration
    column_mapping: Dict[str, str] = field(
        default_factory=lambda: {
            "target": "name",
            "aa_substitutions": "mut_info",
            "log10Ka": "label",
            "variant_class": "variant_class",
            "n_aa_substitutions": "n_aa_substitutions",
        }
    )

    drop_na_columns: List[str] = field(default_factory=lambda: ["name", "label"])

    type_conversions: Dict[str, str] = field(default_factory=lambda: {"label": "float"})

    # Mutation processing configuration
    validate_mut_workers: int = 16

    process_workers: int = 16

    # Label and pipeline configuration
    label_columns: List[str] = field(default_factory=lambda: ["label"])

    primary_label_column: str = "label"

    pipeline_name: str = "RBDACE2 pipeline"

    def validate(self) -> None:
        """Validate RBD ACE2 cleaner configuration values.

        Raises
        ------
        ValueError
            If the configuration is internally inconsistent.
        """

        super().validate()

        if not self.reference_sequences:
            raise ValueError("reference_sequences cannot be empty")
        if not self.label_columns:
            raise ValueError("label_columns cannot be empty")
        if self.primary_label_column not in self.label_columns:
            raise ValueError(f"primary_label_column '{self.primary_label_column}' " f"must be in label_columns {self.label_columns}")
        required_standard_columns = {
            "name",
            "mut_info",
            "label",
            "variant_class",
        }
        missing_standard_columns = required_standard_columns - set(self.column_mapping.values())
        if missing_standard_columns:
            raise ValueError("column_mapping must provide standardized columns " f"{sorted(required_standard_columns)}, missing {sorted(missing_standard_columns)}")

        for target_name, sequence in self.reference_sequences.items():
            sequence_length = len(str(sequence).strip())
            if sequence_length <= 0:
                raise ValueError(f"Reference sequence for target '{target_name}' cannot be empty")


def create_rbd_ace2_cleaner(
    dataset_or_path: Optional[Union[pd.DataFrame, str, Path]] = None,
    config: Optional[Union[RBDACE2CleanerConfig, Dict[str, Any], str, Path]] = None,
) -> Pipeline:
    """Create the RBD ACE2 cleaning pipeline.

    Parameters
    ----------
    dataset_or_path : Optional[Union[pd.DataFrame, str, Path]], default=None
        Raw RBD ACE2 dataframe or input file path.
    config : Optional[Union[RBDACE2CleanerConfig, Dict[str, Any], str, Path]], default=None
        Cleaner configuration object, partial configuration dictionary, JSON path,
        or ``None`` to use the built-in default configuration.

    Returns
    -------
    Pipeline
        Delayed cleaning pipeline.

    Raises
    ------
    TypeError
        If ``dataset_or_path`` or ``config`` uses an unsupported type.
    """

    default_config = RBDACE2CleanerConfig()
    if config is None:
        final_config = default_config
    elif isinstance(config, RBDACE2CleanerConfig):
        final_config = config
    elif isinstance(config, dict):
        final_config = default_config.merge(config)
    elif isinstance(config, (str, Path)):
        final_config = RBDACE2CleanerConfig.from_json(config)
    else:
        raise TypeError(f"config must be RBDACE2CleanerConfig, dict, str, Path or None, got {type(config)}")

    target_name_column = final_config.column_mapping.get("target", "target")
    mutation_column = final_config.column_mapping.get("aa_substitutions", "aa_substitutions")
    variant_class_column = final_config.column_mapping.get("variant_class", "variant_class")

    logger.info(
        "RBD ACE2 dataset will be cleaned with pipeline: %s",
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
            drop_na_columns=final_config.drop_na_columns,
        )
        .delayed_then(
            standardize_rbd_target_names,
            target_name_aliases=final_config.target_name_aliases,
            name_column=target_name_column,
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
            name_columns=[target_name_column, mutation_column],
            label_columns=final_config.label_columns,
        )
        .delayed_then(
            subtract_labels_by_wt,
            name_column=target_name_column,
            label_columns=final_config.label_columns,
            mutation_column=mutation_column,
            wt_identifier="WT",
            in_place=True,
            drop_wt_row=True,
        )
        .delayed_then(
            add_reference_sequences_by_target,
            reference_sequences=final_config.reference_sequences,
            name_column=target_name_column,
            sequence_column="sequence",
        )
        .delayed_then(
            apply_mutations_to_sequences,
            sequence_column="sequence",
            name_column=target_name_column,
            mutation_column=mutation_column,
            mutation_sep=",",
            is_zero_based=True,
            sequence_type="protein",
            num_workers=final_config.process_workers,
        )
        .delayed_then(
            convert_to_mutation_dataset_format,
            name_column=target_name_column,
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


def clean_rbd_ace2_dataset(
    pipeline: Pipeline,
) -> Tuple[Pipeline, MutationDataset]:
    """Execute the RBD ACE2 cleaning pipeline.

    Parameters
    ----------
    pipeline : Pipeline
        Pipeline created by :func:`create_rbd_ace2_cleaner`.

    Returns
    -------
    Tuple[Pipeline, MutationDataset]
        Executed pipeline and cleaned mutation dataset.
    """

    pipeline.execute()
    dataset_df, reference_sequences = pipeline.data
    dataset = MutationDataset.from_dataframe(dataset_df, reference_sequences)
    logger.info(
        "Successfully cleaned RBD ACE2 dataset: %s mutations from %s references",
        len(dataset_df),
        len(reference_sequences),
    )
    return pipeline, dataset
