from __future__ import annotations

import logging
import pandas as pd
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from .base_config import BaseCleanerConfig
from .basic_cleaners import (
    average_labels_by_name,
    apply_mutations_to_sequences,
    convert_data_types,
    convert_to_mutation_dataset_format,
    extract_and_rename_columns,
    filter_and_clean_data,
    validate_mutations,
)
from .codon_dms_substitutions_custom_cleaners import (
    filiter_stop_codon_mutation,
    read_codon_dms_substitutions_dataset,
)
from ..core.alphabet import DNAAlphabet
from ..core.dataset import MutationDataset
from ..core.mutation import CodonMutation
from ..core.pipeline import Pipeline, create_pipeline

if TYPE_CHECKING:
    from typing import Any, Callable, Dict, List, Optional, Tuple, Union

__all__ = [
    "CodonDMSSubstitutionsCleanerConfig",
    "create_codon_dms_substitutions_cleaner",
    "clean_codon_dms_substitutions_dataset",
]


def __dir__() -> List[str]:
    return __all__


logger = logging.getLogger(__name__)


@dataclass
class CodonDMSSubstitutionsCleanerConfig(BaseCleanerConfig):
    """
    Configuration class for the Codon DMS Substitutions Dataset cleaner.

    Parameters
    ----------
    column_mapping : Dict[str, str]
        Mapping from source column names to standardized column names.
    filters : Dict[str, Callable]
        Filtering rules applied during data cleaning.
    type_conversions : Dict[str, str]
        Data type conversion specifications.
    mutation_sep : str
        Separator used between multiple codon mutations.
    is_zero_based : bool
        Whether mutation positions in the source dataset are zero-based.
    validate_mut_workers : int
        Number of workers used for mutation validation.
    process_workers : int
        Number of workers used for mutation application.
    label_columns : List[str]
        Label columns to aggregate.
    primary_label_column : str
        Primary label column used in the final dataset.
    """

    column_mapping: Dict[str, str] = field(
        default_factory=lambda: {
            "name": "name",
            "codon_mutation": "mut_info",
            "label": "label",
            "wt_sequence": "wt_seq",
        }
    )

    filters: Dict[str, Callable] = field(
        default_factory=lambda: {
            "label": lambda s: pd.to_numeric(s, errors="coerce").notna()
        }
    )

    type_conversions: Dict[str, str] = field(
        default_factory=lambda: {"label": "float64"}
    )

    mutation_sep: str = ","
    is_zero_based: bool = True

    validate_mut_workers: int = 16
    process_workers: int = 16

    label_columns: List[str] = field(default_factory=lambda: ["label"])
    primary_label_column: str = "label"

    pipeline_name: str = "Codon DMS Substitutions Pipeline"

    def validate(self) -> None:
        """Validate Codon DMS Substitutions cleaner configuration."""
        super().validate()

        if not self.label_columns:
            raise ValueError("label_columns cannot be empty")

        if self.primary_label_column not in self.label_columns:
            raise ValueError(
                f"primary_label_column '{self.primary_label_column}' "
                f"must be in label_columns {self.label_columns}"
            )

        required_columns = {"name", "codon_mutation", "label", "wt_sequence"}
        missing = required_columns - set(self.column_mapping)
        if missing:
            raise ValueError(f"Missing required column mappings: {missing}")


def create_codon_dms_substitutions_cleaner(
    data_path: Union[str, Path],
    config: Optional[
        Union[CodonDMSSubstitutionsCleanerConfig, Dict[str, Any], str, Path]
    ] = None,
) -> Pipeline:
    """
    Create the Codon DMS Substitutions Dataset cleaning pipeline.

    Parameters
    ----------
    data_path : Union[str, Path]
        Path to the dataset directory or ZIP archive.
    config : CodonDMSSubstitutionsCleanerConfig, Dict[str, Any], str, Path, optional
        Cleaner configuration, partial configuration dictionary, JSON
        configuration path, or None to use defaults.

    Returns
    -------
    Pipeline
        Configured cleaning pipeline.

    Raises
    ------
    FileNotFoundError
        If ``data_path`` does not exist.
    TypeError
        If ``data_path`` is neither a directory nor a ZIP file, or ``config``
        has an invalid type.
    """
    data_path = Path(data_path)

    if not data_path.exists():
        raise FileNotFoundError(f"Data path does not exist: {data_path}")
    if not data_path.is_dir() and not (
        data_path.is_file() and data_path.suffix.lower() == ".zip"
    ):
        raise TypeError(f"Data path must be a directory or ZIP file: {data_path}")

    if config is None:
        final_config = CodonDMSSubstitutionsCleanerConfig()
    elif isinstance(config, CodonDMSSubstitutionsCleanerConfig):
        final_config = config
    elif isinstance(config, dict):
        final_config = CodonDMSSubstitutionsCleanerConfig().merge(config)
    elif isinstance(config, (str, Path)):
        final_config = CodonDMSSubstitutionsCleanerConfig.from_json(config)
    else:
        raise TypeError(
            "config must be CodonDMSSubstitutionsCleanerConfig, "
            f"dict, str, Path or None, got {type(config)}"
        )

    logger.info(
        f"Codon DMS Substitutions Dataset will be cleaned with pipeline: "
        f"{final_config.pipeline_name}"
    )
    logger.debug(f"Configuration:\n{final_config.get_summary()}")

    try:
        pipeline = create_pipeline(data_path, final_config.pipeline_name)

        pipeline = (
            pipeline.delayed_then(read_codon_dms_substitutions_dataset)
            .delayed_then(
                extract_and_rename_columns,
                column_mapping=final_config.column_mapping,
            )
            .delayed_then(
                filter_and_clean_data,
                filters=final_config.filters,
            )
            .delayed_then(
                convert_data_types,
                type_conversions=final_config.type_conversions,
            )
            .delayed_then(
                validate_mutations,
                mutation_column="mut_info",
                mutation_sep=final_config.mutation_sep,
                is_zero_based=final_config.is_zero_based,
                mutation_type=CodonMutation,
                alphabet=DNAAlphabet(),
                num_workers=final_config.validate_mut_workers,
            )
            .delayed_then(
                filiter_stop_codon_mutation,
                mutation_column="mut_info",
                mutation_sep=",",
                is_zero_based=True,
                alphabet=DNAAlphabet(),
            )
            .delayed_then(
                average_labels_by_name,
                name_columns=("name", "mut_info"),
                label_columns=final_config.label_columns,
            )
            .delayed_then(
                apply_mutations_to_sequences,
                mutation_column="mut_info",
                sequence_column="wt_seq",
                sequence_type="dna",
                mutation_type=CodonMutation,
                alphabet=DNAAlphabet(),
                num_workers=final_config.process_workers,
            )
            .delayed_then(
                convert_to_mutation_dataset_format,
                name_column="name",
                mutation_column="mut_info",
                sequence_column="wt_seq",
                mutated_sequence_column="mut_seq",
                sequence_type="dna",
                label_column=final_config.primary_label_column,
                is_zero_based=True,
            )
        )

        return pipeline

    except Exception as e:
        logger.error(f"Error in creating Codon DMS Substitutions cleaning pipeline: {e}")
        raise RuntimeError(
            f"Error in creating Codon DMS Substitutions cleaning pipeline: {e}"
        ) from e


def clean_codon_dms_substitutions_dataset(
    pipeline: Pipeline,
) -> Tuple[Pipeline, MutationDataset]:
    """
    Clean the Codon DMS Substitutions Dataset.

    Parameters
    ----------
    pipeline : Pipeline
        Configured Codon DMS Substitutions cleaning pipeline.

    Returns
    -------
    Tuple[Pipeline, MutationDataset]
        Executed pipeline and cleaned mutation dataset.
    """
    try:
        pipeline.execute()

        dataset_df, reference_sequences = pipeline.data
        dataset = MutationDataset.from_dataframe(dataset_df, reference_sequences)

        logger.info(
            f"Successfully cleaned Codon DMS Substitutions Dataset: "
            f"{len(dataset_df)} mutations from {len(reference_sequences)} sequences"
        )

        return pipeline, dataset

    except Exception as e:
        logger.error(
            f"Error in running Codon DMS Substitutions cleaning pipeline: {e}"
        )
        raise RuntimeError(
            f"Error in running Codon DMS Substitutions cleaning pipeline: {e}"
        ) from e