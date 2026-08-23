# mutcleaner/cleaners/TrpB_cleaner.py
from __future__ import annotations

import pandas as pd
from typing import TYPE_CHECKING
from dataclasses import dataclass, field
from pathlib import Path
import logging

from .base_config import BaseCleanerConfig
from .basic_cleaners import (
    read_dataset,
    add_columns,
    extract_and_rename_columns,
    filter_and_clean_data,
    validate_mutations,
    convert_data_types,
    subtract_labels_by_wt,
    apply_mutations_to_sequences,
    convert_to_mutation_dataset_format,
    average_labels_by_name,
)
from .trpb_custom_cleaner import standardize_trpb_mutation_dataset

from ..core.dataset import MutationDataset
from ..core.pipeline import Pipeline, create_pipeline

if TYPE_CHECKING:
    from typing import Any, Dict, List, Literal, Optional, Tuple, Union, Callable

__all__ = [
    "TrpBCleanerConfig",
    "create_trpb_cleaner",
    "clean_trpb_dataset",
]


def __dir__() -> List[str]:
    return __all__


# Create module logger
logger = logging.getLogger(__name__)


@dataclass
class TrpBCleanerConfig(BaseCleanerConfig):
    """
    Configuration class for TrpB dataset cleaner.
    Inherits from BaseCleanerConfig and adds TrpB-specific configuration options.
    Simply run `mutcleaner.download_trpb_source_file()` to download the dataset.

    Alternatively, the raw TrpB file can be obtained from:
    
    - Hugging Face: https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/TrpB_Epistasis_Dataset/TrpB_Epistasis_Dataset.csv

    Attributes
    ----------
    column_mapping : Dict[str, str]
        Mapping from source to target column names
    filters : Dict[str, Any]
        Filter conditions for data cleaning
    type_conversions : Dict[str, str]
        Data type conversion specifications
    is_zero_based : bool
        Whether mutation positions are zero-based
    validation_workers : int
        Number of workers for mutation validation
    infer_wt_workers : int
        Number of workers for wildtype sequence inference
    handle_multiple_wt : Literal["error", "first", "separate"]
        Strategy for handling multiple wildtype sequences
    label_columns : List[str]
        List of score columns to process
    primary_label_column : str
        Primary score column for the dataset
    """

    # Column mapping configuration
    column_mapping: Dict[str, str] = field(
        default_factory=lambda: {
            "mutation_name": "mut_info",
            "fitness": "label",
        }
    )

    # Data filtering configuration - no specific filters needed for TrpB
    filters: Dict[str, Callable] = field(
        default_factory=lambda: {
            "label": lambda s: pd.to_numeric(s, errors="coerce").notna()
        }
    )

    # Type conversion configuration
    type_conversions: Dict[str, str] = field(default_factory=lambda: {"label": "float"})

    # Obtained from article
    wt_sequence: str = (
        "MKGYFGPYGGQYVPEILMGALEELEAAYEGIMKDESFWKEFNDLLRDYAGRPTPLYFARRLSEKYGARVYLKREDLLHTGAHKINNAIGQVLLAKLMGKTRIIAETGAGQHGVATATAAALFGMECVIYMGEEDTIRQKLNVERMKLLGAKVVPVKSGSRTLKDAIDEALRDWITNLQTTYYVFGSVVGPHPYPIIVRNFQKVIGEETKKQIPEKEGRLPDYIVACVSGGSNAAGIFYPFIDSGVKLIGVEAGGEGLETGKHAASLLKGKIGYLHGSKTFVLQDDWGQVQVSHSVSAGLDYSGVGPEHAYWRETGKVLYDAVTDEEALDAFIELSRLEGIIPALESSHALAYLKKINIKGKVVVVNLSGRGDKDLESVLNHPYVRERIRLEHHHHHH"
    )

    # workers configuration
    num_workers: int = 16

    # Score columns configuration
    label_columns: List[str] = field(default_factory=lambda: ["label"])
    primary_label_column: str = "label"

    # Override default pipeline name
    pipeline_name: str = "TrpB Cleaning Pipeline"

    def validate(self) -> None:
        """Validate TrpB-specific configuration parameters

        Raises
        ------
        ValueError
            If configuration is invalid
        """
        # Call parent validation
        super().validate()

        # Validate score columns
        if not self.label_columns:
            raise ValueError("label_columns cannot be empty")

        if self.primary_label_column not in self.label_columns:
            raise ValueError(
                f"primary_label_column '{self.primary_label_column}' must be in label_columns {self.label_columns}"
            )

        # Validate column mapping
        required_mappings = {"mutation_name", "fitness"}
        missing = required_mappings - set(self.column_mapping.keys())
        if missing:
            raise ValueError(f"Missing required column mappings: {missing}")


def create_trpb_cleaner(
    dataset_or_path: Union[str, Path],
    config: Optional[Union[TrpBCleanerConfig, Dict[str, Any], str, Path]] = None,
) -> Pipeline:
    """Create TrpB dataset cleaning pipeline

    Parameters
    ----------
    dataset_or_path : Optional[Union[pd.DataFrame, str, Path]], default=None
        Raw dataset DataFrame or file path to TrpB dataset.
    config : Optional[Union[TrpBCleanerConfig, Dict[str, Any], str, Path]]
        Configuration for the cleaning pipeline. Can be:
        - TrpBCleanerConfig object
        - Dictionary with configuration parameters (merged with defaults)
        - Path to JSON configuration file (str or Path)
        - None (uses default configuration)

    Returns
    -------
    Pipeline
        Pipeline: The cleaning pipeline used

    Raises
    ------
    TypeError
        If config has invalid type
    ValueError
        If configuration validation fails
    """

    # Handle configuration parameter
    if config is None:
        final_config = TrpBCleanerConfig()
    elif isinstance(config, TrpBCleanerConfig):
        final_config = config
    elif isinstance(config, dict):
        # Partial configuration - merge with defaults
        default_config = TrpBCleanerConfig()
        final_config = default_config.merge(config)
    elif isinstance(config, (str, Path)):
        # Load from file
        final_config = TrpBCleanerConfig.from_json(config)
    else:
        raise TypeError(
            f"config must be TrpBCleanerConfig, dict, str, Path or None, got {type(config)}"
        )

    # Log configuration summary
    logger.info(
        f"TrpB dataset will cleaning with pipeline: {final_config.pipeline_name}"
    )
    logger.debug(f"Configuration:\n{final_config.get_summary()}")

    try:
        # Create pipeline
        pipeline = create_pipeline(dataset_or_path, final_config.pipeline_name)

        # Add cleaning steps
        pipeline = (
            pipeline.delayed_then(
                standardize_trpb_mutation_dataset,
            )
            .delayed_then(
                extract_and_rename_columns,
                column_mapping=final_config.column_mapping,
            )
            .delayed_then(filter_and_clean_data, filters=final_config.filters)
            .delayed_then(
                convert_data_types, type_conversions=final_config.type_conversions
            )
            .delayed_then(
                add_columns,
                columns_to_add={"name": "TrpB", "wt_seq": final_config.wt_sequence},
            )
            .delayed_then(
                validate_mutations,
                mutation_column=final_config.column_mapping.get(
                    "mutation_name", "mutation_name"
                ),
                is_zero_based=True,
                num_workers=final_config.num_workers,
                exclude_patterns=["WT"],
            )
            .delayed_then(
                average_labels_by_name,
                name_columns=final_config.column_mapping.get(
                    "mutation_name", "mutation_name"
                ),
                label_columns=final_config.label_columns,
            )
            .delayed_then(
                subtract_labels_by_wt,
                name_column="name",
                label_columns=final_config.label_columns,
                mutation_column=final_config.column_mapping.get(
                    "mutation_name", "mutation_name"
                ),
                wt_identifier="WT",
                drop_wt_row=True,
            )
            .delayed_then(
                apply_mutations_to_sequences,
                sequence_column="wt_seq",
                mutation_column=final_config.column_mapping.get(
                    "mutation_name", "mutation_name"
                ),
                is_zero_based=True,
                num_workers=final_config.num_workers,
            )
            .delayed_then(
                convert_to_mutation_dataset_format,
                name_column="name",
                mutation_column=final_config.column_mapping.get(
                    "mutation_name", "mutation_name"
                ),
                sequence_column="wt_seq",
                mutated_sequence_column="mut_seq",
                label_column=final_config.primary_label_column,
                is_zero_based=True,
            )
        )

        # Create pipeline based on dataset_or_path type
        if isinstance(dataset_or_path, (str, Path)):
            pipeline.add_delayed_step(read_dataset, 0)
        elif not isinstance(dataset_or_path, pd.DataFrame):
            raise TypeError(
                f"dataset_or_path must be pd.DataFrame or str/Path, got {type(dataset_or_path)}"
            )

        return pipeline

    except Exception as e:
        logger.error(f"Error in creating TrpB cleaning pipeline: {str(e)}")
        raise RuntimeError(f"Error in creating TrpB cleaning pipeline: {str(e)}")


def clean_trpb_dataset(
    pipeline: Pipeline,
) -> Tuple[Pipeline, MutationDataset]:
    """Clean TrpB dataset using configurable pipeline

    Parameters
    ----------
    pipeline : Pipeline
        TrpB dataset cleaning pipeline

    Returns
    -------
    Tuple[Pipeline, MutationDataset]
        - Pipeline: The cleaned pipeline
        - MutationDataset: The cleaned TrpB dataset

    Examples
    --------
    >>> pipeline = create_TrpB_cleaner(df)  # df is raw TrpB dataset file
    Use default configuration:

    >>> pipeline, dataset = clean_TrpB_dataset(pipeline)

    Use partial configuration:

    >>> pipeline, dataset = clean_TrpB_dataset(df, config={
    ...     "validate_mut_workers": 8,
    ... })

    Load configuration from file:

    >>> pipeline, dataset = clean_TrpB_dataset(df, config="config.json")
    """
    try:
        # Run pipeline
        pipeline.execute()

        # Extract results
        TrpB_dataset_df, TrpB_ref_seq = pipeline.data
        TrpB_dataset = MutationDataset.from_dataframe(TrpB_dataset_df, TrpB_ref_seq)

        logger.info(
            f"Successfully cleaned TrpB dataset: {len(TrpB_dataset_df)} mutations from {len(TrpB_ref_seq)} proteins"
        )

        return pipeline, TrpB_dataset
    except Exception as e:
        logger.error(f"Error in running TrpB dataset cleaning pipeline: {str(e)}")
        raise RuntimeError(f"Error in running TrpB dataset cleaning pipeline: {str(e)}")
