# mutcleaner/cleaners/antitoxin_pard3_cleaner.py
from __future__ import annotations

import pandas as pd
from typing import TYPE_CHECKING
from dataclasses import dataclass, field
from pathlib import Path
import logging

from .base_config import BaseCleanerConfig
from .basic_cleaners import (
    read_dataset,
    extract_and_rename_columns,
    filter_and_clean_data,
    convert_data_types,
    convert_to_mutation_dataset_format,
    validate_mutations,
    add_columns,
    apply_mutations_to_sequences,
    average_labels_by_name,
    subtract_labels_by_wt,
)
from .antitoxin_pard3_custom_cleaners import (
    simplify_mutations,
)

from ..core.dataset import MutationDataset
from ..core.pipeline import Pipeline, create_pipeline

if TYPE_CHECKING:
    from typing import Any, Callable, Dict, List, Optional, Tuple, Union

__all__ = [
    "AntitoxinParD3CleanerConfig",
    "create_antitoxin_pard3_cleaner",
    "clean_antitoxin_pard3_dataset",
]


def __dir__() -> List[str]:
    return __all__


# Create module logger
logger = logging.getLogger(__name__)


@dataclass
class AntitoxinParD3CleanerConfig(BaseCleanerConfig):
    """
    Configuration class for Antitoxin dataset cleaner.
    Inherits from BaseCleanerConfig and adds Antitoxin-specific configuration options.

    Simply run `mutcleaner.download_antitoxin_source_file()` to download the dataset.

    Alternatively, the raw Antitoxin file can be obtained from:

    - Hugging Face: https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/Antitoxin_ParD3_Epistasis_Dataset/Antitoxin_ParD3_Epistasis_Dataset.csv

    Attributes
    ----------
    column_mapping : Dict[str, str]
        Mapping from source to target column names
    filters : Dict[str, Callable]
        Filter conditions for data cleaning
    wt_sequence : str
        Wildtype sequence for the dataset, used for mutation validation
    type_conversions : Dict[str, str]
        Data type conversion specifications
    validate_mut_workers : int
        Number of workers for mutation validation, set to -1 to use all available CPUs
    process_workers : int
        Number of workers for applying mutations to sequences, set to -1 to use all available CPUs
    label_columns : List[str]
        List of score columns to process
    primary_label_column : str
        Primary score column for the dataset
    pipeline_name : str
        Name of the cleaning pipeline
    """

    # Column mapping configuration
    column_mapping: Dict[str, str] = field(
        default_factory=lambda: {
            "mutation": "mut_info",
            "label": "label",
        }
    )

    # Data filtering configuration
    filters: Dict[str, Callable] = field(default_factory=lambda: {"label": lambda s: pd.to_numeric(s, errors="coerce").notna()})

    # obtained from the article
    wt_sequence = "MANVEKMSVAVTPQQAAVMREAVEAGEYATASEIVREAVRDWLAKRELRHDDIRRLRQLWDEGKASGRPEPVDFDALRKEARQKLTEVPPNGR"

    # Type conversion configuration
    type_conversions: Dict[str, str] = field(default_factory=lambda: {"label": "float"})

    # Mutation validation parameters
    validate_mut_workers: int = 16

    process_workers: int = 16

    # Score columns configuration
    label_columns: List[str] = field(default_factory=lambda: ["label"])

    primary_label_column: str = "label"

    # Override default pipeline name
    pipeline_name: str = "Antitoxin Pipeline"

    def validate(self) -> None:
        """Validate Antitoxin-specific configuration parameters

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
            raise ValueError(f"primary_label_column '{self.primary_label_column}' must be in label_columns {self.label_columns}")

        # Validate column mapping
        required_mappings = {"mutation"}
        missing = required_mappings - set(self.column_mapping.keys())
        if missing:
            raise ValueError(f"Missing required column mappings: {missing}")


def create_antitoxin_pard3_cleaner(
    dataset_or_path: Optional[Union[pd.DataFrame, str, Path]] = None,
    config: Optional[Union[AntitoxinParD3CleanerConfig, Dict[str, Any], str, Path]] = None,
) -> Pipeline:
    """Create Antitoxin dataset cleaning pipeline

    Parameters
    ----------
    dataset_or_path : Optional[Union[pd.DataFrame, str, Path]], default=None
        Raw dataset DataFrame or file path to Antitoxin dataset.
    config : Optional[Union[AntitoxinCleanerConfig, Dict[str, Any], str, Path]]
        Configuration for the cleaning pipeline. Can be:
        - AntitoxinCleanerConfig object
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
        final_config = AntitoxinParD3CleanerConfig()
    elif isinstance(config, AntitoxinParD3CleanerConfig):
        final_config = config
    elif isinstance(config, dict):
        # Partial configuration - merge with defaults
        default_config = AntitoxinParD3CleanerConfig()
        final_config = default_config.merge(config)
    elif isinstance(config, (str, Path)):
        # Load from file
        final_config = AntitoxinParD3CleanerConfig.from_json(config)
    else:
        raise TypeError(f"config must be AntitoxinParD3CleanerConfig, dict, str, Path or None, got {type(config)}")

    # Log configuration summary
    logger.info(f"Antitoxin dataset will be cleaned with pipeline: {final_config.pipeline_name}")
    logger.debug(f"Configuration:\n{final_config.get_summary()}")

    try:
        # Create pipeline
        pipeline = create_pipeline(dataset_or_path, final_config.pipeline_name)

        # Add cleaning steps
        pipeline = (
            pipeline.delayed_then(
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
                add_columns,
                columns_to_add={
                    "name": "antitoxin",
                    "wt_seq": final_config.wt_sequence,
                },
            )
            .delayed_then(
                simplify_mutations,
                mutation_column=final_config.column_mapping.get("mutation", "mutation"),
                mutation_sep=":",
            )
            .delayed_then(
                validate_mutations,
                mutation_column=final_config.column_mapping.get("mutation", "mutation"),
                mutation_sep=",",
                is_zero_based=True,
                exclude_patterns="WT",
                num_workers=final_config.validate_mut_workers,
            )
            .delayed_then(
                average_labels_by_name,
                name_columns=(
                    "name",
                    final_config.column_mapping.get("mutation", "mutation"),
                ),
                label_columns=final_config.primary_label_column,
            )
            .delayed_then(
                subtract_labels_by_wt,
                name_column="name",
                label_columns=final_config.primary_label_column,
                mutation_column=final_config.column_mapping.get("mutation", "mutation"),
                wt_identifier="WT",
                in_place=True,
            )
            .delayed_then(
                apply_mutations_to_sequences,
                sequence_column="wt_seq",
                name_column="name",
                mutation_column=final_config.column_mapping.get("mutation", "mutation"),
                is_zero_based=True,
                sequence_type="protein",
                num_workers=final_config.process_workers,
            )
            .delayed_then(
                convert_to_mutation_dataset_format,
                name_column="name",
                mutation_column=final_config.column_mapping.get("mutation", "mutation"),
                sequence_column="wt_seq",
                label_column=final_config.primary_label_column,
                is_zero_based=True,
            )
        )

        # Create pipeline based on dataset_or_path type
        if isinstance(dataset_or_path, (str, Path)):
            pipeline.add_delayed_step(read_dataset, 0, file_format="csv")
        elif not isinstance(dataset_or_path, pd.DataFrame):
            raise TypeError(f"dataset_or_path must be pd.DataFrame or str/Path, got {type(dataset_or_path)}")

        return pipeline

    except Exception as e:
        logger.error(f"Error in creating aav capsid cleaning pipeline: {str(e)}")
        raise RuntimeError(f"Error in creating aav capsid cleaning pipeline: {str(e)}")


def clean_antitoxin_pard3_dataset(
    pipeline: Pipeline,
) -> Tuple[Pipeline, MutationDataset]:
    """Clean Antitoxin dataset using configurable pipeline

    Parameters
    ----------
    pipeline : Pipeline
        Antitoxin dataset cleaning pipeline

    Returns
    -------
    Tuple[Pipeline, MutationDataset]
        - Pipeline: The cleaned pipeline
        - MutationDataset: The cleaned Antitoxin dataset

    Examples
    --------
    Use default configuration:

    >>> pipeline = create_antitoxin_cleaner(df)  # df is raw Antitoxin dataset file

    Use partial configuration:

    >>> pipeline = create_antitoxin_cleaner(df, config={
    ...     "validate_mut_workers": 8,
    ... })

    Load configuration from file:

    >>> pipeline = create_antitoxin_cleaner(df, config="config.json")
    >>> pipeline, dataset = clean_antitoxin_dataset(pipeline)
    """
    try:
        # Run pipeline
        pipeline.execute()

        # Extract results
        antitoxin_dataset_df, antitoxin_ref_seq = pipeline.data
        antitoxin_dataset = MutationDataset.from_dataframe(antitoxin_dataset_df, antitoxin_ref_seq)

        logger.info(f"Successfully cleaned antitoxin dataset: {len(antitoxin_dataset_df)} mutations from {len(antitoxin_ref_seq)} proteins")

        return pipeline, antitoxin_dataset
    except Exception as e:
        logger.error(f"Error in running antitoxin dataset cleaning pipeline: {str(e)}")
        raise RuntimeError(f"Error in running antitoxin dataset cleaning pipeline: {str(e)}")
