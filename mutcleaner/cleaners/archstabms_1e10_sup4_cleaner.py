# mutcleaner/cleaners/archstabms_1e10_sup4_cleaner.py
from __future__ import annotations

import pandas as pd
from typing import TYPE_CHECKING
from dataclasses import dataclass, field
from pathlib import Path
import logging

from .basic_cleaners import (
    read_dataset,
    extract_and_rename_columns,
    filter_and_clean_data,
    convert_data_types,
    convert_to_mutation_dataset_format,
)
from .archstabms_1e10_custom_cleaners import compute_mutations

from .base_config import BaseCleanerConfig

from ..core.pipeline import Pipeline, create_pipeline
from ..core.dataset import MutationDataset

if TYPE_CHECKING:
    from typing import Callable, Optional, Tuple, Dict, Union, Any, List

__all__ = [
    "ArchStabMS1E10CleanerSup4Config",
    "create_archstabms_1e10_sup4_cleaner",
    "clean_archstabms_1e10_sup4_dataset",
]


def __dir__() -> List[str]:
    return __all__


# Create module logger
logger = logging.getLogger(__name__)


@dataclass
class ArchStabMS1E10CleanerSup4Config(BaseCleanerConfig):
    """
    Configuration class for ArchStabMS1E10 dataset cleaner.
    Inherits from BaseCleanerConfig and adds ArchStabMS1E10-specific configuration options.

    Simply run `mutcleaner.download_archstabms1e10_source_file()` to download the dataset.

    Alternatively, the raw archstabms1e10 file can be obtained from:

    - Hugging Face: https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Dataset.csv

    Attributes
    ----------
    column_mapping : Dict[str, str]
        Mapping from source to target column names
    filters : Dict[str, Callable]
        Filter conditions for ndata cleaning
    type_conversions : Dict[str, str]
        Data type conversion specifications
    label_columns : List[str]
        List of score columns to process
    primary_label_column : str
        Primary score column for the dataset
    """

    # column mapping configuration
    column_mapping: Dict[str, str] = field(
        default_factory=lambda: {
            "name": "name",
            "WT": "WT",
            "full_aa_seq": "mut_seq",
            "fitness": "fitness",
        }
    )

    # Data filtering configuration
    filters: Dict[str, Callable] = field(default_factory=lambda: {"name": lambda x: x != "1_Abundance"})

    # Type conversion configuration
    type_conversions: Dict[str, str] = field(default_factory=lambda: {"fitness": "float"})

    # Score columns configuration
    label_columns: List[str] = field(default_factory=lambda: ["fitness"])
    primary_label_column: str = "fitness"

    # Override default pipeline name
    pipeline_name: str = "archstabms1e10_cleaner"

    def validate(self) -> None:
        """Validate ArchStabMS1E10CleanerConfig

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
        required_mappings = {"name", "full_aa_seq", "fitness", "WT"}
        missing = required_mappings - set(self.column_mapping.keys())
        if missing:
            raise ValueError(f"Missing required column mappings: {missing}")


def create_archstabms_1e10_sup4_cleaner(
    dataset_or_path: Optional[Union[pd.DataFrame, str, Path]] = None,
    config: Optional[Union[ArchStabMS1E10CleanerSup4Config, Dict[str, Any], str, Path]] = None,
) -> Pipeline:
    """Create ArchStabMS1E10 dataset cleaning piipeline

    Parameters
    ----------
    dataset_or_path : Optional[Union[pd.DataFrame, str, Path]], default=None
        Raw dataset DataFrame or file path to archstams 1e10 dataset
    config : Optional[Union[ArchStabMS1E10CleanerSup4Config, Dict[str, Any], str, Path]]. default=None
        Configuration for the cleaing pipeline. Can be:
        - CDNAProtelysisCleanerConfig object
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
        final_config = ArchStabMS1E10CleanerSup4Config()
    elif isinstance(config, ArchStabMS1E10CleanerSup4Config):
        final_config = config
    elif isinstance(config, dict):
        default_config = ArchStabMS1E10CleanerSup4Config()
        final_config = default_config.merge(config)
    elif isinstance(config, (str, Path)):
        # Load from file
        final_config = ArchStabMS1E10CleanerSup4Config.from_json(config)
    else:
        raise TypeError(f"config must be ArchStabMS1E10CleanerSup4Config, dict, str, Path or None, got {type(config)}")

    # Log configuration summary
    logger.info(f"archstabms 1e10 dataset will be cleaned with pipeline: {final_config.pipeline_name}")
    logger.debug(f"Configuration:\n{final_config.get_summary()}")

    try:
        # Create pipeline
        pipeline = create_pipeline(dataset_or_path, final_config.pipeline_name)

        # Add cleaning steps
        pipeline = (
            pipeline.delayed_then(
                filter_and_clean_data,
                filters=final_config.filters,
            )
            .delayed_then(
                extract_and_rename_columns,
                column_mapping=final_config.column_mapping,
            )
            .delayed_then(
                convert_data_types,
                type_conversions=final_config.type_conversions,
            )
            .delayed_then(
                compute_mutations,
                WT_column=final_config.column_mapping.get("WT", "WT"),
                mut_seq=final_config.column_mapping.get("mut_seq", "mut_seq"),
                name_column=final_config.column_mapping.get("name", "name"),
            )
            .delayed_then(
                convert_to_mutation_dataset_format,
                mutated_sequence_column=final_config.column_mapping.get("full_aa_seq", "full_aa_seq"),
                label_column=final_config.primary_label_column,
                is_zero_based=True,
            )
        )

        # Create pipeline based on dataset_or_path type
        if isinstance(dataset_or_path, (str, Path)):
            pipeline.add_delayed_step(read_dataset, 0)
        elif not isinstance(dataset_or_path, pd.DataFrame):
            raise TypeError(f"dataset_or_path must be pd.DataFrame or str/Path, got {type(dataset_or_path)}")

        return pipeline

    except Exception as e:
        logger.error(f"Error in creating archstabms cleaning pipeline: {str(e)}")
        raise RuntimeError(f"Error in creating archstabms cleaning pipeline: {str(e)}")


def clean_archstabms_1e10_sup4_dataset(
    pipeline: Pipeline,
) -> Tuple[Pipeline, MutationDataset]:
    """Clean ArchStabMS1E10 dataset using configurable pipeline

    Parameters
    ----------
    pipeline : Pipeline
        ArchStabMS1E10 dataset cleaning pipeline

    Returns
    -------
    Tuple[Pipeline, MutationDataset]
        - Pipeline: The cleaned pipeline
        - MutationDataset: The cleaned ArchStabMS1E10 dataset

    Examples
    --------
    Use default configuration:

    >>> pipeline = create_archstabms_1e10_sup4_cleaner(df)  # df is raw ArchStabMS1E10 dataset file

    Use partial configuration:

    >>> pipeline = create_archstabms_1e10_sup4_cleaner(df, config={
    ...     "column_mapping": {
    ...         "name": "name_column",
    ...         "WT":"wt",
    ...         "full_aa_seq": "mut_seq",
    ...         "fitness": "fitness",
    ... }})

    Load configuration from file:

    >>> pipeline = create_archstabms_1e10_sup4_cleaner(df, config="config.json")
    >>> pipeline, dataset = clean_archstabms_1e10_sup4_dataset(pipeline)
    """
    try:
        # Run pipeline
        pipeline.execute()

        # Extract results
        archstabms_1e10_dataset_df, archstabms_1e10_ref_seq = pipeline.data
        archstabms_1e10_dataset = MutationDataset.from_dataframe(archstabms_1e10_dataset_df, archstabms_1e10_ref_seq)

        logger.info(f"Successfully cleaned archstabms1e10 dataset: {len(archstabms_1e10_dataset_df)} mutations from {len(archstabms_1e10_ref_seq)} proteins")

        return pipeline, archstabms_1e10_dataset
    except Exception as e:
        logger.error(f"Error in running archstabms1e10 dataset cleaning pipeline: {str(e)}")
        raise RuntimeError(f"Error in running archstabms1e10 dataset cleaning pipeline: {str(e)}")
