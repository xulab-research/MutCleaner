# mutcleaner/cleaners/ddg_dtm_cleaners.py
from __future__ import annotations

import pandas as pd
from typing import TYPE_CHECKING
from dataclasses import dataclass, field
from pathlib import Path
import logging

from .base_config import BaseCleanerConfig
from .basic_cleaners import (
    read_dataset,
    split_columns,
    merge_columns,
    extract_and_rename_columns,
    infer_mutations_from_sequences,
    convert_data_types,
    aggregate_labels_by_name,
    convert_to_mutation_dataset_format,
)

from ..core.dataset import MutationDataset
from ..core.pipeline import Pipeline, create_pipeline

if TYPE_CHECKING:
    from typing import Any, Dict, List, Literal, Optional, Tuple, Union

__all__ = [
    "DdgDtmCleanerConfig",
    "create_ddg_dtm_cleaner",
    "clean_ddg_dtm_dataset",
]


def __dir__() -> List[str]:
    return __all__


# Create module logger
logger = logging.getLogger(__name__)


@dataclass
class DdgDtmCleanerConfig(BaseCleanerConfig):
    """
    Configuration class for ddG-dTm dataset cleaner.
    Inherits from BaseCleanerConfig and adds ddG-dTm-specific configuration options.

    Simply run `mutcleaner.download_ddg_dtm_source_file()` to download the dataset.

    Alternatively, the raw ddG-dTm files can be obtained from:

    - Hugging Face: https://huggingface.co/datasets/xulab-research/MutCleaner/tree/main/ddG_Dataset
    - Hugging Face: https://huggingface.co/datasets/xulab-research/MutCleaner/tree/main/dTm_Dataset

    Attributes
    ----------
    column_mapping : Dict[str, str]
        Mapping from source to target column names
    type_conversions : Dict[str, str]
        Data type conversion specifications
    infer_mut_workers : int
        Number of workers for mutation inference, set to -1 to use all available CPUs
    aggregation_strategy : Literal["mean", "first", "nearest"]
        Aggregate labels by name, see `aggregate_labels_by_name` for details
    nearest_by : List[Tuple[str, float]]
        Keep mutation by distance, see `aggregate_labels_by_name` for details
    label_columns : List[str]
        List of score columns to process
    primary_label_column : str
        Primary score column for the dataset
    """

    # Column mapping configuration
    column_mapping: Dict[str, str] = field(
        default_factory=lambda: {
            "name": "name",
            "wt_seq": "wt_seq",
            "mut_seq": "mut_seq",
            "pH": "pH",
            # 'temp' and 'label' cols are added due to dTm or ddG in `create_ddg_dtm_cleaner`
        }
    )

    # Type conversion configuration
    type_conversions: Dict[str, str] = field(default_factory=lambda: {"label": "float"})

    # Mutation inference parameters
    infer_mut_workers: int = 16

    # Score configuration
    aggregation_strategy: Literal["mean", "first", "nearest"] = "nearest"
    nearest_by: List[Tuple[str, float]] = field(default_factory=list)

    # Score columns configuration
    label_columns: List[str] = field(default_factory=lambda: ["label"])
    primary_label_column: str = "label"

    # Override default pipeline name
    pipeline_name: str = "ddG-dTm"

    def __post_init__(self):
        # If user didn't provide nearest_by, set a sensible default based on column_mapping
        if self.aggregation_strategy == "nearest" and not self.nearest_by:
            self.nearest_by = [(self.column_mapping.get("pH", "pH"), 7.0)]
        # Normalize types: ensure str,float and tuple form
        self.nearest_by = [(str(col), float(target)) for col, target in self.nearest_by]

    def validate(self) -> None:
        """Validate ddG-dTm-specific configuration parameters

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
        required_mappings = {"name", "wt_seq", "mut_seq"}
        missing = required_mappings - set(self.column_mapping.keys())
        if missing:
            raise ValueError(f"Missing required column mappings: {missing}")


def create_ddg_dtm_cleaner(
    dataset_or_path: Optional[Union[pd.DataFrame, str, Path]] = None,
    config: Optional[Union[DdgDtmCleanerConfig, Dict[str, Any], str, Path]] = None,
) -> Pipeline:
    """Create ddG-dTm dataset cleaning pipeline

    Parameters
    ----------
    dataset_or_path : Optional[Union[pd.DataFrame, str, Path]], default=None
        Raw dataset DataFrame or file path to ddG-dTm dataset.
    config : Optional[Union[DdgDtmCleanerConfig, Dict[str, Any], str, Path]]
        Configuration for the cleaning pipeline. Can be:
        - DdgDtmCleanerConfig object
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

    Notes
    -----
    Label columns (dTm or ddG) are automatically detected and added to the pipeline.

    Examples
    --------

    """
    # Handle configuration parameter
    if config is None:
        final_config = DdgDtmCleanerConfig()
    elif isinstance(config, DdgDtmCleanerConfig):
        final_config = config
    elif isinstance(config, dict):
        # Partial configuration - merge with defaults
        default_config = DdgDtmCleanerConfig()
        final_config = default_config.merge(config)
    elif isinstance(config, (str, Path)):
        # Load from file
        final_config = DdgDtmCleanerConfig.from_json(config)
    else:
        raise TypeError(f"config must be DdgDtmCleanerConfig, dict, str, Path or None, got {type(config)}")

    # Log configuration summary
    logger.info(f"ddG-dTm dataset will cleaning with pipeline: {final_config.pipeline_name}")
    logger.debug(f"Configuration:\n{final_config.get_summary()}")

    def _detect_label_columns(data: pd.DataFrame) -> str:
        colnames = data.columns
        if "dTm" in colnames or "ddG" in colnames:
            label_col = "dTm" if "dTm" in colnames else "ddG"
            return label_col
        else:
            raise ValueError("No dTm or ddG columns found in the dataset")

    try:
        # Create pipeline
        pipeline = create_pipeline(dataset_or_path, final_config.pipeline_name)

        # Detect label columns
        if isinstance(dataset_or_path, (str, Path)):
            pipeline.then(read_dataset)
            if pipeline.data is not None:
                label_col = _detect_label_columns(pipeline.data)
            else:
                raise ValueError("No data found in the dataset")
        elif isinstance(dataset_or_path, pd.DataFrame):
            label_col = _detect_label_columns(dataset_or_path)
        else:
            raise TypeError(f"dataset_or_path must be pd.DataFrame or str/Path, got {type(dataset_or_path)}")
        final_config.column_mapping.update({label_col: "label"})
        if label_col == "ddG":
            # Add temp configuration for ddG
            final_config.column_mapping.update({"temp": "temp"})
            final_config.nearest_by.append(("temp", 25))

        # Add cleaning steps
        pipeline = (
            pipeline.delayed_then(
                extract_and_rename_columns,
                column_mapping=final_config.column_mapping,
            )
            .delayed_then(
                split_columns,
                column_to_split=final_config.column_mapping.get("name", "name"),
                new_column_names=["__rcsb", "__accession", "__chain", "__other"],
                separator="_",
                max_splits=3,
                drop_original=True,
            )
            .delayed_then(
                merge_columns,
                columns_to_merge=["__rcsb", "__accession", "__chain"],
                new_column_name=final_config.column_mapping.get("name", "name"),
                separator="_",
                drop_original=True,
            )
            .delayed_then(
                infer_mutations_from_sequences,
                wt_sequence_column=final_config.column_mapping.get("wt_seq", "wt_seq"),
                mut_sequence_column=final_config.column_mapping.get("mut_seq", "mut_seq"),
                num_workers=final_config.infer_mut_workers,
            )
            .delayed_then(convert_data_types, type_conversions=final_config.type_conversions)
            .delayed_then(
                aggregate_labels_by_name,
                name_columns=[
                    final_config.column_mapping.get("name", "name"),
                    "inferred_mutations",
                ],
                label_columns=final_config.label_columns,
                remove_origin_columns=True,
                strategy=final_config.aggregation_strategy,
                nearest_by=final_config.nearest_by,
            )
            .delayed_then(
                convert_to_mutation_dataset_format,
                name_column=final_config.column_mapping.get("name", "name"),
                mutation_column="inferred_mutations",
                sequence_column=final_config.column_mapping.get("wt_seq", "wt_seq"),
                mutated_sequence_column=final_config.column_mapping.get("mut_seq", "mut_seq"),
                label_column=final_config.primary_label_column,
                is_zero_based=True,
            )
        )

        return pipeline

    except Exception as e:
        logger.error(f"Error in creating ddG-dTm cleaning pipeline: {str(e)}")
        raise RuntimeError(f"Error in creating ddG-dTm cleaning pipeline: {str(e)}")


def clean_ddg_dtm_dataset(
    pipeline: Pipeline,
) -> Tuple[Pipeline, MutationDataset]:
    """Clean ddG-dTm dataset using configurable pipeline

    Parameters
    ----------
    pipeline : Pipeline
        ddG-dTm dataset cleaning pipeline

    Returns
    -------
    Tuple[Pipeline, MutationDataset]
        - Pipeline: The cleaned pipeline
        - MutationDataset: The cleaned ddG-dTm dataset

    Examples
    --------
    >>> pipeline = create_ddg_dtm_cleaner(df)  # df is raw ddG-dTm dataset file

    Use default configuration:

    >>> pipeline, dataset = clean_ddg_dtm_dataset(pipeline)

    Use partial configuration:

    >>> pipeline, dataset = clean_ddg_dtm_dataset(df, config={
    ...     "infer_mut_workers": 8,
    ... })

    Load configuration from file:

    >>> pipeline, dataset = clean_ddg_dtm_dataset(df, config="config.json")
    """
    try:
        # Run pipeline
        pipeline.execute()

        # Extract results
        ddg_dtm_dataset_df, ddg_dtm_ref_seq = pipeline.data
        ddg_dtm_dataset = MutationDataset.from_dataframe(ddg_dtm_dataset_df, ddg_dtm_ref_seq)

        logger.info(f"Successfully cleaned ddG-dTm dataset: {len(ddg_dtm_dataset_df)} mutations from {len(ddg_dtm_ref_seq)} proteins")

        return pipeline, ddg_dtm_dataset
    except Exception as e:
        logger.error(f"Error in running ddG-dTm dataset cleaning pipeline: {str(e)}")
        raise RuntimeError(f"Error in running ddG-dTm dataset cleaning pipeline: {str(e)}")
