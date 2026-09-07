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
    add_columns,
    validate_mutations,
    average_labels_by_name,
    convert_to_mutation_dataset_format,
    apply_mutations_to_sequences,
)
from .ctxm_custom_cleaner import expand_mutations, map_ambler_positions
from ..core.dataset import MutationDataset
from ..core.pipeline import Pipeline, create_pipeline

if TYPE_CHECKING:
    from typing import Any, Callable, Dict, List, Optional, Tuple, Union

__all__ = [
    "CTXMCleanerConfig",
    "create_ctxm_cleaner",
    "clean_ctxm_dataset",
]


def __dir__() -> List[str]:
    return __all__


# Create module logger
logger = logging.getLogger(__name__)


@dataclass
class CTXMCleanerConfig(BaseCleanerConfig):
    """
    Configuration class for CTXM dataset cleaner.
    Inherits from BaseCleanerConfig and adds CTXM-specific configuration options.

    Simply run `mutcleaner.download_ctxm_source_file()` to download the dataset.

    Alternatively, the raw CTXM file can be obtained from:

    - Hugging Face: https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/CTXM_Epistasis_Dataset/CTXM_Ampicillin_Epistasis_Dataset.csv
    - Hugging Face: https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/CTXM_Epistasis_Dataset/CTXM_Cefotaxime_Epistasis_Dataset.csv

    Attributes
    ----------
    column_mapping : Dict[str, str]
        Mapping from source to target column names
    filters : Dict[str, Callable]
        Filter conditions for data cleaning
    type_conversions : Dict[str, str]
        Data type conversion specifications
    validate_mut_workers : int
        Number of workers for mutation validation, set to -1 to use all available CPUs
    validate_wt_workers : int
        Number of workers for wildtype sequence validation, set to -1 to use all available CPUs
    label_columns : List[str]
        List of score columns to process
    primary_label_column : str
        Primary score column for the dataset
    """

    # Column mapping configuration
    column_mapping: Dict[str, str] = field(
        default_factory=lambda: {
            "mut_info": "mut_info",
            "fitness": "fitness",
        }
    )

    # Data filtering configuration
    filters: Dict[str, Callable] = field(
        default_factory=lambda: {
            "is.reads0": lambda x: x == True,
        }
    )

    # Type conversion configuration
    type_conversions: Dict[str, str] = field(
        default_factory=lambda: {"fitness": "float"}
    )

    # Wildtype sequence obtained from article
    wt_sequence = "RMMFAAAACIPLLLGSAPLYAQTSAVQQKLAALEKSSGGRLGVALIDTADNTQVLYRGDERFPMCSTSKVMAAAAVLKQSETQKQLLNQPVEIKPADLVNYNPIAEKHVNGTMTLAELSAAALQYSDNTAMNKLIAQLGGPGGVTAFARAIGDETFRLDRTEPTLNTAIPGDPRDTTTPRAMAQTLRQLTLGHALGETQRAQLVTWLKGNTTGAASIRAGLPTSWTVGDKTGSGDYGTTNDIAVIWPQGRAPLVLVTYFTQPQQNAESRRDVLASAARIIAEGL"

    # process parameters
    process_workers: int = 16

    # validation parameters
    validate_mut_workers: int = 16

    # Score columns configuration
    label_columns: List[str] = field(default_factory=lambda: ["fitness"])

    ambler_to_seq_mapping: Dict[str, str] = field(
        default_factory=lambda: {
            "S70": "S66",
            "K73": "K69",
            "N104": "N100",
            "Y105": "Y101",
            "N106": "N102",
            "S130": "S126",
            "N132": "N128",
            "E166": "E162",
            "P167": "P163",
            "N170": "N166",
            "K234": "K230",
            "T235": "T231",
            "G236": "G232",
            "S237": "S233",
            "G238": "G234",
            "D240": "D235",
            "R276": "R270",
        }
    )

    # wt name
    wt_name: str = "CTXM_ampicillin"

    primary_label_column: str = "fitness"

    # Override default pipeline name
    pipeline_name: str = "CTXM Cleaning Pipeline"

    def validate(self) -> None:
        """Validate CTXM-specific configuration parameters

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
        required_mappings = {"mut_info", "fitness"}
        missing = required_mappings - set(self.column_mapping.keys())
        if missing:
            raise ValueError(f"Missing required column mappings: {missing}")


def create_ctxm_cleaner(
    dataset_or_path: Optional[Union[pd.DataFrame, str, Path]] = None,
    config: Optional[Union[CTXMCleanerConfig, Dict[str, Any], str, Path]] = None,
) -> Pipeline:
    """Create CTXM dataset cleaning pipeline

    Parameters
    ----------
    dataset_or_path : Optional[Union[pd.DataFrame, str, Path]], default=None
        Raw dataset DataFrame or file path to CTXM dataset.
    config : Optional[Union[CTXMCleanerConfig, Dict[str, Any], str, Path]]
        Configuration for the cleaning pipeline. Can be:
        - CTXMCleanerConfig object
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
        final_config = CTXMCleanerConfig()
    elif isinstance(config, CTXMCleanerConfig):
        final_config = config
    elif isinstance(config, dict):
        # Partial configuration - merge with defaults
        default_config = CTXMCleanerConfig()
        final_config = default_config.merge(config)
    elif isinstance(config, (str, Path)):
        # Load from file
        final_config = CTXMCleanerConfig.from_json(config)
    else:
        raise TypeError(
            f"config must be CTXMCleanerConfig, dict, str, Path or None, got {type(config)}"
        )

    # Log configuration summary
    logger.info(
        f"CTXM dataset will cleaning with pipeline: {final_config.pipeline_name}"
    )
    logger.debug(f"Configuration:\n{final_config.get_summary()}")

    try:
        # Create pipeline
        pipeline = create_pipeline(dataset_or_path, final_config.pipeline_name)

        # Add cleaning steps
        pipeline = (
            pipeline.delayed_then(
                read_dataset,
                sep=r"\s+",
            )
            .delayed_then(filter_and_clean_data, filters=final_config.filters)
            .delayed_then(
                expand_mutations,
            )
            .delayed_then(
                extract_and_rename_columns,
                column_mapping=final_config.column_mapping,
            )
            .delayed_then(
                convert_data_types, type_conversions=final_config.type_conversions
            )
            .delayed_then(
                add_columns,
                columns_to_add={
                    "name": final_config.wt_name,
                    "wt_seq": final_config.wt_sequence,
                },
            )
            .delayed_then(
                map_ambler_positions,
                mapping=final_config.ambler_to_seq_mapping,
                mutation_column=final_config.column_mapping.get("mut_info", "mut_info"),
            )
            .delayed_then(
                validate_mutations,
                mutation_column=final_config.column_mapping.get("mut_info", "mut_info"),
                num_workers=final_config.validate_mut_workers,
            )
            .delayed_then(
                average_labels_by_name,
                name_columns=final_config.column_mapping.get("mut_info", "mut_info"),
                label_columns=final_config.label_columns,
            )
            .delayed_then(
                apply_mutations_to_sequences,
                sequence_column="wt_seq",
                name_column="name",
                mutation_column=final_config.column_mapping.get("mut_info", "mut_info"),
                is_zero_based=True,
                num_workers=final_config.process_workers,
            )
            .delayed_then(
                convert_to_mutation_dataset_format,
                name_column="name",
                mutation_column=final_config.column_mapping.get("mut_info", "mut_info"),
                sequence_column="wt_seq",
                mutated_sequence_column="mut_seq",
                label_column=final_config.primary_label_column,
                is_zero_based=True,
            )
        )

        return pipeline

    except Exception as e:
        logger.error(f"Error in creating CTXM cleaning pipeline: {str(e)}")
        raise RuntimeError(f"Error in creating CTXM cleaning pipeline: {str(e)}")


def clean_ctxm_dataset(
    pipeline: Pipeline,
) -> Tuple[Pipeline, MutationDataset]:
    """Clean CTXM dataset using configurable pipeline

    Parameters
    ----------
    pipeline : Pipeline
        CTXM dataset cleaning pipeline

    Returns
    -------
    Tuple[Pipeline, MutationDataset]
        - Pipeline: The cleaned pipeline
        - MutationDataset: The cleaned CTXM dataset

    Examples
    --------
    >>> pipeline = create_CTXM_cleaner(df)  # df is raw CTXM dataset file
    Use default configuration:

    >>> pipeline, dataset = clean_CTXM_dataset(pipeline)

    Use partial configuration:

    >>> pipeline, dataset = clean_CTXM_dataset(df, config={
    ...     "process_workers": 8,
    ... })

    Load configuration from file:

    >>> pipeline, dataset = clean_CTXM_dataset(df, config="config.json")
    """
    try:
        # Run pipeline
        pipeline.execute()

        # Extract results
        CTXM_dataset_df, CTXM_ref_seq = pipeline.data
        CTXM_dataset = MutationDataset.from_dataframe(CTXM_dataset_df, CTXM_ref_seq)

        logger.info(
            f"Successfully cleaned CTXM dataset: {len(CTXM_dataset_df)} mutations from {len(CTXM_ref_seq)} proteins"
        )

        return pipeline, CTXM_dataset
    except Exception as e:
        logger.error(f"Error in running CTXM dataset cleaning pipeline: {str(e)}")
        raise RuntimeError(f"Error in running CTXM dataset cleaning pipeline: {str(e)}")
