# mutcleaner/cleaners/human_domainome_sup2_cleaner.py
from __future__ import annotations

import logging
import pandas as pd
from typing import TYPE_CHECKING
from dataclasses import dataclass, field
from pathlib import Path

from .base_config import BaseCleanerConfig
from .basic_cleaners import (
    read_dataset,
    extract_and_rename_columns,
    filter_and_clean_data,
    convert_data_types,
    convert_to_mutation_dataset_format,
    validate_mutations,
    infer_wildtype_sequences,
)
from .human_domainome_custom_cleaners import generate_mutation_strings
from ..core.dataset import MutationDataset
from ..core.pipeline import Pipeline, create_pipeline

if TYPE_CHECKING:
    from typing import Any, Callable, Dict, List, Literal, Optional, Tuple, Union

__all__ = [
    "HumanDomainomeSup2CleanerConfig",
    "create_human_domainome_sup2_cleaner",
    "clean_human_domainome_sup2_dataset",
]


def __dir__() -> List[str]:
    return __all__


# Create module logger
logger = logging.getLogger(__name__)


@dataclass(kw_only=True)
class HumanDomainomeSup2CleanerConfig(BaseCleanerConfig):
    """
    Configuration class for HumanDomainome dataset cleaner - SupplementaryTable2.
    Inherits from BaseCleanerConfig and adds HumanDomainome-specific configuration options.

    Simply run `mutcleaner.download_human_domainome_source_file()` to download the dataset.

    Alternatively, the raw HumanDomainome file and the wild type fasta file can be obtained from:

    - Hugging Face: https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/Human_Domainome_Dataset/SupplementaryTable2.txt
    - Hugging Face: https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/Human_Domainome_Dataset/wild_type.fasta

    Attributes
    ----------
    column_mapping : Dict[str, str]
        Mapping from source to target column names
    filters : Dict[str, Callable]
        Filter conditions for data cleaning
    type_conversions : Dict[str, str]
        Data type conversion specifications
    drop_na_columns: List[str]
        List of column names where null values should be dropped
    validation_workers : int
        Number of workers for mutations validation, set to -1 to use all available CPUs
    infer_wt_workers : int
        Number of workers for wildtype sequences inference, set to -1 to use all available CPUs
    label_columns : List[str]
        List of score columns to process
    primary_label_column : str
        Primary score column for the dataset
    """

    # Column mapping configuration
    column_mapping: Dict[str, str] = field(
        default_factory=lambda: {
            "domain_ID": "name",
            "aa_seq": "mut_seq",
            "wt_aa": "wt_aa",
            "mut_aa": "mut_aa",
            "position": "pos",
            "normalized_fitness": "label_humanDomainome",
        }
    )

    # Exclude nonsense mutations by default
    filters: Dict[str, Callable] = field(default_factory=lambda: {"mut_aa": lambda x: x != "*"})

    # columns to perfrom dropping NA
    drop_na_columns: List = field(
        default_factory=lambda: [
            "name",
            "mut_seq",
            "wt_aa",
            "mut_aa",
            "pos",
            "label_humanDomainome",
        ]
    )

    # Type conversion configuration
    type_conversions: Dict[str, str] = field(default_factory=lambda: {"label_humanDomainome": "float"})

    # Mutation validation parameters
    validation_workers: int = 16

    # Wildtype inference parameters
    infer_wt_workers: int = 16
    handle_multiple_wt: Literal["error", "first", "separate"] = "error"

    # Score columns configuration
    label_columns: List[str] = field(default_factory=lambda: ["label_humanDomainome"])
    primary_label_column: str = "label_humanDomainome"

    # Override default pipeline name
    pipeline_name: str = "human_domainome_cleaner"

    def __post_init__(self):
        self.type_conversions.update({"pos": "int", "mut_rel_pos": "int"})
        return super().__post_init__()

    def validate(self) -> None:
        """Validate HumanDomainome-specific configuration parameters

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
        required_mappings = set(self.column_mapping.keys())
        missing = required_mappings - set(self.column_mapping.keys())
        if missing:
            raise ValueError(f"Missing required column mappings: {missing}")


def create_human_domainome_sup2_cleaner(
    dataset_or_path: Union[str, Path, pd.DataFrame],
    config: Optional[Union[HumanDomainomeSup2CleanerConfig, Dict[str, Any], str, Path]] = None,
) -> Pipeline:
    """Create HumanDomainome ledataset cleaning pipeline - SupplementaryTable2

    Parameters
    ----------
    dataset_or_path : Union[pd.DataFrame, str, Path]
        Raw HumanDomainome dataset DataFrame or file path to HumanDomainome
        - File: `SupplementaryTable2.txt` from the article
          'Site-saturation mutagenesis of 500 human protein domains'
    config : Optional[Union[HumanDomainomeCleanerConfig, Dict[str, Any], str, Path]]
        Configuration for the cleaning pipeline. Can be:
        - HumanDomainomeCleanerConfig object
        - Dictionary with configuration parameters (merged with defaults)
        - Path to JSON configuration file (str or Path)
        - None (uses default configuration)

    Returns
    -------
    Pipeline
        The cleaning pipeline

    Raises
    ------
    FileNotFoundError
        If data file or sequence dictionary file not found
    TypeError
        If config has invalid type
    ValueError
        If configuration validation fails

    Examples
    --------
    Basic usage:

    >>> pipeline = create_human_domainome_sup2_cleaner(
    ...     "human_domainome.csv"
    ... )
    >>> pipeline, dataset = clean_human_domainome_dataset(pipeline)

    Custom configuration:

    >>> config = {
    ...     "process_workers": 8,
    ...     "type_conversions": {"label_humanDomainome": "float32"}
    ... }
    >>> pipeline = create_human_domainome_sup2_cleaner(
    ...     "human_domainome.csv"
    ...     config=config
    ... )

    Load configuration from file:

    >>> pipeline = create_human_domainome_sup2_cleaner(
    ...     "data.csv",
    ...     config="config.json"
    ... )
    """
    # Handle configuration parameter
    if config is None:
        final_config = HumanDomainomeSup2CleanerConfig()
    elif isinstance(config, HumanDomainomeSup2CleanerConfig):
        final_config = config
    elif isinstance(config, dict):
        # Partial configuration - merge with defaults
        default_config = HumanDomainomeSup2CleanerConfig()
        final_config = default_config.merge(config)
    elif isinstance(config, (str, Path)):
        # Load from file
        final_config = HumanDomainomeSup2CleanerConfig.from_json(config)
    else:
        raise TypeError(f"config must be HumanDomainomeSup2CleanerConfig, dict, str, Path or None, got {type(config)}")

    # Log configuration summary
    logger.info(f"HumanDomainome dataset (SupplementaryTable2) will be cleaned with pipeline: {final_config.pipeline_name}")
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
                drop_na_columns=final_config.drop_na_columns,
            )
            .delayed_then(
                convert_data_types,
                type_conversions=final_config.type_conversions,
            )
            .delayed_then(
                generate_mutation_strings,
                name_column=final_config.column_mapping.get("domain_ID", "domain_ID"),
                wt_aa_column=final_config.column_mapping.get("wt_aa", "wt_aa"),
                mut_aa_column=final_config.column_mapping.get("mut_aa", "mut_aa"),
                aa_pos_column=final_config.column_mapping.get("position", "pos"),
            )
            .delayed_then(
                validate_mutations,
                mutation_column="mut_info",
                format_mutations=False,
                is_zero_based=True,
                num_workers=final_config.validation_workers,
            )
            .delayed_then(
                infer_wildtype_sequences,
                name_column=final_config.column_mapping.get("domain_ID", "domain_ID"),
                mutation_column="mut_info",
                sequence_column=final_config.column_mapping.get("aa_seq", "aa_seq"),
                label_columns=final_config.label_columns,
                handle_multiple_wt=final_config.handle_multiple_wt,
                is_zero_based=True,
                num_workers=final_config.infer_wt_workers,
            )
            .delayed_then(
                convert_to_mutation_dataset_format,
                name_column=final_config.column_mapping.get("domain_ID", "domain_ID"),
                mutation_column="mut_info",
                mutated_sequence_column=final_config.column_mapping.get("aa_seq", "aa_seq"),
                label_column=final_config.primary_label_column,
                is_zero_based=True,
            )
        )

        # Create pipeline based on dataset_or_path type
        if isinstance(dataset_or_path, (str, Path)):
            pipeline.add_delayed_step(read_dataset, 0, file_format="tsv")
        elif not isinstance(dataset_or_path, pd.DataFrame):
            raise TypeError(f"dataset_or_path must be pd.DataFrame or str/Path, got {type(dataset_or_path)}")

        return pipeline

    except Exception as e:
        logger.error(f"Error in creating HumanDomainome cleaning pipeline: {str(e)}")
        raise RuntimeError(f"Error in creating HumanDomainome cleaning pipeline: {str(e)}")


def clean_human_domainome_sup2_dataset(
    pipeline: Pipeline,
) -> Tuple[Pipeline, MutationDataset]:
    """Clean HumanDomainome dataset using configurable pipeline - SupplementaryTable2

    Parameters
    ----------
    pipeline : Pipeline
        HumanDomainome dataset cleaning pipeline

    Returns
    -------
    Tuple[Pipeline, MutationDataset]
        - Pipeline: The cleaned pipeline
        - MutationDataset: The cleaned HumanDomainome dataset

    Raises
    ------
    RuntimeError
        If pipeline execution fails
    """
    try:
        # Run pipeline
        pipeline.execute()

        # Extract results
        dataset_df, ref_sequences = pipeline.data
        human_domainome_dataset = MutationDataset.from_dataframe(dataset_df, ref_sequences)

        logger.info(f"Successfully cleaned HumanDomainome dataset: {len(dataset_df)} mutations from {len(ref_sequences)} proteins")

        return pipeline, human_domainome_dataset

    except Exception as e:
        logger.error(f"Error in running HumanDomainome cleaning pipeline: {str(e)}")
        raise RuntimeError(f"Error in running HumanDomainome cleaning pipeline: {str(e)}")
