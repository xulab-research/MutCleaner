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
    filter_stop_codon_mutations,
    convert_to_mutation_dataset_format,
    add_columns,
    average_labels_by_name,
    validate_mutations,
    apply_mutations_to_sequences,
)

from ..core.alphabet import DNAAlphabet
from ..core.dataset import MutationDataset
from ..core.mutation import CodonMutation
from ..core.pipeline import Pipeline, create_pipeline

if TYPE_CHECKING:
    from typing import Any, Callable, Dict, List, Optional, Tuple, Union

__all__ = [
    "CodonHumanMyoglobinCleanerConfig",
    "create_codon_human_myoglobin_cleaner",
    "clean_codon_human_myoglobin_dataset",
]


def __dir__() -> List[str]:
    return __all__


# Create module logger
logger = logging.getLogger(__name__)


@dataclass
class CodonHumanMyoglobinCleanerConfig(BaseCleanerConfig):
    """
    Configuration class for the Codon Human Myoglobin Epistasis Dataset cleaner.
    Inherits from BaseCleanerConfig and adds hMb-specific configuration options.

    Simply run `mutcleaner.download_codon_human_myoglobin_source_file()` to download the dataset.

    Alternatively, the raw codon hMb file can be obtained from:

    - Hugging Face: https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/Codon_Human_Myoglobin_Epistasis_Dataset/Codon_Human_Myoglobin_Epistasis_Dataset.csv

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
    process_workers : int
        Number of workers for parallel processing
    label_columns : List[str]
        List of score columns to process
    primary_label_column : str
        Primary score column for the dataset
    """

    # Column mapping configuration
    column_mapping: Dict[str, str] = field(
        default_factory=lambda: {
            "COD": "mut_info",
            "fitness": "label",
        }
    )

    # Data filtering configuration
    filters: Dict[str, Callable] = field(
        default_factory=lambda: {
            "label": lambda s: pd.to_numeric(s, errors="coerce").notna()
        }
    )

    # Type conversion configuration
    type_conversions: Dict[str, str] = field(default_factory=lambda: {"label": "float"})

    # obtained from the article
    wt_sequence = ("ATGGGTTTATCGGATGGAGAATGGCAGTTGGTACTTAATGTATGGGGTAAAGTTGAGGCAGACATCCCAGGACATGGCCA"
                   "GGAAGTATTAATCCGTTTATTTAAAGGACACCCAGAGACGCTAGAGAAATTCGATAAATTTAAGCATCTTAAATCCGAAG"
                   "ACGAGATGAAGGCTTCTGAGGACTTAAAGAAACACGGGGCTACTGTGTTGACTGCATTAGGTGGTATTCTAAAGAAGAAA"
                   "GGTCACCACGAGGCCGAAATAAAGCCACTAGCCCAGTCCCATGCTACAAAACACAAAATTCCCGTAAAATATCTAGAGTT"
                   "TATTTCAGAGTGCATAATTCAGGTTTTGCAATCTAAACACCCAGGCGACTTCGGAGCCGACGCTCAGGGTGCGATGAACA"
                   "AAGCTTTAGAATTGTTTAGGAAGGACATGGCCTCTAATTACAAGGAGCTAGGCTTCCAGGGC")

    # Mutation validation parameters
    validate_mut_workers: int = 16

    # Processing parameters
    process_workers: int = 16

    # Score columns configuration
    label_columns: List[str] = field(default_factory=lambda: ["label"])
    primary_label_column: str = "label"

    # Override default pipeline name
    pipeline_name: str = "Codon hMb Cleaning Pipeline"

    def validate(self) -> None:
        """Validate hMb-specific configuration parameters

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
        required_mappings = {"COD", "fitness"}
        missing = required_mappings - set(self.column_mapping.keys())
        if missing:
            raise ValueError(f"Missing required column mappings: {missing}")


def create_codon_human_myoglobin_cleaner(
    dataset_or_path: Optional[Union[pd.DataFrame, str, Path]] = None,
    config: Optional[
        Union[CodonHumanMyoglobinCleanerConfig, Dict[str, Any], str, Path]
    ] = None,
) -> Pipeline:
    """Create codon human myoglobin dataset cleaning pipeline

    Parameters
    ----------
    dataset_or_path : Optional[Union[pd.DataFrame, str, Path]], default=None
        Raw dataset DataFrame or file path to codon human myoglobin dataset.
    config : Optional[Union[CodonHumanMyoglobinCleanerConfig, Dict[str, Any], str, Path]]
        Configuration for the cleaning pipeline. Can be:
        - CodonHumanMyoglobinCleanerConfig object
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
        final_config = CodonHumanMyoglobinCleanerConfig()
    elif isinstance(config, CodonHumanMyoglobinCleanerConfig):
        final_config = config
    elif isinstance(config, dict):
        # Partial configuration - merge with defaults
        default_config = CodonHumanMyoglobinCleanerConfig()
        final_config = default_config.merge(config)
    elif isinstance(config, (str, Path)):
        # Load from file
        final_config = CodonHumanMyoglobinCleanerConfig.from_json(config)
    else:
        raise TypeError(
            f"config must be CodonHumanMyoglobinCleanerConfig, dict, str, Path or None, got {type(config)}"
        )

    # Log configuration summary
    logger.info(
        f"Human myoglobin dataset will be cleaned with pipeline: {final_config.pipeline_name}"
    )
    logger.debug(f"Configuration:\n{final_config.get_summary()}")

    try:
        # Create cleaning pipeline
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
                    "name": "hMb",
                    "wt_seq": final_config.wt_sequence,
                },
            )
            .delayed_then(
                validate_mutations,
                mutation_column=final_config.column_mapping.get("COD", "COD"),
                num_workers=final_config.validate_mut_workers,
                mutation_sep=",",
                mutation_type=CodonMutation,
                alphabet=DNAAlphabet(),
            )
            .delayed_then(
                filter_stop_codon_mutations,
                mutation_column=final_config.column_mapping.get("COD", "COD"),
                mutation_sep=",",
                is_zero_based=True,
                alphabet=DNAAlphabet(),
            )
            .delayed_then(
                apply_mutations_to_sequences,
                mutation_column=final_config.column_mapping.get("COD", "COD"),
                sequence_column="wt_seq",
                sequence_type="dna",
                mutation_type=CodonMutation,
                alphabet=DNAAlphabet(),
                num_workers=final_config.process_workers,
            )
            .delayed_then(
                average_labels_by_name,
                name_columns=final_config.column_mapping.get("COD", "COD"),
                label_columns=final_config.primary_label_column,
            )
            .delayed_then(
                convert_to_mutation_dataset_format,
                name_column="name",
                mutation_column=final_config.column_mapping.get("COD", "COD"),
                sequence_column="wt_seq",
                mutated_sequence_column="mut_seq",
                sequence_type="dna",
                label_column=final_config.primary_label_column,
                is_zero_based=True,
            )
        )

        # Create pipeline based on dataset_or_path type
        if isinstance(dataset_or_path, (str, Path)):
            pipeline.add_delayed_step(read_dataset, 0, file_format="csv")
        elif not isinstance(dataset_or_path, pd.DataFrame):
            raise TypeError(
                f"dataset_or_path must be pd.DataFrame or str/Path, got {type(dataset_or_path)}"
            )

        return pipeline

    except Exception as e:
        logger.error(f"Error in creating human myoglobin cleaning pipeline: {str(e)}")
        raise RuntimeError(
            f"Error in creating human myoglobin cleaning pipeline: {str(e)}"
        )


def clean_codon_human_myoglobin_dataset(
    pipeline: Pipeline,
) -> Tuple[Pipeline, MutationDataset]:
    """Clean human myoglobin dataset using configurable pipeline

    Parameters
    ----------
    pipeline : Pipeline
        Human myoglobin dataset cleaning pipeline

    Returns
    -------
    Tuple[Pipeline, MutationDataset]
        - Pipeline: The cleaned pipeline
        - MutationDataset: The cleaned human myoglobin dataset

    Examples
    --------
    Use default configuration:

    >>> pipeline = create_codon_human_myoglobin_cleaner(df)  # df is raw human myoglobin dataset file

    Use partial configuration:

    >>> pipeline = create_codon_human_myoglobin_cleaner(df, config={
    ...     "validate_mut_workers": 8,
    ... })

    Load configuration from file:

    >>> pipeline = create_codon_human_myoglobin_cleaner(df, config="config.json")
    >>> pipeline, dataset = clean_codon_human_myoglobin_dataset(pipeline)
    """
    try:
        # Run pipeline
        pipeline.execute()

        # Extract results
        human_myoglobin_dataset_df, human_myoglobin_ref_seq = pipeline.data
        human_myoglobin_dataset = MutationDataset.from_dataframe(
            human_myoglobin_dataset_df, human_myoglobin_ref_seq
        )

        logger.info(
            f"Successfully cleaned human myoglobin dataset: {len(human_myoglobin_dataset_df)} mutations from {len(human_myoglobin_ref_seq)} proteins"
        )

        return pipeline, human_myoglobin_dataset
    except Exception as e:
        logger.error(
            f"Error in running human myoglobin dataset cleaning pipeline: {str(e)}"
        )
        raise RuntimeError(
            f"Error in running human myoglobin dataset cleaning pipeline: {str(e)}"
        )
