# mutcleaner/cleaners/cdna_proteolysis_cleaner.py
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
    validate_mutations,
    average_labels_by_name,
    add_sequences_to_dataset,
    convert_to_mutation_dataset_format,
    apply_mutations_to_sequences,
)
from ..core.alphabet import DNAAlphabet
from ..core.mutation import CodonMutation
from ..core.dataset import MutationDataset
from ..core.pipeline import Pipeline, create_pipeline

if TYPE_CHECKING:
    from typing import Any, Callable, Dict, List, Optional, Tuple, Union

__all__ = [
    "CodoncDNAProteolysisCleanerConfig",
    "create_codon_cdna_proteolysis_cleaner",
    "clean_codon_cdna_proteolysis_dataset",
]


def __dir__() -> List[str]:
    return __all__


# Create module logger
logger = logging.getLogger(__name__)


@dataclass(kw_only=True)
class CodoncDNAProteolysisCleanerConfig(BaseCleanerConfig):
    """
    Configuration class for cDNAProteolysis dataset cleaner.
    Inherits from BaseCleanerConfig and adds cDNAProteolysis-specific configuration options.

    Simply run `mutcleaner.download_cdna_proteolysis_source_file()` to download the dataset.

    Alternatively, the raw cDNAProteolysis file can be obtained from:

    - Hugging Face: https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/cDNA_Proteolysis_Dataset/Tsuboyama2023_Dataset2_Dataset3_20230416.csv

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
    # Path to sequence source or Dict with sequence data
    sequence_source: Union[Dict[str, str], str, Path]

    # Column mapping configuration
    column_mapping: Dict[str, str] = field(
        default_factory=lambda: {
            "protein": "name",
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
    type_conversions: Dict[str, str] = field(default_factory=lambda: {"label": "float64"})

    # Mutation validation parameters
    validate_mut_workers: int = 16

    # apply mutations to sequences parameters
    process_workers: int = 16

    # Score columns configuration
    label_columns: List[str] = field(default_factory=lambda: ["label"])
    primary_label_column: str = "label"

    # Override default pipeline name
    pipeline_name: str = "Codon cDNAProteolysis Pipeline"

    def validate(self) -> None:
        """Validate cDNAProteolysis-specific configuration parameters

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
        required_mappings = {"protein", "COD", "fitness"}
        missing = required_mappings - set(self.column_mapping.keys())
        if missing:
            raise ValueError(f"Missing required column mappings: {missing}")


def create_codon_cdna_proteolysis_cleaner(
    dataset_or_path: Optional[Union[pd.DataFrame, str, Path]],
    sequence_source: Union[Dict[str, str], str, Path],
    config: Optional[
        Union[CodoncDNAProteolysisCleanerConfig, Dict[str, Any], str, Path]
    ] = None,
) -> Pipeline:
    """Create cDNAProteolysis dataset cleaning pipeline

    Parameters
    ----------
    dataset_or_path : Optional[Union[pd.DataFrame, str, Path]]
        Raw dataset DataFrame or file path to cDNAProteolysis dataset.
    config : Optional[Union[CodoncDNAProteolysisCleanerConfig, Dict[str, Any], str, Path]]
        Configuration for the cleaning pipeline. Can be:
        - CodoncDNAProteolysisCleanerConfig object
        - Dictionary with configuration parameters (merged with defaults)
        - Path to JSON configuration file (str or Path)
        - None (uses default configuration)

    Returns
    -------
    Pipeline
        The cleaning pipeline used

    Raises
    ------
    TypeError
        If config has invalid type
    ValueError
        If configuration validation fails
    """
    seq_path_obj = Path(sequence_source)
    if not seq_path_obj.exists():
        raise FileNotFoundError(
            f"Sequence dictionary file does not exist: {sequence_source}"
        )

    # Handle configuration parameter
    if config is None:
        final_config = CodoncDNAProteolysisCleanerConfig(
            sequence_source=sequence_source
        )
    elif isinstance(config, CodoncDNAProteolysisCleanerConfig):
        final_config = config
        # Override sequence_source if not set
        if final_config.sequence_source is None:
            final_config.sequence_source = sequence_source
    elif isinstance(config, dict):
        # Partial configuration - merge with defaults
        default_config = CodoncDNAProteolysisCleanerConfig(
            sequence_source=sequence_source
        )
        final_config = default_config.merge(config)
    elif isinstance(config, (str, Path)):
        # Load from file
        final_config = CodoncDNAProteolysisCleanerConfig.from_json(config)
        # Override sequence_source if not set
        if final_config.sequence_source is None:
            final_config.sequence_source = sequence_source
    else:
        raise TypeError(
            f"config must be CodoncDNAProteolysisCleanerConfig, dict, str, Path or None, got {type(config)}"
        )

    # Log configuration summary
    logger.info(
        f"Codon cDNAProteolysis dataset will clean with pipeline: {final_config.pipeline_name}"
    )
    logger.debug(f"Configuration:\n{final_config.get_summary()}")

    mutation_column=final_config.column_mapping.get("COD", "COD")
    name_column=final_config.column_mapping.get("protein", "protein")

    try:
        # Create pipeline
        pipeline = create_pipeline(dataset_or_path, final_config.pipeline_name)

        # Add cleaning steps
        pipeline = (
            pipeline.delayed_then(
                extract_and_rename_columns,
                column_mapping=final_config.column_mapping,
            )
            .delayed_then(filter_and_clean_data, filters=final_config.filters)
            .delayed_then(
                convert_data_types, type_conversions=final_config.type_conversions
            )
            .delayed_then(
                validate_mutations,
                mutation_column=mutation_column,
                mutation_sep=",",
                is_zero_based=True,
                mutation_type=CodonMutation,
                alphabet=DNAAlphabet(),
                num_workers=final_config.validate_mut_workers,
            )
            .delayed_then(
                add_sequences_to_dataset,
                name_column=name_column,
                sequence_source=final_config.sequence_source,
            )
            .delayed_then(
                average_labels_by_name,
                name_columns=(
                    mutation_column,
                    name_column,
                ),
                label_columns=final_config.label_columns,
            )
            .delayed_then(
                apply_mutations_to_sequences,
                mutation_column=mutation_column,
                sequence_column="sequence",
                sequence_type="dna",
                mutation_type=CodonMutation,
                alphabet=DNAAlphabet(),
                num_workers=final_config.process_workers,
            )
            .delayed_then(
                convert_to_mutation_dataset_format,
                name_column=name_column,
                mutation_column=mutation_column,
                sequence_column="sequence",
                mutated_sequence_column="mut_seq",
                sequence_type="dna",
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
        logger.error(f"Error in creating cDNAProteolysis cleaning pipeline: {str(e)}")
        raise RuntimeError(
            f"Error in creating cDNAProteolysis cleaning pipeline: {str(e)}"
        )


def clean_codon_cdna_proteolysis_dataset(
    pipeline: Pipeline,
) -> Tuple[Pipeline, MutationDataset]:
    """Clean cDNAProteolysis dataset using configurable pipeline

    Parameters
    ----------
    pipeline : Pipeline
        cDNAProteolysis dataset cleaning pipeline

    Returns
    -------
    Tuple[Pipeline, MutationDataset]
        - Pipeline: The cleaned pipeline
        - MutationDataset: The cleaned cDNAProteolysis dataset

    Examples
    --------
    Use default configuration:

    >>> pipeline = create_codon_cdna_proteolysis_cleaner(df)  # df is raw cDNAProteolysis dataset file

    Use partial configuration:

    >>> pipeline = create_codon_cdna_proteolysis_cleaner(df, config={
    ...     "validate_mut_workers": 8,
    ... })

    Load configuration from file:

    >>> pipeline = create_codon_cdna_proteolysis_cleaner(df, config="config.json")
    >>> pipeline, dataset = clean_codon_cdna_proteolysis_dataset(pipeline)
    """
    try:
        # Run pipeline
        pipeline.execute()

        # Extract results
        cdna_proteolysis_dataset_df, cdna_proteolysis_ref_seq = pipeline.data
        cdna_proteolysis_dataset = MutationDataset.from_dataframe(
            cdna_proteolysis_dataset_df, cdna_proteolysis_ref_seq
        )

        logger.info(
            f"Successfully cleaned cDNAProteolysis dataset:{len(cdna_proteolysis_dataset_df)} mutations from {len(cdna_proteolysis_ref_seq)} proteins"
        )

        return pipeline, cdna_proteolysis_dataset
    except Exception as e:
        logger.error(
            f"Error in running cDNAProteolysis dataset cleaning pipeline: {str(e)}"
        )
        raise RuntimeError(
            f"Error in running cDNAProteolysis dataset cleaning pipeline: {str(e)}"
        )
