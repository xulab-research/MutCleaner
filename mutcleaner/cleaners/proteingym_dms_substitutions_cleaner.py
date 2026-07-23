# mutcleaner/cleaners/proteingym_dms_substitutions_cleaner.py
from __future__ import annotations

from typing import TYPE_CHECKING
from dataclasses import dataclass, field
from pathlib import Path
import logging

from .base_config import BaseCleanerConfig
from .basic_cleaners import (
    extract_and_rename_columns,
    filter_and_clean_data,
    convert_data_types,
    validate_mutations,
    infer_wildtype_sequences,
    convert_to_mutation_dataset_format,
)
from .proteingym_dms_substitutions_custom_cleaners import read_proteingym_dms_substitutions_data
from ..core.dataset import MutationDataset
from ..core.pipeline import Pipeline, create_pipeline

if TYPE_CHECKING:
    from typing import Any, Dict, List, Literal, Optional, Tuple, Union

__all__ = [
    "ProteinGymCleanerConfig",
    "create_proteingym_dms_substitutions_cleaner",
    "clean_proteingym_dms_substitutions_dataset",
]


def __dir__() -> List[str]:
    return __all__


# Create module logger
logger = logging.getLogger(__name__)


@dataclass
class ProteinGymCleanerConfig(BaseCleanerConfig):
    """
    Configuration class for ProteinGym dataset cleaner.
    Inherits from BaseCleanerConfig and adds ProteinGym-specific configuration options.

    Simply run `mutcleaner.download_protein_gym_source_file()` to download the dataset.

    Alternatively, the raw ProteinGym file can be obtained from:

    - Hugging Face: https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/ProteinGym_DMS_Substitutions_Dataset/DMS_ProteinGym_substitutions.zip

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
            "name": "name",
            "mutated_sequence": "mut_seq",
            "mutant": "mut_info",
            "DMS_score": "DMS_score",
        }
    )

    # Data filtering configuration - no specific filters needed for ProteinGym
    filters: Dict[str, Any] = field(default_factory=dict)

    # Type conversion configuration
    type_conversions: Dict[str, str] = field(
        default_factory=lambda: {"DMS_score": "float"}
    )

    # Mutation validation parameters
    validation_workers: int = 16

    # Wildtype inference parameters
    infer_wt_workers: int = 16
    handle_multiple_wt: Literal["error", "first", "separate"] = (
        "error"  # ProteinGym default
    )

    # Score columns configuration
    label_columns: List[str] = field(default_factory=lambda: ["DMS_score"])
    primary_label_column: str = "DMS_score"

    # Override default pipeline name
    pipeline_name: str = "ProteinGym Pipeline"

    def validate(self) -> None:
        """Validate ProteinGym-specific configuration parameters

        Raises
        ------
        ValueError
            If configuration is invalid
        """
        # Call parent validation
        super().validate()

        # Validate handle_multiple_wt parameter
        valid_strategies = ["error", "first", "separate"]
        if self.handle_multiple_wt not in valid_strategies:
            raise ValueError(
                f"handle_multiple_wt must be one of {valid_strategies}, got {self.handle_multiple_wt}"
            )

        # Validate score columns
        if not self.label_columns:
            raise ValueError("label_columns cannot be empty")

        if self.primary_label_column not in self.label_columns:
            raise ValueError(
                f"primary_label_column '{self.primary_label_column}' must be in label_columns {self.label_columns}"
            )

        # Validate column mapping
        required_mappings = {"name", "mutated_sequence", "mutant", "DMS_score"}
        missing = required_mappings - set(self.column_mapping.keys())
        if missing:
            raise ValueError(f"Missing required column mappings: {missing}")


def create_proteingym_dms_substitutions_cleaner(
    data_path: Union[str, Path],
    config: Optional[Union[ProteinGymCleanerConfig, Dict[str, Any], str, Path]] = None,
) -> Pipeline:
    """Create ProteinGym dataset cleaning pipeline

    Parameters
    ----------
    data_path : Union[str, Path]
        Path to directory containing ProteinGym CSV files or path to zip file
        - Download from: https://proteingym.org/download
        - File: DMS_ProteinGym_substitutions.zip
    config : Optional[Union[ProteinGymCleanerConfig, Dict[str, Any], str, Path]]
        Configuration for the cleaning pipeline. Can be:
        - ProteinGymCleanerConfig object
        - Dictionary with configuration parameters (merged with defaults)
        - Path to JSON configuration file (str or Path)
        - None (uses default configuration)

    Returns
    -------
    Pipeline
        The cleaning pipeline

    Raises
    ------
    TypeError
        If config has invalid type
    ValueError
        If configuration validation fails

    Examples
    --------
    Process directory of ProteinGym CSV files:

    >>> pipeline = create_protein_gym_cleaner("DMS_ProteinGym_substitutions/")
    >>> pipeline, dataset = clean_protein_gym_dataset(pipeline)

    Process zip file:

    >>> pipeline = create_protein_gym_cleaner("DMS_ProteinGym_substitutions.zip")
    >>> pipeline, dataset = clean_protein_gym_dataset(pipeline)

    Custom configuration:

    >>> config = {
    ...     "validation_workers": 8,
    ...     "handle_multiple_wt": "first"
    ... }
    >>> pipeline = create_protein_gym_cleaner("data/", config=config)

    Load configuration from file:

    >>> pipeline = create_protein_gym_cleaner("data/", config="config.json")
    """
    # Validate input path
    path_obj = Path(data_path)
    if not path_obj.exists():
        raise FileNotFoundError(f"Data path does not exist: {data_path}")

    # ProteinGym only supports directory or zip file input
    if not (path_obj.is_dir() or path_obj.suffix.lower() == ".zip"):
        raise TypeError(
            f"ProteinGym cleaner only supports directory or zip file input, got: {data_path}"
        )

    # Handle configuration parameter
    if config is None:
        final_config = ProteinGymCleanerConfig()
    elif isinstance(config, ProteinGymCleanerConfig):
        final_config = config
    elif isinstance(config, dict):
        # Partial configuration - merge with defaults
        default_config = ProteinGymCleanerConfig()
        final_config = default_config.merge(config)
    elif isinstance(config, (str, Path)):
        # Load from file
        final_config = ProteinGymCleanerConfig.from_json(config)
    else:
        raise TypeError(
            f"config must be ProteinGymCleanerConfig, dict, str, Path or None, got {type(config)}"
        )

    # Log configuration summary
    logger.info(
        f"ProteinGym dataset will be cleaned with pipeline: {final_config.pipeline_name}"
    )
    logger.debug(f"Configuration:\n{final_config.get_summary()}")

    try:
        # Create pipeline
        pipeline = create_pipeline(data_path, final_config.pipeline_name)

        # Add cleaning steps using basic_cleaners functions
        pipeline = (
            pipeline.delayed_then(
                read_proteingym_dms_substitutions_data,
            )
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
                mutation_column=final_config.column_mapping.get("mutant", "mutant"),
                mutation_sep=":",
                is_zero_based=False,
                num_workers=final_config.validation_workers,
            )
            .delayed_then(
                infer_wildtype_sequences,
                name_column=final_config.column_mapping.get("name", "name"),
                mutation_column=final_config.column_mapping.get("mutant", "mutant"),
                sequence_column=final_config.column_mapping.get(
                    "mutated_sequence", "mutated_sequence"
                ),
                label_columns=final_config.label_columns,
                mutation_sep=",",
                is_zero_based=True,  # Always True after validate_mutations
                handle_multiple_wt=final_config.handle_multiple_wt,
                num_workers=final_config.infer_wt_workers,
            )
            .delayed_then(
                convert_to_mutation_dataset_format,
                name_column=final_config.column_mapping.get("name", "name"),
                mutation_column=final_config.column_mapping.get("mutant", "mutant"),
                mutated_sequence_column=final_config.column_mapping.get(
                    "mutated_sequence", "mutated_sequence"
                ),
                label_column=final_config.primary_label_column,
                is_zero_based=True,
            )
        )

        return pipeline

    except Exception as e:
        logger.error(f"Error in creating ProteinGym cleaning pipeline: {str(e)}")
        raise RuntimeError(f"Error in creating ProteinGym cleaning pipeline: {str(e)}")


def clean_proteingym_dms_substitutions_dataset(pipeline: Pipeline) -> Tuple[Pipeline, MutationDataset]:
    """Clean ProteinGym dataset using configurable pipeline

    Parameters
    ----------
    pipeline : Pipeline
        ProteinGym dataset cleaning pipeline

    Returns
    -------
    Tuple[Pipeline, MutationDataset]
        - Pipeline: The cleaned pipeline
        - MutationDataset: The cleaned ProteinGym dataset
    """
    try:
        # Run pipeline
        pipeline.execute()

        # Extract results
        proteingym_dms_substitutions_dataset_df, proteingym_dms_substitutions_ref_seq = pipeline.data
        proteingym_dms_substitutions_dataset = MutationDataset.from_dataframe(
            proteingym_dms_substitutions_dataset_df, proteingym_dms_substitutions_ref_seq
        )

        logger.info(
            f"Successfully cleaned ProteinGym dataset: {len(proteingym_dms_substitutions_dataset_df)} mutations from {len(proteingym_dms_substitutions_ref_seq)} proteins"
        )

        return pipeline, proteingym_dms_substitutions_dataset

    except Exception as e:
        logger.error(f"Error in running ProteinGym cleaning pipeline: {str(e)}")
        raise RuntimeError(f"Error in running ProteinGym cleaning pipeline: {str(e)}")
