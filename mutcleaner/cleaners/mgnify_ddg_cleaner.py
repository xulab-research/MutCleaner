# mutcleaner/cleaners/mgnify_ddg_cleaner.py
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from .base_config import BaseCleanerConfig
from .basic_cleaners import (
    aggregate_labels_by_name,
    convert_data_types,
    convert_to_mutation_dataset_format,
    extract_and_rename_columns,
    infer_mutations_from_sequences,
    read_dataset,
)
from ..core.dataset import MutationDataset
from ..core.pipeline import Pipeline, create_pipeline

if TYPE_CHECKING:
    from typing import Any, Dict, List, Literal, Optional, Tuple, Union

__all__ = [
    "MGnifyddGCleanerConfig",
    "create_mgnify_ddg_cleaner",
    "clean_mgnify_ddg_dataset",
]


def __dir__() -> List[str]:
    return __all__


logger = logging.getLogger(__name__)


@dataclass
class MGnifyddGCleanerConfig(BaseCleanerConfig):
    """
    Configuration class for MGnify protein stability dataset cleaner.
    Inherits from BaseCleanerConfig and adds MGnify-specific configuration options.

    This configuration is specifically optimized for the MGnify cluster-subset
    schema which pre-calculates ddG values from sequence pairs without extra
    environmental columns (such as pH or temperature).

    Attributes
    ----------
    column_mapping : Dict[str, str]
        Mapping from source columns to MutCleaner standard columns.
        Default targets are: 'WT_name' -> 'name', 'wt_seq' -> 'wt_seq',
        'mut_seq' -> 'mut_seq', 'ddG' -> 'label'.
    type_conversions : Dict[str, str]
        Data type conversion specifications. Default forces 'label' to 'float32'.
    infer_mut_workers : int
        Number of parallel workers used for mutation inference. Set to -1 to use all CPUs.
    aggregation_strategy : Literal["mean", "first", "nearest"]
        Strategy to aggregate labels sharing identical mutation keys. Defaults to "mean".
    nearest_by : List[Tuple[str, float]]
        Keep mutation by distance metrics. Unused for pure computational subsets.
    label_columns : List[str]
        List of score/label columns to process.
    primary_label_column : str
        Primary score column for the dataset package.
    """

    column_mapping: Dict[str, str] = field(
        default_factory=lambda: {
            "WT_name": "name",
            "wt_seq": "wt_seq",
            "mut_seq": "mut_seq",
            "ddG": "label",
        }
    )

    type_conversions: Dict[str, str] = field(default_factory=lambda: {"label": "float32"})

    infer_mut_workers: int = 16

    aggregation_strategy: Literal["mean", "first", "nearest"] = "mean"
    nearest_by: List[Tuple[str, float]] = field(default_factory=list)

    label_columns: List[str] = field(default_factory=lambda: ["label"])
    primary_label_column: str = "label"

    pipeline_name: str = "MGnify-ddG"

    def __post_init__(self):
        self.nearest_by = [(str(col), float(target)) for col, target in self.nearest_by]

    def validate(self) -> None:
        """
        Validate MGnify-specific configuration parameters.

        Raises
        ------
        ValueError
            If label columns are misconfigured or required mappings are missing.
        """
        super().validate()

        if not self.label_columns:
            raise ValueError("label_columns cannot be empty")

        if self.primary_label_column not in self.label_columns:
            raise ValueError(f"primary_label_column '{self.primary_label_column}' " f"must be contained in label_columns {self.label_columns}")

        required_mappings = {"name", "wt_seq", "mut_seq", "label"}
        missing = required_mappings - set(self.column_mapping.values())
        if missing:
            raise ValueError(f"Missing required target column mappings inside column_mapping: {missing}")


def create_mgnify_ddg_cleaner(
    dataset_or_path: Optional[Union[pd.DataFrame, str, Path]] = None,
    config: Optional[Union[MGnifyddGCleanerConfig, Dict[str, Any], str, Path]] = None,
) -> Pipeline:
    """
    Create MGnify protein stability dataset cleaning pipeline.

    This function pieces together modular basic cleaners to standardize the MGnify
    subsets into regularized training matrices for deep learning architectures.

    Parameters
    ----------
    dataset_or_path : Optional[Union[pd.DataFrame, str, Path]], default=None
        Raw dataset DataFrame or file path to the raw MGnify_ddG_Dataset.csv file.
    config : Optional[Union[MGnifyddGCleanerConfig, Dict[str, Any], str, Path]]
        Configuration instance, parameter dictionary, or path to a JSON config.

    Returns
    -------
    Pipeline
        The constructed cleaning pipeline instance populated with delayed steps.
    """
    if config is None:
        final_config = MGnifyddGCleanerConfig()
    elif isinstance(config, MGnifyddGCleanerConfig):
        final_config = config
    elif isinstance(config, dict):
        default_config = MGnifyddGCleanerConfig()
        final_config = default_config.merge(config)
    elif isinstance(config, (str, Path)):
        final_config = MGnifyddGCleanerConfig.from_json(config)
    else:
        raise TypeError(f"config must be MGnifyddGCleanerConfig, dict, str, Path or None, got {type(config)}")

    logger.info(f"MGnify_ddG dataset will be cleaned with pipeline: {final_config.pipeline_name}")
    logger.debug(f"Configuration summary:\n{final_config.get_summary()}")

    try:
        pipeline = create_pipeline(dataset_or_path, final_config.pipeline_name)
        if isinstance(dataset_or_path, (str, Path)):
            pipeline.then(read_dataset)
        pipeline = (
            pipeline.delayed_then(
                extract_and_rename_columns,
                column_mapping=final_config.column_mapping,
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
        logger.error(f"Failed to initialize MGnify-ddG cleaning pipeline: {str(e)}")
        raise RuntimeError(f"Failed to initialize MGnify-ddG cleaning pipeline: {str(e)}")


def clean_mgnify_ddg_dataset(
    pipeline: Pipeline,
) -> Tuple[Pipeline, MutationDataset]:
    """
    Execute the MGnify_ddG dataset cleaning pipeline and package into a MutationDataset.

    Parameters
    ----------
    pipeline : Pipeline
        The pre-constructed MGnify_ddG dataset cleaning pipeline.

    Returns
    -------
    Tuple[Pipeline, MutationDataset]
        - pipeline: The executed pipeline instance.
        - mgnify_ddg_dataset: Standardized MutationDataset object containing data.csv,
          wt.fasta, and metadata.json.
    """
    try:
        pipeline.execute()
        formatted_df, ref_dict = pipeline.data

        mgnify_dataset = MutationDataset.from_dataframe(formatted_df, reference_sequences=ref_dict)

        logger.info(f"Successfully executed MGnify-ddG pipeline: {len(formatted_df)} mutations across {len(ref_dict)} unique proteins.")

        return pipeline, mgnify_dataset

    except Exception as e:
        logger.error(f"Error encountered during MGnify-ddG pipeline execution: {str(e)}")
        raise RuntimeError(f"Error encountered during MGnify-ddG pipeline execution: {str(e)}")
