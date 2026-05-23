from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from .base_config import BaseCleanerConfig
from .basic_cleaners import (
    apply_mutations_to_sequences,
    convert_to_mutation_dataset_format,
)
from .Chitosanase_custom_cleaners import parse_chitosanase_raw_file
from ..core.dataset import MutationDataset
from ..core.pipeline import Pipeline, create_pipeline

if TYPE_CHECKING:
    from typing import Any, Callable, Dict, List, Literal, Optional, Tuple, Union

__all__ = [
    "ChitosanaseCleanerConfig",
    "create_chitosanase_cleaner",
    "clean_chitosanase_dataset",
]


def __dir__() -> List[str]:
    return __all__


# Create module logger
logger = logging.getLogger(__name__)


@dataclass
class ChitosanaseCleanerConfig(BaseCleanerConfig):
    """
    Configuration class for Chitosanase dataset cleaner.

    This configuration holds dataset-specific settings used by the
    Chitosanase cleaning pipeline. The pipeline expects raw Chitosanase
    files where a small CSV is followed by a wild-type sequence
    separated by the token defined in ``wt_separator``.

    Attributes
    ----------
    infer_mut_workers : int
        Number of workers to use when inferring/applying mutations (default 16).
    pipeline_name : str
        Human-readable name for the pipeline used in logging and artifacts.
    wt_separator : str
        Substring token separating the CSV and wild-type sequence
        inside each Chitosanase raw file.
    """

    infer_mut_workers: int = 16
    pipeline_name: str = "Chitosanase"
    wt_separator: str = '">wt'

    def validate(self) -> None:
        super().validate()


def create_chitosanase_cleaner(
    dataset_or_path: Optional[Union[str, Path]] = None,
    config: Optional[Union[ChitosanaseCleanerConfig, Dict[str, Any], str, Path]] = None,
) -> Pipeline:
    """
    Create Chitosanase dataset cleaning pipeline

    Parameters
    ----------
    dataset_or_path : Optional[Union[pd.DataFrame, str, Path]]
        Path to a raw Chitosanase input file (or a DataFrame in other callers).
        The raw file must contain a CSV followed by the WT sequence.
    config : Optional[Union[ChitosanaseCleanerConfig, Dict[str, Any], str, Path]]
        Pipeline configuration. May be an instance of
        `ChitosanaseCleanerConfig`, a dict to merge with defaults, or a
        path to a JSON config file.

    Returns
    -------
    Pipeline
        Configured pipeline instance which can be executed with
        :meth:`Pipeline.execute()` to perform cleaning.

    Raises
    ------
    TypeError
        If `config` has an unsupported type.
    RuntimeError
        If pipeline creation fails for any reason.
    """
    # Handle config
    if config is None:
        final_config = ChitosanaseCleanerConfig()
    elif isinstance(config, ChitosanaseCleanerConfig):
        final_config = config
    elif isinstance(config, dict):
        final_config = ChitosanaseCleanerConfig().merge(config)
    elif isinstance(config, (str, Path)):
        final_config = ChitosanaseCleanerConfig.from_json(config)
    else:
        raise TypeError(f"config has invalid type: {type(config)}")

    if dataset_or_path is None:
        raise TypeError("dataset_or_path must be a Chitosanase file path")

    try:
        pipeline = create_pipeline(dataset_or_path, final_config.pipeline_name)

        # Add cleaning steps
        pipeline = (
            pipeline.delayed_then(parse_chitosanase_raw_file, wt_separator=final_config.wt_separator)
            .delayed_then(
                apply_mutations_to_sequences,
                sequence_column="sequence",
                mutation_column="mut_info",
                sequence_type="protein",
                is_zero_based=False,
                num_workers=final_config.infer_mut_workers,
            )
            .delayed_then(
                convert_to_mutation_dataset_format,
                name_column="name",
                mutation_column="mut_info",
                sequence_column="sequence",
                mutated_sequence_column="mut_seq",
                label_column="dTm",
                is_zero_based=False,
            )
        )
        return pipeline
    except Exception as e:
        logger.error(f"Error in creating Chitosanase cleaning pipeline: {str(e)}")
        raise RuntimeError(f"Error in creating Chitosanase cleaning pipeline: {str(e)}")


def clean_chitosanase_dataset(
    pipeline: Pipeline,
) -> Tuple[Pipeline, MutationDataset]:
    """
    Execute the cleaning pipeline and return the formatted dataset.

    This helper runs the provided :class:`Pipeline`, converts the resulting
    formatted DataFrame into a :class:`MutationDataset` and returns both the
    executed pipeline and the dataset object. Artifacts and diagnostic logs
    are produced by the pipeline and can be saved with
    ``pipeline.save_artifacts(path)``.

    Parameters
    ----------
    pipeline : Pipeline
        The configured cleaning pipeline to execute.

    Returns
    -------
    Tuple[Pipeline, MutationDataset]
        The executed pipeline and the resulting :class:`MutationDataset`.

    Raises
    ------
    RuntimeError
        When pipeline execution fails.
    """
    try:
        pipeline.execute()

        formatted_df, ref_dict = pipeline.data
        chitosanase_dataset = MutationDataset.from_dataframe(formatted_df, reference_sequences=ref_dict)

        logger.info(f"Successfully cleaned Chitosanase dataset: " f"{len(formatted_df)} mutations from {len(ref_dict)} proteins")
        return pipeline, chitosanase_dataset
    except Exception as e:
        logger.error(f"Error in running Chitosanase cleaning pipeline: {str(e)}")
        raise RuntimeError(f"Error in running Chitosanase cleaning pipeline: {str(e)}")
