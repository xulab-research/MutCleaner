from __future__ import annotations

import io
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from .base_config import BaseCleanerConfig
from .basic_cleaners import (
    apply_mutations_to_sequences,
    convert_to_mutation_dataset_format,
)
from ..core.dataset import MutationDataset
from ..core.pipeline import Pipeline, create_pipeline

if TYPE_CHECKING:
    from typing import Any, Dict, Optional, Tuple, Union

__all__ = [
    "ChitosanaseCleanerConfig",
    "create_chitosanase_cleaner",
    "clean_chitosanase_dataset",
]

# Create module logger
logger = logging.getLogger(__name__)


@dataclass
class ChitosanaseCleanerConfig(BaseCleanerConfig):
    """
    Configuration class for Chitosanase dataset cleaner.
    """

    infer_mut_workers: int = 16
    pipeline_name: str = "Chitosanase"
    wt_separator: str = '">wt'

    def validate(self) -> None:
        super().validate()


def parse_chitosanase_raw_file(file_path: Union[str, Path], wt_separator: str = '">wt') -> pd.DataFrame:
    """
    Extract WT sequence and generate intermediate DataFrame.
    """
    with open(file_path, "r") as f:
        content = f.read()

    if wt_separator in content:
        parts = content.split(wt_separator)
        csv_text = parts[0].strip()
        wt_seq = parts[1].replace('"', "").replace(",", "").strip()
        wt_seq = "".join(wt_seq.split())
    else:
        raise ValueError(f"Cannot find WT sequence separator '{wt_separator}' in the expected format.")

    df = pd.read_csv(io.StringIO(csv_text))
    df["aa_mut"] = df["aa_mut"].astype(str).str.replace('"', "").str.strip()
    df = df.dropna(subset=["Tm"])

    wt_mask = df["aa_mut"] == "WT"
    if wt_mask.any():
        wt_tm = float(df[wt_mask]["Tm"].iloc[0])
        df["dTm"] = df["Tm"].astype(float) - wt_tm
        df = df[~wt_mask].copy()
    else:
        df["dTm"] = df["Tm"].astype(float)

    df["name"] = "Chitosanase"
    df["mut_info"] = df["aa_mut"]
    df["wt_seq"] = wt_seq
    df["sequence"] = wt_seq

    return df


def create_chitosanase_cleaner(
    dataset_or_path: Optional[Union[pd.DataFrame, str, Path]] = None,
    config: Optional[Union[ChitosanaseCleanerConfig, Dict[str, Any], str, Path]] = None,
) -> Pipeline:
    """Create Chitosanase dataset cleaning pipeline"""
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

    # Parse and read data
    if isinstance(dataset_or_path, (str, Path)):
        df_clean = parse_chitosanase_raw_file(dataset_or_path, wt_separator=final_config.wt_separator)
    elif isinstance(dataset_or_path, pd.DataFrame):
        df_clean = dataset_or_path
    else:
        raise TypeError("dataset_or_path must be pd.DataFrame or str/Path")

    try:
        pipeline = create_pipeline(df_clean, final_config.pipeline_name)

        # Add cleaning steps
        pipeline = pipeline.delayed_then(
            apply_mutations_to_sequences,
            sequence_column="sequence",
            mutation_column="mut_info",
            sequence_type="protein",
            is_zero_based=False,
            num_workers=final_config.infer_mut_workers,
        ).delayed_then(
            convert_to_mutation_dataset_format,
            name_column="name",
            mutation_column="mut_info",
            sequence_column="sequence",
            mutated_sequence_column="mut_seq",
            label_column="dTm",
            is_zero_based=False,
        )
        return pipeline
    except Exception as e:
        logger.error(f"Error in creating Chitosanase cleaning pipeline: {str(e)}")
        raise RuntimeError(f"Error in creating Chitosanase cleaning pipeline: {str(e)}")


def clean_chitosanase_dataset(
    pipeline: Pipeline,
) -> Tuple[Pipeline, MutationDataset]:
    """Clean Chitosanase dataset using configurable pipeline"""
    try:
        pipeline.execute()

        formatted_df, ref_dict = pipeline.data
        chitosanase_dataset = MutationDataset.from_dataframe(formatted_df, reference_sequences=ref_dict)

        logger.info(f"Successfully cleaned Chitosanase dataset: " f"{len(formatted_df)} mutations from {len(ref_dict)} proteins")
        return pipeline, chitosanase_dataset
    except Exception as e:
        logger.error(f"Error in running Chitosanase cleaning pipeline: {str(e)}")
        raise RuntimeError(f"Error in running Chitosanase cleaning pipeline: {str(e)}")
