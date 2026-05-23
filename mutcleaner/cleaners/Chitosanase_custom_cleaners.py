from __future__ import annotations

import io
from pathlib import Path
import pandas as pd
from typing import TYPE_CHECKING

from ..core.pipeline import pipeline_step

if TYPE_CHECKING:
    from typing import Any, Dict, List, Sequence, Tuple, Union

__all__ = ["parse_chitosanase_raw_file"]


def __dir__() -> List[str]:
    return __all__


@pipeline_step
def parse_chitosanase_raw_file(file_path: str | Path, wt_separator: str = '">wt') -> pd.DataFrame:
    """
    Parse the raw Chitosanase text block and generate an intermediate DataFrame.

    The raw input contains a CSV followed by a wild-type sequence,
    separated by a custom token. This helper reads the file content, extracts
    the CSV and WT sequence, then returns the intermediate dataframe.

    Parameters
    ----------
    file_path : str or pathlib.Path
        Path to the Chitosanase input file.
    wt_separator : str, default '">wt'
        Substring separating the CSV block from the WT sequence block.

    Returns
    -------
    pd.DataFrame
        A DataFrame with standard intermediate columns ready for the pipeline.
    """
    with open(file_path, "r") as f:
        raw_text = f.read()

    if wt_separator in raw_text:
        parts = raw_text.split(wt_separator)
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
