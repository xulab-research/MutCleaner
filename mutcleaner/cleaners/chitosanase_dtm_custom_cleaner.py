# cleaners/chitosanase_dtm_custom_cleaner.py
from __future__ import annotations

import io
from pathlib import Path
import pandas as pd
from typing import TYPE_CHECKING

from ..core.pipeline import pipeline_step

if TYPE_CHECKING:
    from typing import List

__all__ = ["parse_chitosanase_raw_file"]


def __dir__() -> List[str]:
    return __all__


@pipeline_step
def parse_chitosanase_raw_file(file_path: str | Path, wt_separator: str = '">wt') -> pd.DataFrame:
    """Parse a raw Chitosanase input file and return the raw DataFrame.

    The raw file contains a CSV block followed by a wild-type sequence.
    This step only reads the file, splits the two blocks, parses the CSV into
    a DataFrame, and attaches the WT sequence as a constant column so that the
    downstream pipeline steps can do the actual cleaning and formatting.

    Parameters
    ----------
    file_path : str | pathlib.Path
        Path to the raw Chitosanase input file. The file should contain a
        CSV section then the WT sequence separated by ``wt_separator``.
    wt_separator : str, optional
        Substring that separates CSV and WT sequence blocks (default '\">wt').

    Returns
    -------
    pd.DataFrame
        Raw DataFrame parsed from the CSV block. It keeps the original input
        columns and adds ``wt_seq`` so downstream pipeline steps can create
        the sequence-related columns.

    Raises
    ------
    ValueError
        If the expected WT separator is not found.

    Examples
    --------
    >>> from pathlib import Path
    >>> df = parse_chitosanase_dtm_raw_file(Path("/path/to/Chitosanase_dTm_Dataset.csv"))
    >>> sorted(df.columns)
    ['Tm', 'aa_mut', 'wt_seq']
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
    df["wt_seq"] = wt_seq
    return df
