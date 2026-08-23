# mutcleaner/cleaners/trpb_custom_cleaner.py

from __future__ import annotations

import pandas as pd

from typing import TYPE_CHECKING

from ..core.pipeline import pipeline_step

if TYPE_CHECKING:
    from typing import List

__all__ = ["standardize_trpb_mutation"]


def __dir__() -> List[str]:
    return __all__


@pipeline_step
def standardize_trpb_mutation_dataset(df: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize TrpB mutation information into unified mutation format.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw TrpB mutation dataset containing amino acid states in
        columns AA1-AA4.

    Returns
    -------
    pandas.DataFrame
        DataFrame with standardized mutation_name column and removed
        stop/WT mutations.
    """
    df = df.copy()

    positions = {
        "AA1": ("V", 182),
        "AA2": ("F", 183),
        "AA3": ("V", 226),
        "AA4": ("S", 227),
    }

    def build_mutation(row):
        mutations = []

        for col, (wt, pos) in positions.items():
            mut = row[col]

            if mut != wt:
                mutations.append(f"{wt}{pos}{mut}")

        return ",".join(mutations) if mutations else "WT"

    df["mutation_name"] = df.apply(build_mutation, axis=1)

    df = df.dropna(subset=["mutation_name"])
    return df
