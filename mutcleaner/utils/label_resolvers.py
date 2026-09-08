# label_resolvers.py

r"""
Label resolvers for aggregating per-group target columns.

This module provides small, composable **resolvers**—callables that, given a
pandas ``DataFrame`` *group* and a list of *label columns*, return a
``Series`` of resolved label values. Resolvers are intended to be used inside
grouped operations (e.g., ``DataFrameGroupBy.apply``) where multiple rows per
entity must be collapsed to a single, consistent set of labels.

The public entry point :func:`make_resolver` constructs a resolver from a
strategy name (e.g., ``"mean"``, ``"first"``, ``"nearest"``) or accepts a
user-supplied callable. The :func:`nearest_resolver_factory` enables
lexicographic, weighted "nearest row" selection across multiple numeric
criteria columns.

Notes
-----
- **Resolver signature**: ``Resolver = Callable[[pd.DataFrame, list[str]], pd.Series]``.
  The returned ``Series`` must align with the provided ``label_cols`` order.
- **Missing values**: For the ``"nearest"`` strategy, missing (NaN) values in
  criterion columns are treated as ``+inf`` distance so such rows never win.
- **Column order**: When using mappings (``dict``) for criteria, Python 3.7+
  preserves insertion order. That order determines lexicographic priority.
- **Type expectations**:
  - ``"mean"`` requires numeric label columns; non-numeric values raise
    ``ValueError``.
  - ``"nearest"`` requires numeric criterion columns; non-numeric values raise
    ``ValueError``.
- **Error handling**: Missing required columns raise ``KeyError``; unknown
  strategy names raise ``ValueError``.

See Also
--------
pandas.DataFrame.groupby : Grouping rows for split-apply-combine workflows.
numpy.lexsort : Related concept for lexicographic ordering.

"""

from __future__ import annotations

from typing import Callable, List, Mapping, cast, TYPE_CHECKING
import pandas as pd
import numpy as np

if TYPE_CHECKING:
    from typing import Dict, Sequence, Tuple

__all__ = ["make_resolver"]


def __dir__() -> List[str]:
    return __all__


# Public type for resolvers
Resolver = Callable[[pd.DataFrame, List[str]], pd.Series]


def mean_resolver(group: pd.DataFrame, label_cols: List[str]) -> pd.Series:
    """
    Compute the numeric mean for each label column within a group.

    Parameters
    ----------
    group : pandas.DataFrame
        A single group of rows (as provided by ``DataFrameGroupBy.apply``).
    label_cols : list of str
        Column names whose values will be averaged within the group.

    Returns
    -------
    pandas.Series
        A Series containing the per-label mean values. The index matches
        ``label_cols`` order.

    Raises
    ------
    ValueError
        If any label column is non-numeric.

    Notes
    -----
    - NaN values are ignored by ``DataFrame.mean`` (`numeric_only=True`).
    - The result is reindexed to preserve the input order of ``label_cols``.
    """
    for c in label_cols:
        if not pd.api.types.is_numeric_dtype(group[c]):
            raise ValueError(f"label column '{c}' must be numeric for strategy='mean'")
    s = group[label_cols].mean(numeric_only=True)
    return s.reindex(label_cols)


def first_resolver(group: pd.DataFrame, label_cols: List[str]) -> pd.Series:
    """
    Take label values from the first row of the group.

    Parameters
    ----------
    group : pandas.DataFrame
        A single group of rows (as provided by ``DataFrameGroupBy.apply``).
    label_cols : list of str
        Column names to extract from the first row.

    Returns
    -------
    pandas.Series
        A Series containing the first-row values for the given labels.
        The index matches ``label_cols`` order.

    Notes
    -----
    - No type checks are performed; values are taken as-is from the first row.
    """
    return group.iloc[0][label_cols]


def nearest_resolver_factory(
    criteria: Mapping[str, float] | Sequence[Tuple[str, float]],
    weights: Mapping[str, float] | Sequence[Tuple[str, float]] | None = None,
) -> Resolver:
    """
    Create a resolver that selects the row minimizing a lexicographic,
    weighted distance across multiple numeric columns.

    The produced resolver compares rows by building a distance tuple:
    ``(w1 * |c1 - t1|, w2 * |c2 - t2|, ...)`` for the specified criteria
    and then picks the row with the smallest tuple (lexicographic order).
    Missing values (NaN) are treated as ``+inf`` so those rows never win.

    Parameters
    ----------
    criteria : mapping[str, float] or sequence of (str, float)
        Target values for the distance computation. Keys (or first elements
        of pairs) are column names, values are target scalars. If a mapping
        is provided, insertion order determines priority (Python 3.7+).
    weights : mapping[str, float] or sequence of (str, float), optional
        Optional per-column weights. If omitted, all weights default to 1.0.
        Must reference the same columns as ``criteria`` if provided.

    Returns
    -------
    Resolver
        A callable of signature ``(group: DataFrame, label_cols: list[str]) -> Series``
        that returns the label values from the nearest row in the group.

    Raises
    ------
    KeyError
        If a required criterion column is not present in the group.
    ValueError
        If a criterion column is non-numeric.

    Examples
    --------
    >>> factory = nearest_resolver_factory({"temperature": 25.0, "pH": 7.4})
    >>> # later inside a groupby-apply:
    >>> # factory(group_df, ["ddG", "dTm"])
    """
    # normalize criteria into ordered list
    if isinstance(criteria, Mapping):
        crit_pairs = list(criteria.items())
    else:
        crit_pairs = list(criteria)

    # normalize weights into dict
    weight_map: Dict[str, float] = {}
    if weights is not None:
        weight_map = dict(weights.items()) if isinstance(weights, Mapping) else dict(weights)

    def _resolver(group: pd.DataFrame, label_cols: List[str]) -> pd.Series:
        # validate columns & numeric types
        for col, _ in crit_pairs:
            if col not in group.columns:
                raise KeyError(f"nearest_by column '{col}' not found in group")
            if not pd.api.types.is_numeric_dtype(group[col]):
                raise ValueError(f"nearest_by column '{col}' must be numeric")

        # build a distance tuple for each row
        def dist_tuple(idx) -> Tuple[float, ...]:
            row = group.loc[idx]
            parts: List[float] = []
            for col, target in crit_pairs:
                w = float(weight_map.get(col, 1.0))
                val = row[col]
                parts.append(np.inf if pd.isna(val) else abs(float(val) - float(target)) * w)
            return tuple(parts)

        distances = group.index.to_series().apply(dist_tuple)
        best_idx = distances.idxmin()
        return cast(pd.Series, group.loc[best_idx].loc[list(label_cols)])

    return _resolver


# --- Registry (string -> resolver) ---
RESOLVERS: Dict[str, Resolver] = {
    "mean": mean_resolver,
    "first": first_resolver,
    # 'nearest' is created via factory, see make_resolver() below
}


def make_resolver(
    strategy: str | Resolver,
    *,
    nearest_by: Mapping[str, float] | Sequence[Tuple[str, float]] | None = None,
    nearest_weights: Mapping[str, float] | Sequence[Tuple[str, float]] | None = None,
) -> Resolver:
    """
    Create or return a label resolver from a strategy name or callable.

    Parameters
    ----------
    strategy : {"mean", "first", "nearest"} or Callable
        Strategy identifier or a custom resolver. If a callable is provided,
        it must have signature ``(group: DataFrame, label_cols: list[str]) -> Series``.
    nearest_by : mapping[str, float] or sequence of (str, float), optional
        Criteria for the ``"nearest"`` strategy. Required when
        ``strategy="nearest"``. See :func:`nearest_resolver_factory`.
    nearest_weights : mapping[str, float] or sequence of (str, float), optional
        Weights for the ``"nearest"`` strategy. See :func:`nearest_resolver_factory`.

    Returns
    -------
    Resolver
        A resolver callable that can be passed into higher-level aggregation code.

    Raises
    ------
    ValueError
        If ``strategy`` is unknown or required parameters are missing for
        the ``"nearest"`` strategy.

    Examples
    --------
    >>> res = make_resolver("mean")
    >>> # or nearest:
    >>> res = make_resolver("nearest", nearest_by={"temperature": 25.0})
    >>> # or custom:
    >>> def pick_row(g, labels):
    ...     return g.loc[g["score"].idxmax(), labels]
    >>> res = make_resolver(pick_row)
    """
    if callable(strategy):
        return strategy

    key = str(strategy).lower()
    if key == "nearest":
        if nearest_by is None:
            raise ValueError("strategy='nearest' requires `nearest_by`")
        return nearest_resolver_factory(nearest_by, nearest_weights)

    if key in RESOLVERS:
        return RESOLVERS[key]

    raise ValueError(f"Unknown strategy: {strategy!r}")
