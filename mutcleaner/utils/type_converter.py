"""
Type conversion utilities for pandas DataFrames.

This module provides flexible and efficient data type conversion capabilities
with support for pandas, numpy, and Python built-in types.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from tqdm import tqdm
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Callable, Dict, List, Tuple, Type, Union

__all__ = [
    "convert_data_types",
    "convert_data_types_batch",
    "convert_to_boolean",
    "get_conversion_function",
    "normalize_type_conversions",
]


def __dir__() -> List[str]:
    return __all__


def normalize_type_conversions(type_conversions: Dict[str, Union[str, Type, np.dtype]], optimize_memory: bool) -> Dict[str, Tuple[str, Callable]]:
    """
    Normalize type conversion mapping to standardized format.

    Returns
    -------
    Dict[str, Tuple[str, Callable]]
        Mapping of {column: (type_name, conversion_function)}
    """
    normalized = {}

    for col, target_type in type_conversions.items():
        type_name, conversion_func = get_conversion_function(target_type, optimize_memory)
        normalized[col] = (type_name, conversion_func)

    return normalized


def get_conversion_function(target_type: Union[str, Type, np.dtype], optimize_memory: bool) -> Tuple[str, Callable]:
    """
    Get appropriate conversion function for target type.

    Parameters
    ----------
    target_type : Union[str, Type, np.dtype]
        Target data type
    optimize_memory : bool
        Whether to use memory-optimized types

    Returns
    -------
    Tuple[str, Callable]
        Type name and conversion function
    """
    # Handle numpy types
    if isinstance(target_type, type) and issubclass(target_type, np.number):
        return str(target_type), lambda series, errors: pd.to_numeric(series, errors=errors).astype(target_type)

    if isinstance(target_type, np.dtype):
        return str(target_type), lambda series, errors: pd.to_numeric(series, errors=errors).astype(target_type)

    # Handle string type names
    if isinstance(target_type, str):
        target_type_lower = target_type.lower()

        # Numeric types
        if target_type_lower in ["float", "float64"]:
            dtype = np.float32 if optimize_memory else np.float64
            type_name = "float32" if optimize_memory else "float64"
            return type_name, lambda series, errors: pd.to_numeric(series, errors=errors).astype(dtype)

        elif target_type_lower in ["int", "int64"]:
            return "Int64", lambda series, errors: pd.to_numeric(series, errors=errors, downcast="integer" if optimize_memory else None).astype("Int64")

        elif target_type_lower == "int32":
            return "Int32", lambda series, errors: pd.to_numeric(series, errors=errors).astype("Int32")

        # Pandas extension types
        elif target_type in ["Int64", "Int32", "Int16", "Int8"]:
            return target_type, lambda series, errors: pd.to_numeric(series, errors=errors).astype(target_type)

        elif target_type in ["Float64", "Float32"]:
            return target_type, lambda series, errors: pd.to_numeric(series, errors=errors).astype(target_type)

        elif target_type_lower in ["str", "string"]:
            return "string", lambda series, _: series.astype("string")

        elif target_type_lower in ["bool", "boolean"]:
            return "boolean", lambda series, errors: convert_to_boolean(series, errors)

        elif target_type_lower == "category":
            return "category", lambda series, _: series.astype("category")

        elif target_type_lower == "datetime":
            return "datetime64[ns]", lambda series, errors: pd.to_datetime(series, errors=errors)

        # Direct pandas astype
        else:
            return target_type, lambda series, _: series.astype(target_type)

    # Handle Python built-in types
    elif target_type == float:
        dtype = np.float32 if optimize_memory else np.float64
        type_name = "float32" if optimize_memory else "float64"
        return type_name, lambda series, errors: pd.to_numeric(series, errors=errors).astype(dtype)

    elif target_type == int:
        return "Int64", lambda series, errors: pd.to_numeric(series, errors=errors).astype("Int64")

    elif target_type == str:
        return "string", lambda series, _: series.astype("string")

    elif target_type == bool:
        return "boolean", lambda series, errors: convert_to_boolean(series, errors)

    else:
        return str(target_type), lambda series, _: series.astype(target_type)


def convert_to_boolean(series: pd.Series, errors: str) -> pd.Series:
    """
    Intelligent boolean conversion handling various boolean representations.

    Parameters
    ----------
    series : pd.Series
        Input series to convert
    errors : str
        Error handling strategy

    Returns
    -------
    pd.Series
        Boolean series
    """
    if errors == "coerce":
        bool_map = {
            "true": True,
            "false": False,
            "yes": True,
            "no": False,
            "1": True,
            "0": False,
            "y": True,
            "n": False,
            "t": True,
            "f": False,
        }

        str_series = series.astype(str).str.lower().str.strip()
        result = str_series.map(bool_map)

        # Handle numeric values
        numeric_mask = pd.to_numeric(series, errors="coerce").notna()
        result.loc[numeric_mask] = pd.to_numeric(series.loc[numeric_mask], errors="coerce").astype(bool)

        return result.astype("boolean")
    else:
        return series.astype("boolean")


def convert_data_types(
    dataset: pd.DataFrame,
    type_conversions: Dict[str, Union[str, Type, np.dtype]],
    handle_errors: str = "coerce",
    optimize_memory: bool = True,
) -> pd.DataFrame:
    """
    Convert data types for specified columns with enhanced type support.

    Parameters
    ----------
    dataset : pd.DataFrame
        Input dataset with columns to be converted
    type_conversions : Dict[str, Union[str, Type, np.dtype]]
        Type conversion mapping in format {column_name: target_type}
        Supported formats:
        - String types: 'float', 'int', 'str', 'category', 'bool', 'datetime'
        - Numpy types: np.float32, np.float64, np.int32, np.int64, etc.
        - Pandas types: 'Int64', 'Float64', 'string', 'boolean'
        - Python types: float, int, str, bool
    handle_errors : str, default='coerce'
        Error handling strategy: 'raise', 'coerce', or 'ignore'
    optimize_memory : bool, default=True
        Whether to automatically optimize memory usage

    Returns
    -------
    pd.DataFrame
        Dataset with converted data types
    """
    normalized_conversions = normalize_type_conversions(type_conversions, optimize_memory)

    missing_columns = set(normalized_conversions.keys()) - set(dataset.columns)
    if missing_columns:
        tqdm.write(f"Warning: Columns not found: {missing_columns}")
        normalized_conversions = {k: v for k, v in normalized_conversions.items() if k not in missing_columns}

    if not normalized_conversions:
        tqdm.write("No valid columns to convert")
        return dataset.copy()

    result = dataset.copy()
    conversion_results = {}

    for col, (type_name, conversion_func) in normalized_conversions.items():
        try:
            converted_series = conversion_func(result[col], handle_errors)
            conversion_results[col] = converted_series
            tqdm.write(f"Converted '{col}' to {type_name}")
        except Exception as e:
            if handle_errors == "raise":
                raise ValueError(f"Failed to convert column '{col}' to {type_name}: {e}")
            else:
                tqdm.write(f"Warning: Failed to convert column '{col}' to {type_name}: {e}")

    for col, converted_series in conversion_results.items():
        result[col] = converted_series

    return result


def convert_data_types_batch(
    dataset: pd.DataFrame,
    type_conversions: Dict[str, Union[str, Type, np.dtype]],
    handle_errors: str = "coerce",
    optimize_memory: bool = True,
    chunk_size: int = 10000,
) -> pd.DataFrame:
    """
    Batch conversion version for large datasets.

    Parameters
    ----------
    dataset : pd.DataFrame
        Input dataset
    type_conversions : Dict[str, Union[str, Type, np.dtype]]
        Type conversion mapping
    handle_errors : str, default='coerce'
        Error handling strategy
    optimize_memory : bool, default=True
        Whether to optimize memory usage
    chunk_size : int, default=10000
        Chunk size for processing large datasets

    Returns
    -------
    pd.DataFrame
        Dataset with converted data types
    """
    if len(dataset) <= chunk_size:
        return convert_data_types(dataset, type_conversions, handle_errors, optimize_memory)

    tqdm.write(f"Converting data types in chunks of {chunk_size}...")

    chunks = []
    for i in range(0, len(dataset), chunk_size):
        chunk = dataset.iloc[i : i + chunk_size]
        converted_chunk = convert_data_types(chunk, type_conversions, handle_errors, optimize_memory)
        chunks.append(converted_chunk)

    result = pd.concat(chunks, ignore_index=True)
    return result
