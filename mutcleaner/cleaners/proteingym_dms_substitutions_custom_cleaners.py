# mutcleaner/cleaners/proteingym_dms_substitutions_custom_cleaners.py
from __future__ import annotations

import pandas as pd
from pathlib import Path
from typing import TYPE_CHECKING

from ..core.pipeline import multiout_step

if TYPE_CHECKING:
    from typing import List, Tuple, Union

__all__ = ["read_proteingym_dms_substitutions_data"]


def __dir__() -> List[str]:
    return __all__


# Protein Gym data reader function
@multiout_step(main="success", failed="failed")
def read_proteingym_dms_substitutions_data(
    data_path: Union[str, Path],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Read and combine multiple ProteinGym datasets from a directory or zip file.

    ProteinGym datasets are stored as individual CSV files, one per protein.
    This function combines them into a single DataFrame for unified processing.
    Each file contains columns: mutant, mutated_sequence, DMS_score, and various prediction methods.

    Parameters
    ----------
    data_path : Union[str, Path]
        Path to directory containing ProteinGym CSV files or path to zip file

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        (success_dataframe, failed_dataframe) - successfully processed data and failed file info

    Raises
    ------
    FileNotFoundError
        If data_path does not exist
    ValueError
        If no CSV files found or required columns missing

    Examples
    --------
    Process directory of ProteinGym CSV files:
    >>> success_df, failed_df = read_proteingym_dms_substitutions_data("DMS_ProteinGym_substitutions/")

    Process zip file:
    >>> success_df, failed_df = read_proteingym_dms_substitutions_data("DMS_ProteinGym_substitutions.zip")
    """
    import shutil
    import tempfile
    import zipfile

    from tqdm import tqdm

    data_path = Path(data_path)

    if not data_path.exists():
        raise FileNotFoundError(f"Data path does not exist: {data_path}")

    temp_dir = None

    # Handle zip file extraction
    if data_path.suffix.lower() == ".zip":
        tqdm.write(f"Extracting ProteinGym zip file: {data_path}")

        # Create temporary directory
        temp_dir = Path(tempfile.mkdtemp(prefix="proteingym_"))

        try:
            # Extract zip file
            with zipfile.ZipFile(data_path, "r") as zip_ref:
                zip_ref.extractall(temp_dir)

            # Find the actual data directory in extracted content
            # Look for directories containing CSV files
            data_dirs = []
            for item in temp_dir.rglob("*"):
                if item.is_dir() and any(item.glob("*.csv")):
                    data_dirs.append(item)

            if not data_dirs:
                raise ValueError("No directories with CSV files found in zip")

            # Use the directory with most CSV files (main dataset directory)
            working_dir = max(data_dirs, key=lambda d: len(list(d.glob("*.csv"))))
            tqdm.write(f"Using directory: {working_dir.name}")

        except Exception as e:
            if temp_dir and temp_dir.exists():
                shutil.rmtree(temp_dir)
            raise RuntimeError(f"Failed to extract zip file: {e}")

    else:
        # Direct directory processing
        working_dir = data_path

    if not working_dir.is_dir():
        raise ValueError(f"Working directory is not a directory: {working_dir}")

    # Find all CSV files
    csv_files = list(working_dir.glob("*.csv"))
    if not csv_files:
        error_msg = f"No CSV files found in {working_dir}"
        if temp_dir:
            shutil.rmtree(temp_dir)
        raise ValueError(error_msg)

    tqdm.write(f"Found {len(csv_files)} ProteinGym CSV files to process")

    combined_data = []
    failed_data = []

    try:
        for csv_file in tqdm(csv_files, desc="Processing ProteinGym files"):
            try:
                # Extract protein name from filename (without extension)
                protein_name = csv_file.stem

                # Read CSV file
                df = pd.read_csv(csv_file)

                # Check if required ProteinGym columns exist
                required_cols = ["mutant", "mutated_sequence", "DMS_score"]
                missing_cols = [col for col in required_cols if col not in df.columns]
                if missing_cols:
                    error_msg = f"Missing columns {missing_cols}"
                    tqdm.write(f"Warning: {csv_file.name} {error_msg}, skipping")
                    failed_data.append(
                        {
                            "filename": csv_file.name,
                            "protein_name": protein_name,
                            "error_type": "missing_columns",
                            "error_message": error_msg,
                            "missing_columns": str(missing_cols),
                        }
                    )
                    continue

                # Add protein name column
                df["name"] = protein_name

                # Reorder columns: put standard columns first
                standard_columns = ["name", "mutant", "mutated_sequence", "DMS_score"]
                other_columns = [col for col in df.columns if col not in standard_columns]
                final_columns = standard_columns + other_columns
                df = df[final_columns]

                combined_data.append(df)

            except Exception as e:
                error_msg = f"{type(e).__name__}: {str(e)}"
                tqdm.write(f"Error processing {csv_file.name}: {error_msg}")
                failed_data.append(
                    {
                        "filename": csv_file.name,
                        "protein_name": csv_file.stem,
                        "error_type": type(e).__name__,
                        "error_message": error_msg,
                        "missing_columns": None,
                    }
                )
                continue

    finally:
        # Cleanup temporary directory if created
        if temp_dir and temp_dir.exists():
            tqdm.write(f"Cleaning up temporary directory: {temp_dir}")
            shutil.rmtree(temp_dir)

    if not combined_data:
        raise ValueError("No data could be processed from any CSV files")

    # Combine all successful dataframes
    success_df = pd.concat(combined_data, ignore_index=True)

    # Create failed dataframe
    failed_df = pd.DataFrame(failed_data) if failed_data else pd.DataFrame()

    # Log processing results
    total_files = len(csv_files)
    success_files = len(combined_data)
    failed_files = len(failed_data)

    tqdm.write(f"Successfully processed {success_files}/{total_files} files")
    tqdm.write(f"Combined ProteinGym dataset shape: {success_df.shape}")

    if failed_files > 0:
        tqdm.write(f"Failed to process {failed_files} files")
        if not failed_df.empty:
            error_types = failed_df["error_type"].value_counts()
            tqdm.write(f"Error types: {dict(error_types)}")

    return success_df, failed_df
