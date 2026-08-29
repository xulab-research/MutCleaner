from __future__ import annotations

import pandas as pd
from tqdm import tqdm
from pathlib import Path
from typing import TYPE_CHECKING

from ..utils.sequence_io import load_sequences
from ..core.pipeline import pipeline_step

if TYPE_CHECKING:
    from typing import Union


@pipeline_step
def read_codon_dms_substitutions_dataset(
    data_path: Union[str, Path],
) -> pd.DataFrame:
    """
    Read and combine all assays in the Codon DMS Substitutions Dataset.

    The input can be either a directory or a ZIP archive containing MaveDB
    assay subdirectories. Each assay directory must contain ``data.csv`` and
    ``wt.fasta``. The wild-type sequence and assay directory name are added
    to each dataset before all assays are concatenated.

    Parameters
    ----------
    data_path : Union[str, Path]
        Path to the dataset directory or ZIP archive.

    Returns
    -------
    pd.DataFrame
        Combined dataset containing all assays with ``name`` and
        ``wt_sequence`` columns.

    Raises
    ------
    FileNotFoundError
        If ``data_path`` does not exist.
    ValueError
        If no valid assay directories are found or an assay FASTA file does
        not contain exactly one sequence.
    """
    import tempfile
    import zipfile

    data_path = Path(data_path)

    if not data_path.exists():
        raise FileNotFoundError(f"Data path does not exist: {data_path}")

    temp_dir = None

    try:
        if data_path.suffix.lower() == ".zip":
            tqdm.write(f"Extracting Codon DMS Substitutions Dataset: {data_path}")
            temp_dir = tempfile.TemporaryDirectory(prefix="codon_dms_")
            working_dir = Path(temp_dir.name)

            with zipfile.ZipFile(data_path, "r") as zip_ref:
                zip_ref.extractall(working_dir)
        elif data_path.is_dir():
            working_dir = data_path
        else:
            raise ValueError(f"Data path must be a directory or ZIP file: {data_path}")

        assay_dirs = sorted(
            path
            for path in working_dir.rglob("*")
            if path.is_dir()
            and (path / "data.csv").exists()
            and (path / "wt.fasta").exists()
        )

        if not assay_dirs:
            raise ValueError(f"No assay directories found in {data_path}")

        tqdm.write(f"Found {len(assay_dirs)} MaveDB assays to process")

        datasets = []

        for assay_dir in tqdm(assay_dirs, desc="Reading Codon DMS assays"):
            wt_sequences = load_sequences(assay_dir / "wt.fasta")

            if len(wt_sequences) != 1:
                raise ValueError(
                    f"Expected one wild-type sequence in {assay_dir / 'wt.fasta'}, "
                    f"found {len(wt_sequences)}"
                )

            df = pd.read_csv(assay_dir / "data.csv")
            df["wt_sequence"] = next(iter(wt_sequences.values()))
            df["name"] = assay_dir.name
            datasets.append(df)

        dataset = pd.concat(datasets, ignore_index=True)

        tqdm.write(
            f"Loaded Codon DMS Substitutions Dataset: "
            f"{len(dataset)} mutation records from {len(assay_dirs)} assays"
        )

        return dataset

    finally:
        if temp_dir is not None:
            temp_dir.cleanup()
