# mutcleaner/utils/raw_data_downloader.py
from __future__ import annotations

import hashlib
import json
import requests
import sys
import time
from pathlib import Path
from tqdm import tqdm
from typing import TYPE_CHECKING
from urllib.parse import urlparse

from .data_source import DATASETS

if TYPE_CHECKING:
    from typing import Dict, List, Literal, Optional, Sequence, Union

__all__ = [
    "download",
    "download_protein_cdna_proteolysis_source_file",
    "download_proteingym_source_file",
    "download_human_domainome_source_file",
    "download_rbd_antibody_source_file",
    "download_rbd_ace2_source_file",
    "download_antitoxin_pard3_source_file",
    "download_trpb_source_file",
    "download_ctxm_source_file",
    "download_human_myoglobin_source_file",
    "download_archstabms1e10_source_file",
    "download_ddg_dtm_source_file",
    "download_chitosanase_dtm_source_file", 
    "download_mgnify_ddg_source_file",
]


def __dir__() -> List[str]:
    return __all__


HF_ENDPOINTS: Sequence[str] = (
    "https://huggingface.co/",
    "https://hf-mirror.com/",
)

CONFIG_DIR = Path(sys.prefix) / ".mutcleaner"
CONFIG_FILE = CONFIG_DIR / "config.json"
CONFIG_KEY = "hf_endpoint"


def _load_cached_endpoint() -> Optional[str]:
    """Load cached Hugging Face endpoint from configuration."""
    if CONFIG_FILE.exists():
        try:
            with CONFIG_FILE.open() as f:
                config = json.load(f)
                return config.get(CONFIG_KEY)
        except json.JSONDecodeError:
            return None
    return None


def _save_cached_endpoint(endpoint: str) -> None:
    """Save Hugging Face endpoint to configuration."""
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    config = {}
    if CONFIG_FILE.exists():
        try:
            with CONFIG_FILE.open() as f:
                config = json.load(f)
        except json.JSONDecodeError:
            config = {}
    config[CONFIG_KEY] = endpoint
    CONFIG_FILE.write_text(json.dumps(config, indent=2), encoding="utf-8")


def _reachable(url: str, *, timeout: int = 4) -> bool:
    try:
        response = requests.head(url, timeout=timeout)
        response.raise_for_status()
        return True
    except requests.exceptions.RequestException:
        return False


def resolve_hf_endpoint() -> str:
    """Resolve the Hugging Face endpoint to use based on cached configuration."""
    cached_endpoint = _load_cached_endpoint()
    if cached_endpoint is not None and _reachable(cached_endpoint):
        return cached_endpoint

    for endpoint in HF_ENDPOINTS:
        if _reachable(endpoint):
            _save_cached_endpoint(endpoint)
            return endpoint

    # return huggingface.co as default if no reachable endpoint found
    return HF_ENDPOINTS[0]


def download(
    url: str,
    local_path: Union[str, Path],
    overwrite: bool = False,
    chunk_size: int = 8192,
    timeout: int = 30,
    max_retries: int = 3,
    retry_delay: float = 1.0,
    headers: Optional[Dict[str, str]] = None,
    verify_ssl: bool = True,
    expected_checksum: Optional[str] = None,
    checksum_algorithm: str = "md5",
    show_progress: bool = True,
    create_dirs: bool = True,
) -> Path:
    """
    Download data from a URL and save to local path with progress tracking.

    Parameters
    ----------
    url : str
        URL to download data from
    local_path : Union[str, Path]
        Local path where the downloaded file will be saved
    overwrite : bool, default=False
        Whether to overwrite existing files
    chunk_size : int, default=8192
        Size of chunks to download at a time (in bytes)
    timeout : int, default=30
        Request timeout in seconds
    max_retries : int, default=3
        Maximum number of retry attempts if download fails
    retry_delay : float, default=1.0
        Delay between retry attempts in seconds
    headers : Optional[Dict[str, str]], default=None
        Additional HTTP headers to send with the request
    verify_ssl : bool, default=True
        Whether to verify SSL certificates
    expected_checksum : Optional[str], default=None
        Expected checksum of the downloaded file for verification
    checksum_algorithm : str, default="md5"
        Algorithm to use for checksum verification ("md5", "sha1", "sha256")
    show_progress : bool, default=True
        Whether to show download progress bar
    create_dirs : bool, default=True
        Whether to create parent directories if they don't exist

    Returns
    -------
    Path
        Path object pointing to the downloaded file

    Raises
    ------
    ValueError
        If URL is invalid or checksum verification fails
    FileExistsError
        If file exists and overwrite=False
    requests.RequestException
        If download fails after all retries

    Examples
    --------
    Basic usage:

    >>> file_path = download_origin_data(
    ...     "https://example.com/data.csv",
    ...     "data/raw_data.csv"
    ... )
    >>> print(f"Downloaded to: {file_path}")
    Downloaded to: data/raw_data.csv

    With checksum verification:

    >>> file_path = download_origin_data(
    ...     "https://example.com/important_data.xlsx",
    ...     "data/important_data.xlsx",
    ...     expected_checksum="5d41402abc4b2a76b9719d911017c592",
    ...     checksum_algorithm="md5"
    ... )

    With custom headers and retry settings:

    >>> headers = {"User-Agent": "MyApp/1.0"}
    >>> file_path = download_origin_data(
    ...     "https://api.example.com/dataset.json",
    ...     "data/dataset.json",
    ...     headers=headers,
    ...     max_retries=5,
    ...     retry_delay=2.0
    ... )

    Download without progress bar:

    >>> file_path = download_origin_data(
    ...     "https://example.com/data.tsv",
    ...     "data/data.tsv",
    ...     show_progress=False,
    ...     overwrite=True
    ... )
    """
    # Convert to Path object
    local_path = Path(local_path)

    # Validate URL
    parsed_url = urlparse(url)
    if not parsed_url.scheme or not parsed_url.netloc:
        raise ValueError(f"Invalid URL: {url}")

    tqdm.write(f"Downloading data from {url}...")
    tqdm.write(f"Target location: {local_path}")

    # Check if file exists
    if local_path.exists() and not overwrite:
        raise FileExistsError(
            f"File already exists: {local_path}. Use overwrite=True to replace it."
        )

    # Create parent directories if needed
    if create_dirs and local_path.parent != Path("."):
        local_path.parent.mkdir(parents=True, exist_ok=True)
        tqdm.write(f"Created directory: {local_path.parent}")

    # Prepare headers
    default_headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
        "(KHTML, like Gecko) Chrome/139.0.0.0 Safari/537.36"
    }
    if headers:
        default_headers.update(headers)

    # Download with retries
    for attempt in range(max_retries):
        try:
            tqdm.write(f"Download attempt {attempt + 1}/{max_retries}...")

            # Make initial request to get file size
            response = requests.head(
                url,
                headers=default_headers,
                timeout=timeout,
                verify=verify_ssl,
                allow_redirects=True,
            )
            response.raise_for_status()

            # Get content length for progress bar
            total_size = int(response.headers.get("content-length", 0))

            # Start actual download
            response = requests.get(
                url,
                headers=default_headers,
                timeout=timeout,
                verify=verify_ssl,
                stream=True,
                allow_redirects=True,
            )
            response.raise_for_status()

            # Initialize progress bar
            progress_bar = None
            if show_progress and total_size > 0:
                progress_bar = tqdm(
                    total=total_size,
                    unit="B",
                    unit_scale=True,
                    unit_divisor=1024,
                    desc=f"Downloading {local_path.name}",
                )

            # Initialize checksum if needed
            checksum_hash = None
            if expected_checksum:
                if checksum_algorithm.lower() == "md5":
                    checksum_hash = hashlib.md5()
                elif checksum_algorithm.lower() == "sha1":
                    checksum_hash = hashlib.sha1()
                elif checksum_algorithm.lower() == "sha256":
                    checksum_hash = hashlib.sha256()
                else:
                    raise ValueError(
                        f"Unsupported checksum algorithm: {checksum_algorithm}"
                    )

            # Download and save file
            with open(local_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=chunk_size):
                    if chunk:  # Filter out keep-alive chunks
                        f.write(chunk)
                        if checksum_hash:
                            checksum_hash.update(chunk)
                        if progress_bar:
                            progress_bar.update(len(chunk))

            if progress_bar:
                progress_bar.close()

            # Verify checksum if provided
            if expected_checksum and checksum_hash:
                calculated_checksum = checksum_hash.hexdigest()
                if calculated_checksum.lower() != expected_checksum.lower():
                    local_path.unlink()  # Remove corrupted file
                    raise ValueError(
                        f"Checksum verification failed. "
                        f"Expected: {expected_checksum}, "
                        f"Got: {calculated_checksum}"
                    )
                tqdm.write(
                    f"Checksum verification passed ({checksum_algorithm.upper()})"
                )

            file_size = local_path.stat().st_size
            tqdm.write(f"Successfully downloaded {file_size:,} bytes to {local_path}")
            return local_path

        except requests.RequestException as e:
            tqdm.write(f"Download attempt {attempt + 1} failed: {str(e)}")
            if attempt < max_retries - 1:
                tqdm.write(f"Retrying in {retry_delay} seconds...")
                time.sleep(retry_delay)
            else:
                tqdm.write("All download attempts failed.")
                raise

        except Exception as e:
            # Clean up partial download on unexpected errors
            if local_path.exists():
                local_path.unlink()
            raise
    # should never reached because already deal with exception in loop
    raise RuntimeError(f"Failed to download {url} after {max_retries} attempts.")


def download_source_file_from_huggingface(
    dataset_name: str,
    dir: str,
    *,
    overwrite: bool = False,
    sub_dataset: Optional[str] = None,
) -> Dict[str, str]:
    """
    Download the source file for a given dataset.

    This function retrieves the file URL from the `DATASETS` registry by dataset name,
    downloads it to the specified directory, and saves it under the given filename.

    All source files are downloaded from Hugging Face:
    https://huggingface.co/datasets/xulab-research/MutCleaner/tree/main

    Parameters
    ----------
    dataset_name : str
        The key identifying the dataset in the `DATASETS` registry.
    dir : str
        The target directory where the file will be saved.
    overwrite : bool, default=False
        Whether to overwrite the file if it already exists. Default is False.
    sub_dataset : Optional[str], default=None
        If provided, retrieves the file from the specified sub-dataset within the Hugging Face repository.

    Returns
    -------
    Dict[str]
        key: file name,
        value: The local filesystem paths to the downloaded files.

    Raises
    ------
    ValueError
        If no file URL is found for the given dataset.
    FileExistsError
        If the file already exists and `overwrite` is False.

    Examples
    --------
    Download all files registered for the Human Domainome Dataset:

    >>> paths = download_source_file_from_huggingface(
    ...     "Human Domainome Dataset",
    ...     "data",
    ... )
    >>> paths["SupplementaryTable2.txt"]
    'data/SupplementaryTable2.txt'
    >>> paths["SupplementaryTable4.txt"]
    'data/SupplementaryTable4.txt'
    >>> paths["wild_type.fasta"]
    'data/wild_type.fasta'

    Download files from the Human Domainome Sup2 Dataset:

    >>> paths = download_source_file_from_huggingface(
    ...     "Human Domainome Dataset",
    ...     "data",
    ...     sub_dataset="Human Domainome Sup2 Dataset",
    ... )
    >>> paths["SupplementaryTable2.txt"]
    'data/SupplementaryTable2.txt'

    Download files from the Human Domainome Sup4 Dataset:

    >>> paths = download_source_file_from_huggingface(
    ...     "Human Domainome Dataset",
    ...     "data",
    ...     sub_dataset="Human Domainome Sup4 Dataset",
    ... )
    >>> paths["SupplementaryTable4.txt"]
    'data/SupplementaryTable4.txt'
    >>> paths["wild_type.fasta"]
    'data/wild_type.fasta'
    """
    target_dataset = DATASETS.get(dataset_name, {})
    if sub_dataset is not None:
        target_dataset = target_dataset.get("sub_datasets", {}).get(sub_dataset, {})
    if not target_dataset:
        raise ValueError(f"No dataset found with name: {dataset_name}")
    hf_repos = target_dataset.get("huggingface_repos", [])
    if len(hf_repos) == 0:
        raise ValueError(
            f"No Hugging Face repository found for dataset: {dataset_name}"
        )
    filenames = target_dataset.get("file_name", [])

    local_paths = {}
    for repo, filename in zip(hf_repos, filenames):
        if filename is None:
            raise ValueError("No file name provided and no default found in dataset")

        url = f"{resolve_hf_endpoint()}{repo}"
        local_path = Path(dir) / filename
        if local_path.exists():
            if not overwrite:
                raise FileExistsError(
                    f"File already exists: {local_path}. Use overwrite=True to replace it."
                )
            else:
                local_path.unlink()
        print(url)
        download(url, local_path)
        local_paths[filename] = str(local_path)
    return local_paths


def download_protein_cdna_proteolysis_source_file(
    dir: str, *, overwrite: bool = False
) -> Dict[str, str]:
    """
    Download the source file for Protein cDNA Proteolysis Dataset from the original source.

    Returns
    -------
    Dict[str, str]
        key: file name,
        value: file path pointing to Protein cDNA Proteolysis Dataset source file
    """
    return download_source_file_from_huggingface(
        "Protein cDNA Proteolysis Dataset", dir, overwrite=overwrite
    )


def download_proteingym_source_file(
    dir: str, *, overwrite: bool = False
) -> Dict[str, str]:
    """
    Download the source file for ProteinGym DMS Substitutions Dataset from the original source.

    Returns
    -------
    Dict[str, str]
        key: file name,
        value: file path pointing to the ProteinGym DMS Substitutions Dataset source file
    """
    return download_source_file_from_huggingface("ProteinGym DMS Substitutions Dataset", dir, overwrite=overwrite)


def download_human_domainome_source_file(
    dir: str,
    *,
    overwrite: bool = False,
    sub_dataset: Optional[
        Literal[
            "Human Domainome Sup2 Dataset",
            "Human Domainome Sup4 Dataset",
        ]
    ] = None,
) -> Dict[str, str]:
    """
    Download source files for the Human Domainome Dataset.

    Parameters
    ----------
    dir : str
        The target directory where the files will be saved.
    overwrite : bool, default=False
        Whether to overwrite existing files.
    sub_dataset : Optional[
        Literal[
            "Human Domainome Sup2 Dataset",
            "Human Domainome Sup4 Dataset",
        ]
    ], default=None
        Sub-dataset to download. If None, download all Human Domainome
        source files.

    Returns
    -------
    Dict[str, str]
        Mapping from file names to downloaded local file paths.
    """
    return download_source_file_from_huggingface(
        "Human Domainome Dataset",
        dir,
        overwrite=overwrite,
        sub_dataset=sub_dataset,
    )


def download_ddg_dtm_source_file(
    dir: str,
    *,
    dataset_type: Literal["ddg", "dtm"],
    overwrite: bool = False,
    sub_dataset: Optional[str] = None,
) -> Dict[str, str]:
    """
    Download the source file for ddG or dTm datasets from the original source.

    Parameters
    ----------
    dir : str
        The target directory where the file will be saved.
    dataset_type : Literal["ddg", "dtm"]
        Dataset type to download.
        - "ddg": download ΔΔG datasets
        - "dtm": download ΔTm datasets
    overwrite : bool, default=False
        Whether to overwrite the file if it already exists. Default is False.
    sub_dataset : Optional[str], default=None
        Sub-dataset to download. If None, download all datasets of the selected type.

        Supported options for dataset_type="ddg":
        "M1261", "S461", "S669", "S783", "S8754"

        Supported options for dataset_type="dtm":
        "S4346", "S557"

    Returns
    -------
    Dict[str, str]
        key: file name,
        value: file path pointing to the ddG or dTm dataset source file
    """
    dataset_map = {
        "ddg": "ΔΔG Dataset",
        "dtm": "ΔTm Dataset",
    }

    if dataset_type not in dataset_map:
        raise ValueError(
            f"Unsupported dataset type: {dataset_type}. "
            f'Supported options: "ddg", "dtm".'
        )

    dataset_name = dataset_map[dataset_type]
    supported_sub_datasets = list(DATASETS[dataset_name]["sub_datasets"].keys())

    if sub_dataset is not None and sub_dataset not in supported_sub_datasets:
        raise ValueError(
            f"Unsupported sub-dataset for {dataset_type}. Supported options: "
            f"{', '.join(supported_sub_datasets)}"
        )

    return download_source_file_from_huggingface(
        dataset_name,
        dir,
        overwrite=overwrite,
        sub_dataset=sub_dataset,
    )


def download_archstabms1e10_source_file(
    dir: str, *, overwrite: bool = False, sub_dataset: Optional[str] = None
) -> Dict[str, str]:
    """
    Download the source files for the ArchStabMS1E10 dataset.

    Parameters
    ----------
    dir : str
        The target directory where the files will be saved.
    overwrite : bool, default=False
        Whether to overwrite existing files.
    sub_dataset : Optional[str], default=None
        Name of the sub-dataset to download. If None, all registered
        ArchStabMS1E10 source files are downloaded.

    Returns
    -------
    Dict[str, str]
        Mapping from file names to downloaded local file paths.
    """
    return download_source_file_from_huggingface(
        "ArchStabMS1E10 Epistasis Dataset", dir, overwrite=overwrite, sub_dataset=sub_dataset
    )


def download_human_myoglobin_source_file(
    dir: str, *, overwrite: bool = False
) -> Dict[str, str]:
    """
    Download the source file for HumanMyoglobin dataset from the original source.

    Parameters
    ----------
    dir : str
        The target directory where the file will be saved
    overwrite : bool, default=False
        Whether to overwrite the file if it already exists. Default is False.

    Returns
    -------
    Dict[str, str]
        key: file name,
        value: file path pointing to HumanMyoglobin dataset source file
    """
    return download_source_file_from_huggingface(
        "Human Myoglobin Epistasis Dataset", dir, overwrite=overwrite, sub_dataset="Human Myoglobin Epistasis Dataset"
    )


def download_ctxm_source_file(
    dir: str,
    *, 
    overwrite: bool = False,
    sub_dataset: Optional[Literal["CTXM_Ampicillin_Epistasis_Dataset", "CTXM_Cefotaxime_Epistasis_Dataset"]] = None, 
) -> Dict[str, str]:
    """
    Download the source file for CTX-M dataset from the original source.

    Parameters
    ----------
    dir : str
        The target directory where the file will be saved
    overwrite : bool, default=False
        Whether to overwrite the file if it already exists. Default is False.

    Returns
    -------
    Dict[str, str]
        key: file name,
        value: file path pointing to CTX-M dataset source file
    """
    if sub_dataset is not None and sub_dataset not in ["CTXM_Ampicillin_Epistasis_Dataset", "CTXM_Cefotaxime_Epistasis_Dataset"]:
        raise ValueError("Unsupported sub-dataset. Supported options: CTXM_Ampicillin_Epistasis_Dataset, CTXM_Cefotaxime_Epistasis_Dataset")
    return download_source_file_from_huggingface("CTXM Epistasis Dataset", dir, overwrite=overwrite, sub_dataset=sub_dataset)


def download_trpb_source_file(
    dir: str, *, overwrite: bool = False
) -> Dict[str, str]:
    """
    Download the source file for TrpB dataset from the original source.

    Parameters
    ----------
    dir : str
        The target directory where the file will be saved
    overwrite : bool, default=False
        Whether to overwrite the file if it already exists. Default is False.

    Returns
    -------
    Dict[str, str]
        key: file name,
        value: file path pointing to TrpB dataset source file
    """
    return download_source_file_from_huggingface("TrpB Epistasis Dataset", dir, overwrite=overwrite)


def download_antitoxin_pard3_source_file(
    dir: str, *, overwrite: bool = False
) -> Dict[str, str]:
    """
    Download the source file for Antitoxin_ParD3 dataset from the original source.

    Parameters
    ----------
    dir : str
        The target directory where the file will be saved
    overwrite : bool, default=False
        Whether to overwrite the file if it already exists. Default is False.

    Returns
    -------
    Dict[str, str]
        key: file name,
        value: file path pointing to Antitoxin ParD3 Epistasis Dataset source file
    """
    return download_source_file_from_huggingface("Antitoxin ParD3 Epistasis Dataset", dir, overwrite=overwrite)


def download_rbd_antibody_source_file(
    dir: str,
    *,
    overwrite: bool = False,
    sub_dataset: Optional[
        Literal[
            "Moderna",
            "Rockefeller",
            "Vir_mAbs",
        ]
    ] = None,
) -> Dict[str, str]:
    """
    Download the source file(s) for the curated RBD_Antibody dataset bundle.

    Parameters
    ----------
    dir : str
        The target directory where the file(s) will be saved.
    overwrite : bool, default=False
        Whether to overwrite existing files.
    sub_dataset : Optional[Literal[...]], default=None
        If provided, download only the specified RBD antibody sub-dataset.

    Returns
    -------
    Dict[str, str]
        key: file name,
        value: local file path
    """
    valid_sub_datasets = ["Moderna", "Rockefeller", "Vir_mAbs"]
    if sub_dataset is not None and sub_dataset not in valid_sub_datasets:
        raise ValueError(
            f"Unsupported sub-dataset. Supported options: {valid_sub_datasets}"
        )

    return download_source_file_from_huggingface(
        "RBD Antibody Dataset",
        dir,
        overwrite=overwrite,
        sub_dataset=sub_dataset,
    )


def download_rbd_ace2_source_file(
    dir: str,
    *,
    overwrite: bool = False,
    sub_dataset: Optional[
        Literal[
            "Omicron_EG5_FLip_BA286",
            "Omicron_XBB_BQ",
            "Omicron",
            "DMS_variants",
            "Delta",
        ]
    ] = None,
) -> Dict[str, str]:
    """
    Download the source file(s) for the RBD_ACE2 dataset from Hugging Face.

    Parameters
    ----------
    dir : str
        The target directory where the file(s) will be saved.
    overwrite : bool, default=False
        Whether to overwrite existing files.
    sub_dataset : Optional[Literal[...]], default=None
        If provided, download only the specified RBD ACE2 sub-dataset.

    Returns
    -------
    Dict[str, str]
        key: file name,
        value: local file path
    """
    valid_sub_datasets = [
        "Omicron_EG5_FLip_BA286",
        "Omicron_XBB_BQ",
        "Omicron",
        "DMS_variants",
        "Delta",
    ]
    if sub_dataset is not None and sub_dataset not in valid_sub_datasets:
        raise ValueError(
            f"Unsupported sub-dataset. Supported options: {valid_sub_datasets}"
        )

    return download_source_file_from_huggingface(
        "RBD ACE2 Dataset",
        dir,
        overwrite=overwrite,
        sub_dataset=sub_dataset,
    )


def download_chitosanase_dtm_source_file(dir: str, *, overwrite: bool = False) -> Dict[str, str]:
    """
    Download the source file for Chitosanase dTm Dataset from the original source.

    Parameters
    ----------
    dir : str
        The target directory where the file will be saved
    overwrite : bool, default=False
        Whether to overwrite the file if it already exists. Default is False.

    Returns
    -------
    Dict[str, str]
        key: file name,
        value: file path pointing to Chitosanase dTm Dataset source file
    """
    return download_source_file_from_huggingface("Chitosanase dTm Dataset", dir, overwrite=overwrite)


def download_mgnify_ddg_source_file(dir: str, *, overwrite: bool = False) -> Dict[str, str]:
    """
    Download the source file for MGnify ddG Dataset from the original source.

    Parameters
    ----------
    dir : str
        The target directory where the file will be saved
    overwrite : bool, default=False
        Whether to overwrite the file if it already exists. Default is False.

    Returns
    -------
    Dict[str, str]
        key: file name,
        value: file path pointing to MGnify ddG Dataset source file
    """
    return download_source_file_from_huggingface("MGnify ddG Dataset", dir, overwrite=overwrite)
