# Supported Datasets

## Human Domainome Dataset

### Basic Usage
The following example shows the complete workflow for downloading, cleaning,
saving the cleaned Human Domainome dataset, and exporting the
cleaning artifacts: 
```python
import pickle
from pathlib import Path
from mutcleaner import download_human_domainome_source_file
from mutcleaner.cleaners import (
    create_human_domainome_sup2_cleaner,
    clean_human_domainome_sup2_dataset,
)


def main():
    # Prepare data
    download_human_domainome_source_file("raw_dataset/Human_Domainome_Dataset")

    # File settings
    dataset_file_path = Path(
        "raw_dataset/Human_Domainome_Dataset/SupplementaryTable2.txt"
    )
    artifact_path = Path("logs/Human_Domainome_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/Human_Domainome_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    hd_cleaning_pipeline = create_human_domainome_sup2_cleaner(
        dataset_file_path
    )
    hd_cleaning_pipeline, hd_dataset = clean_human_domainome_sup2_dataset(
        hd_cleaning_pipeline
    )

    # Save data
    hd_dataset.save("cleaned_dataset/cleaned_Human_Domainome_Dataset")
    hd_cleaning_pipeline.save_artifacts(artifact_path)

    # Read artifacts from the pickle file and read the object
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(
            f"{artifact_csv_dir}/{artifact_name}.csv", index=False
        )


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()
```

### Advanced Settings

See {py:class}`mutcleaner.cleaners.HumanDomainomeSup2CleanerConfig` for details.

## ProteinGym DMS Substitutions Dataset

### Basic Usage
The following example shows the complete workflow for downloading, cleaning,
saving the cleaned ProteinGym DMS substitutions dataset, and exporting the
cleaning artifacts:
```python
import pickle
from pathlib import Path
from mutcleaner import download_proteingym_source_file
from mutcleaner.cleaners import (
    create_proteingym_dms_substitutions_cleaner,
    clean_proteingym_dms_substitutions_dataset,
)


def main():
    # Prepare data
    download_proteingym_source_file(
        "raw_dataset/ProteinGym_DMS_Substitutions_Dataset"
    )

    # File settings
    dataset_file_path = Path(
        "raw_dataset/ProteinGym_DMS_Substitutions_Dataset/ProteinGym_DMS_substitutions.zip"
    )
    artifact_path = Path(
        "logs/ProteinGym_DMS_Substitutions_Dataset/artifacts.pkl"
    )
    artifact_csv_dir = Path("logs/ProteinGym_DMS_Substitutions_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    pg_cleaning_pipeline = create_proteingym_dms_substitutions_cleaner(
        dataset_file_path
    )
    pg_cleaning_pipeline, pg_dataset = (
        clean_proteingym_dms_substitutions_dataset(pg_cleaning_pipeline)
    )

    # Save data
    pg_dataset.save(
        "cleaned_dataset/cleaned_ProteinGym_DMS_Substitutions_Dataset"
    )
    pg_cleaning_pipeline.save_artifacts(artifact_path)

    # Read artifacts from the pickle file and read the object
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(
            f"{artifact_csv_dir}/{artifact_name}.csv", index=False
        )


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()
```

### Advanced Settings

See {py:class}`mutcleaner.cleaners.ProteinGymCleanerConfig` for details.

## cDNA Proteolysis Dataset


### ΔΔG as Label (Default Pipeline)

The following example shows the complete workflow for downloading, cleaning,
saving the cleaned cDNA Proteolysis dataset, and exporting the
cleaning artifacts, using ΔΔG as the label:
```python
import pickle
from pathlib import Path
from mutcleaner import download_cdna_proteolysis_source_file
from mutcleaner.cleaners import (
    create_cdna_proteolysis_cleaner,
    clean_cdna_proteolysis_dataset,
)


def main():
    # Prepare dataset
    download_cdna_proteolysis_source_file(
        "raw_dataset/cDNA_Proteolysis_Dataset"
    )

    # File settings
    dataset_file_path = Path(
        "raw_dataset/cDNA_Proteolysis_Dataset/Tsuboyama2023_Dataset2_Dataset3_20230416.csv"
    )
    artifact_path = Path("logs/cDNA_Proteolysis_ddG_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/cDNA_Proteolysis_ddG_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    cdna_cleaning_pipeline = create_cdna_proteolysis_cleaner(dataset_file_path)
    cdna_cleaning_pipeline, cdna_dataset = clean_cdna_proteolysis_dataset(
        cdna_cleaning_pipeline
    )

    # Save data
    cdna_dataset.save("cleaned_dataset/cleaned_cDNA_Proteolysis_ddG_Dataset")
    cdna_cleaning_pipeline.save_artifacts(artifact_path)

    # Read artifacts from the pickle file
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(
            f"{artifact_csv_dir}/{artifact_name}.csv", index=False
        )


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()
```

### ΔG as Label

The following example shows the complete workflow for downloading, cleaning,
saving the cleaned cDNA Proteolysis dataset, and exporting the
cleaning artifacts, using ΔG as the label:
```python
import pickle
from pathlib import Path
from mutcleaner.cleaners import (
    CDNAProteolysisCleanerConfig,
    create_cdna_proteolysis_cleaner,
    clean_cdna_proteolysis_dataset,
)


def main():
    # Set cleaning configs
    cdna_cleaning_config = CDNAProteolysisCleanerConfig()
    cdna_cleaning_config.column_mapping = {
        "WT_name": "name",
        "aa_seq": "mut_seq",
        "mut_type": "mut_info",
        "dG_ML": "label_cDNAProteolysis",
    }

    # File settings
    dataset_file_path = Path(
        "raw_dataset/cDNA_Proteolysis_Dataset/Tsuboyama2023_Dataset2_Dataset3_20230416.csv"
    )
    artifact_path = Path("logs/cDNA_Proteolysis_dG_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/cDNA_Proteolysis_dG_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    cdna_cleaning_pipeline = create_cdna_proteolysis_cleaner(
        dataset_file_path, cdna_cleaning_config
    )
    cdna_cleaning_pipeline, cdna_dataset = clean_cdna_proteolysis_dataset(
        cdna_cleaning_pipeline
    )

    # Save data
    cdna_dataset.save("cleaned_dataset/cleaned_cDNA_Proteolysis_dG_Dataset")
    cdna_cleaning_pipeline.save_artifacts(artifact_path)

    # Read artifacts from the pickle file
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(
            f"{artifact_csv_dir}/{artifact_name}.csv", index=False
        )


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()
```

### Advanced Settings
See {py:class}`mutcleaner.cleaners.CDNAProteolysisCleanerConfig` for details.

## ddG Dataset

### Basic Usage

The following example shows the complete workflow for downloading, cleaning,
saving the cleaned ΔΔG dataset, and exporting the cleaning artifacts:
```python
import pickle
from pathlib import Path
from mutcleaner import download_ddg_dtm_source_file
from mutcleaner.cleaners import (
    create_ddg_dtm_cleaner, 
    clean_ddg_dtm_dataset,
)


def main():
    # Prepare data
    raw_data_dir = Path("raw_dataset/ddG_Dataset")
    download_ddg_dtm_source_file(raw_data_dir, dataset_type="ddg", overwrite=True)

    # File settings
    for dataset_filepath in sorted(raw_data_dir.glob("*.csv")):
        data_file = dataset_filepath.stem
        artifact_path = Path(f"logs/ddG_Dataset/{data_file}/artifacts.pkl")
        artifact_csv_dir = Path(f"logs/ddG_Dataset/{data_file}")

        artifact_csv_dir.mkdir(parents=True, exist_ok=True)

        # Clean data
        ddgdtm_cleaning_pipeline = create_ddg_dtm_cleaner(dataset_filepath)
        ddgdtm_cleaning_pipeline, ddgdtm_dataset = clean_ddg_dtm_dataset(
            ddgdtm_cleaning_pipeline
        )

        # Save data
        ddgdtm_dataset.save(f"cleaned_dataset/cleaned_ddG_Dataset/{data_file}")
        ddgdtm_cleaning_pipeline.save_artifacts(artifact_path)

        # open the pickle file
        with open(artifact_path, "rb") as file:
            artifacts = pickle.load(file)

        for artifact_name, artifact_df in artifacts.items():
            artifact_df.to_csv(artifact_csv_dir / f"{artifact_name}.csv", index=False)


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()
```

### Advanced Settings

See {py:class}`mutcleaner.cleaners.DdgDtmCleanerConfig` for details.

## dTm Dataset

### Basic Usage

The following example shows the complete workflow for downloading, cleaning,
saving the cleaned ΔTm dataset, and exporting the cleaning artifacts:
```python
import pickle
from pathlib import Path
from mutcleaner import download_ddg_dtm_source_file
from mutcleaner.cleaners import (
    create_ddg_dtm_cleaner, 
    clean_ddg_dtm_dataset,
)


def main():
    # Prepare data
    raw_data_dir = Path("raw_dataset/dTm_Dataset")
    download_ddg_dtm_source_file(raw_data_dir, dataset_type="dtm", overwrite=True)

    # File settings
    for dataset_filepath in sorted(raw_data_dir.glob("*.csv")):
        data_file = dataset_filepath.stem
        artifact_path = Path(f"logs/dTm_Dataset/{data_file}/artifacts.pkl")
        artifact_csv_dir = Path(f"logs/dTm_Dataset/{data_file}")

        artifact_csv_dir.mkdir(parents=True, exist_ok=True)

        # Clean data
        ddgdtm_cleaning_pipeline = create_ddg_dtm_cleaner(dataset_filepath)
        ddgdtm_cleaning_pipeline, ddgdtm_dataset = clean_ddg_dtm_dataset(
            ddgdtm_cleaning_pipeline
        )

        # Save data
        ddgdtm_dataset.save(f"cleaned_dataset/cleaned_dTm_Dataset/{data_file}")
        ddgdtm_cleaning_pipeline.save_artifacts(artifact_path)

        # open the pickle file
        with open(artifact_path, "rb") as file:
            artifacts = pickle.load(file)

        for artifact_name, artifact_df in artifacts.items():
            artifact_df.to_csv(artifact_csv_dir / f"{artifact_name}.csv", index=False)


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()
```

### Advanced Settings

See {py:class}`mutcleaner.cleaners.DdgDtmCleanerConfig` for details.

## ArchStabMS1E10 Epistasis Dataset

### Basic Usage

The following example shows the complete workflow for downloading, cleaning,
saving the cleaned ArchStabMS1E10 Epistasis dataset, and exporting the
cleaning artifacts:
```python
import pickle
from pathlib import Path
from mutcleaner import download_archstabms1e10_source_file
from mutcleaner.cleaners import (
    create_archstabms_1e10_cleaner,
    clean_archstabms_1e10_dataset,
)


def main():
    # Prepare dataset
    download_archstabms1e10_source_file(
        "raw_dataset/ArchStabMS1E10_Epistasis_Dataset"
    )

    # File settings
    dataset_file_path = Path(
        "raw_dataset/ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Dataset.csv"
    )
    artifact_path = Path("logs/ArchStabMS1E10_Epistasis_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/ArchStabMS1E10_Epistasis_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    archstabms_cleaning_pipeline = create_archstabms_1e10_cleaner(
        dataset_file_path
    )
    archstabms_cleaning_pipeline, archstabms_dataset = (
        clean_archstabms_1e10_dataset(archstabms_cleaning_pipeline)
    )

    # Save data
    archstabms_dataset.save(
        "cleaned_dataset/cleaned_ArchStabMS1E10_Epistasis_Dataset"
    )
    archstabms_cleaning_pipeline.save_artifacts(artifact_path)

    # Read artifacts from the pickle file
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(
            f"{artifact_csv_dir}/{artifact_name}.csv", index=False
        )


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()
```

### Advanced Settings

See {py:class}`mutcleaner.cleaners.ArchStabMS1E10CleanerConfig` for details.

## Antitoxin ParD3 Epistasis Dataset

### Basic Usage

The following example shows the complete workflow for downloading, cleaning,
saving the cleaned Antitoxin ParD3 Epistasis dataset, and exporting the
cleaning artifacts:
```python
import pickle
from pathlib import Path
from mutcleaner import download_antitoxin_pard3_source_file
from mutcleaner.cleaners import (
    create_antitoxin_pard3_cleaner,
    clean_antitoxin_pard3_dataset,
)


def main():
    # Prepare dataset
    download_antitoxin_pard3_source_file(
        "raw_dataset/Antitoxin_ParD3_Epistasis_Dataset"
    )

    # File settings
    dataset_file_path = Path(
        "raw_dataset/Antitoxin_ParD3_Epistasis_Dataset/Antitoxin_ParD3_Epistasis_Dataset.csv"
    )
    artifact_path = Path("logs/Antitoxin_ParD3_Epistasis_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/Antitoxin_ParD3_Epistasis_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    antitoxin_pard3_cleaning_pipeline = create_antitoxin_pard3_cleaner(
        dataset_file_path
    )
    antitoxin_pard3_cleaning_pipeline, antitoxin_pard3_dataset = (
        clean_antitoxin_pard3_dataset(antitoxin_pard3_cleaning_pipeline)
    )

    # Save data
    antitoxin_pard3_dataset.save(
        "cleaned_dataset/cleaned_Antitoxin_ParD3_Epistasis_Dataset"
    )
    antitoxin_pard3_cleaning_pipeline.save_artifacts(artifact_path)

    # Read artifacts from the pickle file
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(
            f"{artifact_csv_dir}/{artifact_name}.csv", index=False
        )


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()
```

### Advanced Settings

See {py:class}`mutcleaner.cleaners.AntitoxinParD3CleanerConfig` for details.

## TrpB Epistasis Dataset

### Basic Usage

The following example shows the complete workflow for downloading, cleaning,
saving the cleaned TrpB Epistasis dataset, and exporting the
cleaning artifacts:
```python
import pickle
from pathlib import Path
from mutcleaner import download_trpb_source_file
from mutcleaner.cleaners import (
    create_trpb_cleaner,
    clean_trpb_dataset,
)


def main():
    # Prepare data
    download_trpb_source_file("raw_dataset/TrpB_Epistasis_Dataset")

    # File settings
    dataset_file_path = Path(
        "raw_dataset/TrpB_Epistasis_Dataset/TrpB_Epistasis_Dataset.csv"
    )
    artifact_path = Path("logs/TrpB_Epistasis_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/TrpB_Epistasis_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    trpb_cleaning_pipeline = create_trpb_cleaner(dataset_file_path)
    trpb_cleaning_pipeline, trpb_dataset = clean_trpb_dataset(
        trpb_cleaning_pipeline
    )

    # Save data
    trpb_dataset.save("cleaned_dataset/cleaned_TrpB_Epistasis_Dataset")
    trpb_cleaning_pipeline.save_artifacts(artifact_path)

    # Read artifacts from the pickle file and read the object
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(
            f"{artifact_csv_dir}/{artifact_name}.csv", index=False
        )


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()
```

### Advanced Settings

See {py:class}`mutcleaner.cleaners.TrpBCleanerConfig` for details.

## Human Myoglobin Epistasis Dataset

### Basic Usage

The following example shows the complete workflow for downloading, cleaning,
saving the cleaned Human Myoglobin Epistasis dataset, and exporting the
cleaning artifacts:
```python
import pickle
from pathlib import Path
from mutcleaner import download_human_myoglobin_source_file
from mutcleaner.cleaners import (
    create_human_myoglobin_cleaner,
    clean_human_myoglobin_dataset,
)


def main():
    # Prepare data
    download_human_myoglobin_source_file(
        "raw_dataset/Human_Myoglobin_Epistasis_Dataset"
    )

    # File settings
    dataset_file_path = Path(
        "raw_dataset/Human_Myoglobin_Epistasis_Dataset/Human_Myoglobin_Epistasis_Dataset.csv"
    )
    artifact_path = Path("logs/Human_Myoglobin_Epistasis_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/Human_Myoglobin_Epistasis_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    human_myoglobin_cleaning_pipeline = create_human_myoglobin_cleaner(
        dataset_file_path
    )
    human_myoglobin_cleaning_pipeline, human_myoglobin_dataset = (
        clean_human_myoglobin_dataset(human_myoglobin_cleaning_pipeline)
    )

    # Save data
    human_myoglobin_dataset.save(
        "cleaned_dataset/cleaned_Human_Myoglobin_Epistasis_Dataset"
    )
    human_myoglobin_cleaning_pipeline.save_artifacts(artifact_path)

    # Read artifacts from the pickle file and read the object
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(
            f"{artifact_csv_dir}/{artifact_name}.csv", index=False
        )


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()
```

### Advanced Settings

See {py:class}`mutcleaner.cleaners.HumanMyoglobinCleanerConfig` for details.

## CTXM Epistasis Dataset


### CTXM Ampicillin Epistasis Dataset

The following example shows the complete workflow for downloading, cleaning,
saving the cleaned CTXM Ampicillin Epistasis Dataset, and exporting the
cleaning artifacts:
```python
import pickle
from pathlib import Path
from mutcleaner import download_ctxm_source_file
from mutcleaner.cleaners import (
    create_ctxm_cleaner,
    clean_ctxm_dataset,
)


def main():
    # Prepare data
    download_ctxm_source_file("raw_dataset/CTXM_Epistasis_Dataset")

    # File settings
    dataset_file_path = Path(
        "raw_dataset/CTXM_Epistasis_Dataset/CTXM_Ampicillin_Epistasis_Dataset.csv"
    )
    artifact_path = Path("logs/CTXM_Ampicillin_Epistasis_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/CTXM_Ampicillin_Epistasis_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    ctxm_cleaning_pipeline = create_ctxm_cleaner(dataset_file_path)
    ctxm_cleaning_pipeline, ctxm_dataset = clean_ctxm_dataset(
        ctxm_cleaning_pipeline
    )

    # Save data
    ctxm_dataset.save("cleaned_dataset/cleaned_CTXM_Epistasis_Dataset")
    ctxm_cleaning_pipeline.save_artifacts(artifact_path)

    # Read artifacts from the pickle file
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(
            f"{artifact_csv_dir}/{artifact_name}.csv", index=False
        )


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()
```

### CTXM Cefotaxime Epistasis Dataset

The following example shows the complete workflow for downloading, cleaning,
saving the cleaned CTXM Cefotaxime Epistasis Dataset, and exporting the
cleaning artifacts:
```python
import pickle
from pathlib import Path
from mutcleaner.cleaners import (
    CTXMCleanerConfig,
    create_ctxm_cleaner,
    clean_ctxm_dataset,
)


def main():
    # Set cleaning configs
    ctxm_cleaning_config = CTXMCleanerConfig()
    ctxm_cleaning_config.wt_name = "CTXM_cefotaxime"

    # File settings
    dataset_file_path = Path(
        "raw_dataset/CTXM_Epistasis_Dataset/CTXM_Cefotaxime_Epistasis_Dataset.csv"
    )
    artifact_path = Path("logs/CTXM_Cefotaxime_Epistasis_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/CTXM_Cefotaxime_Epistasis_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    ctxm_cleaning_pipeline = create_ctxm_cleaner(
        dataset_file_path, ctxm_cleaning_config
    )
    ctxm_cleaning_pipeline, ctxm_dataset = clean_ctxm_dataset(
        ctxm_cleaning_pipeline
    )

    # Save data
    ctxm_dataset.save("cleaned_dataset/cleaned_CTXM_Epistasis_Dataset")
    ctxm_cleaning_pipeline.save_artifacts(artifact_path)

    # Read artifacts from the pickle file
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(
            f"{artifact_csv_dir}/{artifact_name}.csv", index=False
        )


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()
```

### Advanced Settings

See {py:class}`mutcleaner.cleaners.CTXMCleanerConfig` for details.


## ACE2-RBD Binding Scores Dataset

### Basic Usage

The following example shows the complete workflow for downloading, cleaning,
saving the cleaned ACE2-RBD binding scores dataset, and exporting the
cleaning artifacts from a standalone script in this repository layout:
```python
import pickle
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
OTHER_SOFTWARE_ROOT = REPO_ROOT / "other_software"
if str(OTHER_SOFTWARE_ROOT) not in sys.path:
    sys.path.insert(0, str(OTHER_SOFTWARE_ROOT))

from mutcleaner import download_rbd_ace2_source_file
from mutcleaner.cleaners import (
    create_ace2_binding_scores_cleaner,
    clean_ace2_binding_scores_dataset,
)


def main():
    # Prepare data
    raw_data_dir = Path("raw_dataset/RBD_ACE2_Dataset")
    download_rbd_ace2_source_file(str(raw_data_dir))

    # File settings
    dataset_file_path = (
        raw_data_dir
        / "SARS-CoV-2-RBD_DMS_Omicron-EG5-FLip-BA286_bc_binding.csv"
    )
    artifact_path = Path("logs/RBD_ACE2_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/RBD_ACE2_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    ace2_cleaning_pipeline = create_ace2_binding_scores_cleaner(
        dataset_file_path
    )
    ace2_cleaning_pipeline, ace2_dataset = (
        clean_ace2_binding_scores_dataset(ace2_cleaning_pipeline)
    )

    # Save data
    ace2_dataset.save("cleaned_dataset/cleaned_RBD_ACE2_Dataset")
    ace2_cleaning_pipeline.save_artifacts(artifact_path)

    # Read artifacts from the pickle file and read the object
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(
            f"{artifact_csv_dir}/{artifact_name}.csv", index=False
        )


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()
```

You can also download all source files at once:

```python
from mutcleaner import download_rbd_ace2_source_file

file_paths = download_rbd_ace2_source_file("path/to/target/folder")
```

Supported sub-datasets:
- `Omicron_EG5_FLip_BA286`
- `Omicron_XBB_BQ`
- `Omicron`
- `DMS_variants`
- `Delta`

Alternatively, you can download it from [Hugging Face](https://huggingface.co/datasets/xulab-research/RBD_ACE2).


### Advanced Settings

See {py:class}`mutcleaner.cleaners.ACE2BindingScoresCleanerConfig` for details.

## RBD-Antibody Dataset

### Basic Usage

You can download the source file directly by running (see {py:func}`mutcleaner.utils.download_rbd_antibody_source_file` for details):
```python
from mutcleaner import download_rbd_antibody_source_file

file_paths = download_rbd_antibody_source_file("path/to/target/folder")
```

You can also download and process a specific sub-dataset:

```python
from mutcleaner import download_rbd_antibody_source_file

file_paths = download_rbd_antibody_source_file(
    "path/to/target/folder",
    sub_dataset="AZ_Abs",
)
```

Supported sub-datasets:
- `AZ_Abs`
- `HAARVI_sera`
- `Moderna`
- `Rockefeller`
- `Vir_mAbs`
- `clinical_Abs`

Alternatively, you can download it from [Hugging Face](https://huggingface.co/datasets/Zoey13891350636/RBD_Antibody).



### Advanced Settings

See {py:class}`mutcleaner.cleaners.RBDAntibodyCleanerConfig` for details.
