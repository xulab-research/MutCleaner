# Supported Datasets

## Human Domainome Dataset

### Human Domainome Sup2 Dataset
The following example shows the complete workflow for downloading, cleaning,
saving the cleaned Human Domainome sup2 dataset, and exporting the
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
    download_human_domainome_source_file("raw_dataset/Human_Domainome_Dataset", overwrite=True)
    # File settings
    dataset_filepath = Path("raw_dataset/Human_Domainome_Dataset/SupplementaryTable2.txt")
    artifact_path = Path("logs/Human_Domainome_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/Human_Domainome_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    hd_cleaning_pipeline = create_human_domainome_sup2_cleaner(dataset_filepath)
    hd_cleaning_pipeline, hd_dataset = clean_human_domainome_sup2_dataset(
        hd_cleaning_pipeline
    )

    # Save data
    hd_dataset.save("cleaned_dataset/cleaned_Human_Domainome_Dataset")
    hd_cleaning_pipeline.save_artifacts(artifact_path)

    # open the pickle file and read the object
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(f"{artifact_csv_dir}/{artifact_name}.csv", index=False)


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()
```

### Advanced Settings

See {py:class}`mutcleaner.cleaners.HumanDomainomeSup2CleanerConfig` for details.

### Human Domainome Sup4 Dataset
The following example shows the complete workflow for downloading, cleaning,
saving the cleaned Human Domainome sup4 dataset, and exporting the cleaning artifacts: 
```python
import pickle
from pathlib import Path
from mutcleaner import download_human_domainome_source_file
from mutcleaner.cleaners import (
    create_human_domainome_sup4_cleaner,
    clean_human_domainome_sup4_dataset,
)


def main():
    # File settings
    download_human_domainome_source_file("raw_dataset/Human_Domainome_Dataset", overwrite=True, sub_dataset="Human Domainome Sup4 Dataset")
    artifact_path = Path("logs/Human_Domainome_Dataset/Human_Domainome_Sup4_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/Human_Domainome_Dataset/Human_Domainome_Sup4_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    hd_cleaning_pipeline = create_human_domainome_sup4_cleaner(dataset_or_path="raw_dataset/Human_Domainome_Dataset/SupplementaryTable4.txt", sequence_dict_path="raw_dataset/Human_Domainome_Dataset/wild_type.fasta")
    hd_cleaning_pipeline, hd_dataset = clean_human_domainome_sup4_dataset(
        hd_cleaning_pipeline
    )

    # Save data
    hd_dataset.save("cleaned_dataset/cleaned_Human_Domainome_Dataset/Human_Domainome_Sup4_Dataset")
    hd_cleaning_pipeline.save_artifacts(artifact_path)

    # open the pickle file and read the object
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(f"{artifact_csv_dir}/{artifact_name}.csv", index=False)


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()
```

### Advanced Settings

See {py:class}`mutcleaner.cleaners.HumanDomainomeSup4CleanerConfig` for details.

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
    download_proteingym_source_file("raw_dataset/ProteinGym_DMS_Substitutions_Dataset", overwrite=True)
    # File settings
    dataset_filepath = Path("raw_dataset/ProteinGym_DMS_Substitutions_Dataset/ProteinGym_DMS_substitutions.zip")
    artifact_path = Path("logs/ProteinGym_DMS_Substitutions_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/ProteinGym_DMS_Substitutions_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    pg_cleaning_pipeline = create_proteingym_dms_substitutions_cleaner(dataset_filepath)
    pg_cleaning_pipeline, pg_dataset = clean_proteingym_dms_substitutions_dataset(
        pg_cleaning_pipeline
    )

    # Save data
    pg_dataset.save("cleaned_dataset/cleaned_ProteinGym_DMS_Substitutions_Dataset")
    pg_cleaning_pipeline.save_artifacts(artifact_path)

    # open the pickle file and read the object
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(f"{artifact_csv_dir}/{artifact_name}.csv", index=False)


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
    download_cdna_proteolysis_source_file("raw_dataset/Protein_cDNA_Proteolysis_Dataset", overwrite=True)

    # File settings
    dataset_filepath = Path("raw_dataset/Protein_cDNA_Proteolysis_Dataset/Tsuboyama2023_Dataset2_Dataset3_20230416.csv")
    artifact_path = Path("logs/Protein_cDNA_Proteolysis_ddG_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/Protein_cDNA_Proteolysis_ddG_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    cdna_cleaning_pipeline = create_cdna_proteolysis_cleaner(dataset_filepath)
    cdna_cleaning_pipeline, cdna_dataset = clean_cdna_proteolysis_dataset(
        cdna_cleaning_pipeline
    )

    # Save data
    cdna_dataset.save("cleaned_dataset/cleaned_Protein_cDNA_Proteolysis_ddG_Dataset")
    cdna_cleaning_pipeline.save_artifacts(artifact_path)

    # open the pickle file
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(f"{artifact_csv_dir}/{artifact_name}.csv", index=False)


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
    dataset_filepath = Path("raw_dataset/Protein_cDNA_Proteolysis_Dataset/Tsuboyama2023_Dataset2_Dataset3_20230416.csv")
    artifact_path = Path("logs/Protein_cDNA_Proteolysis_dG_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/Protein_cDNA_Proteolysis_dG_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    cdna_cleaning_pipeline = create_cdna_proteolysis_cleaner(dataset_filepath, cdna_cleaning_config)
    cdna_cleaning_pipeline, cdna_dataset = clean_cdna_proteolysis_dataset(
        cdna_cleaning_pipeline
    )

    # Save data
    cdna_dataset.save("cleaned_dataset/cleaned_Protein_cDNA_Proteolysis_dG_Dataset")
    cdna_cleaning_pipeline.save_artifacts(artifact_path)

    # open the pickle file
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(f"{artifact_csv_dir}/{artifact_name}.csv", index=False)


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
            artifact_df.to_csv(f"{artifact_csv_dir}/{artifact_name}.csv", index=False)


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
            artifact_df.to_csv(f"{artifact_csv_dir}/{artifact_name}.csv", index=False)


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()
```

### Advanced Settings

See {py:class}`mutcleaner.cleaners.DdgDtmCleanerConfig` for details.

## ArchStabMS1E10 Epistasis Dataset

### ArchStabMS1E10 Epistasis Sup4 Dataset

The following example shows the complete workflow for downloading, cleaning,
saving the cleaned ArchStabMS1E10 Epistasis Sup4 dataset, and exporting the
cleaning artifacts:
```python
import pickle
from pathlib import Path
from mutcleaner import download_archstabms1e10_source_file
from mutcleaner.cleaners import (
    create_archstabms_1e10_sup4_cleaner,
    clean_archstabms_1e10_sup4_dataset,
)


def main():
    # Prepare dataset
    download_archstabms1e10_source_file("raw_dataset/ArchStabMS1E10_Epistasis_Dataset", overwrite=True, sub_dataset="ArchStabMS1E10_Epistasis_Sup4_Dataset")

    # File settings
    dataset_filepath = Path("raw_dataset/ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Sup4_Dataset.csv")
    artifact_path = Path("logs/ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Sup4_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Sup4_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    archstabms1e10_sup4_cleaning_pipeline = create_archstabms_1e10_sup4_cleaner(dataset_filepath)
    archstabms1e10_sup4_cleaning_pipeline, archstabms1e10_sup4_dataset = clean_archstabms_1e10_sup4_dataset(
        archstabms1e10_sup4_cleaning_pipeline
    )

    # Save data
    archstabms1e10_sup4_dataset.save("cleaned_dataset/cleaned_ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Sup4_Dataset")
    archstabms1e10_sup4_cleaning_pipeline.save_artifacts(artifact_path)

    # open the pickle file
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(f"{artifact_csv_dir}/{artifact_name}.csv", index=False)


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()
```

### Advanced Settings

See {py:class}`mutcleaner.cleaners.ArchStabMS1E10CleanerSup4Config` for details.

### ArchStabMS1E10 Epistasis Sup5 Dataset
The following example shows the complete workflow for downloading, cleaning,
saving the cleaned ArchStabMS1E10 Epistasis Sup5 dataset, and exporting the
cleaning artifacts:
```python
import pickle
from pathlib import Path
from mutcleaner import download_archstabms1e10_source_file
from mutcleaner.cleaners import (
    create_archstabms_1e10_sup5_cleaner,
    clean_archstabms_1e10_sup5_dataset,
)


def main():
    # Prepare dataset
    download_archstabms1e10_source_file("raw_dataset/ArchStabMS1E10_Epistasis_Dataset", overwrite=True, sub_dataset="ArchStabMS1E10_Epistasis_Sup5_Dataset")

    # File settings
    dataset_filepath = Path("raw_dataset/ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Sup5_Dataset.csv")
    artifact_path = Path("logs/ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Sup5_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Sup5_Dataset/")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    archstabms1e10_sup5_cleaning_pipeline = create_archstabms_1e10_sup5_cleaner(dataset_filepath)
    archstabms1e10_sup5_cleaning_pipeline, archstabms1e10_sup5_dataset = clean_archstabms_1e10_sup5_dataset(
        archstabms1e10_sup5_cleaning_pipeline
    )

    # Save data
    archstabms1e10_sup5_dataset.save("cleaned_dataset/cleaned_ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Sup5_Dataset/")
    archstabms1e10_sup5_cleaning_pipeline.save_artifacts(artifact_path)

    # open the pickle file
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(f"{artifact_csv_dir}/{artifact_name}.csv", index=False)


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()
```

### Advanced Settings

See {py:class}`mutcleaner.cleaners.ArchStabMS1E10CleanerSup5Config` for details.

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
    download_antitoxin_pard3_source_file("raw_dataset/Antitoxin_ParD3_Epistasis_Dataset", overwite=True)
    
    # File settings
    dataset_file_path = Path("raw_dataset/Antitoxin_ParD3_Epistasis_Dataset/Antitoxin_ParD3_Epistasis_Dataset.csv")
    artifact_path = Path("logs/Antitoxin_ParD3_Epistasis_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/Antitoxin_ParD3_Epistasis_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    antitoxin_pard3_cleaning_pipeline = create_antitoxin_pard3_cleaner(dataset_file_path)
    antitoxin_pard3_cleaning_pipeline, antitoxin_pard3_dataset = (
        clean_antitoxin_pard3_dataset(antitoxin_pard3_cleaning_pipeline)
    )

    # Save data
    antitoxin_pard3_dataset.save("cleaned_dataset/cleaned_Antitoxin_ParD3_Epistasis_Dataset")
    antitoxin_pard3_cleaning_pipeline.save_artifacts(artifact_path)

    # open the pickle file
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(f"{artifact_csv_dir}/{artifact_name}.csv", index=False)


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
    download_trpb_source_file("raw_dataset/TrpB_Epistasis_Dataset", overwrite=True)
 
    # File settings
    dataset_filepath = Path("raw_dataset/TrpB_Epistasis_Dataset/TrpB_Epistasis_Dataset.csv")
    artifact_path = Path("logs/TrpB_Epistasis_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/TrpB_Epistasis_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    trpb_cleaning_pipeline = create_trpb_cleaner(dataset_filepath)
    trpb_cleaning_pipeline, trpb_dataset = clean_trpb_dataset(
        trpb_cleaning_pipeline
    )

    # Save data
    trpb_dataset.save("cleaned_dataset/cleaned_TrpB_Epistasis_Dataset")
    trpb_cleaning_pipeline.save_artifacts(artifact_path)

    # open the pickle file and read the object
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(f"{artifact_csv_dir}/{artifact_name}.csv", index=False)


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
    download_human_myoglobin_source_file("raw_dataset/Human_Myoglobin_Epistasis_Dataset", overwrite=True)
    # File settings
    dataset_filepath = Path("raw_dataset/Human_Myoglobin_Epistasis_Dataset/Human_Myoglobin_Epistasis_Dataset.csv")
    artifact_path = Path("logs/Human_Myoglobin_Epistasis_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/Human_Myoglobin_Epistasis_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    human_myoglobin_cleaning_pipeline = create_human_myoglobin_cleaner(dataset_filepath)
    human_myoglobin_cleaning_pipeline, human_myoglobin_dataset_dataset = clean_human_myoglobin_dataset(
        human_myoglobin_cleaning_pipeline
    )

    # Save data
    human_myoglobin_dataset_dataset.save("cleaned_dataset/cleaned_Human_Myoglobin_Epistasis_Dataset")
    human_myoglobin_cleaning_pipeline.save_artifacts(artifact_path)

    # open the pickle file and read the object
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(f"{artifact_csv_dir}/{artifact_name}.csv", index=False)


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
    download_ctxm_source_file("raw_dataset/CTXM_Epistasis_Dataset", overwrite=True, sub_dataset="CTXM_Ampicillin_Epistasis_Dataset")

    # File settings
    dataset_filepath = Path("raw_dataset/CTXM_Epistasis_Dataset/Doubles_A3_processed.txt")
    artifact_path = Path("logs/CTXM_Ampicillin_Epistasis_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/CTXM_Ampicillin_Epistasis_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    ctxm_cleaning_pipeline = create_ctxm_cleaner(dataset_filepath)
    ctxm_cleaning_pipeline, ctxm_dataset = clean_ctxm_dataset(ctxm_cleaning_pipeline)

    # Save data
    ctxm_dataset.save("cleaned_dataset/cleaned_CTXM_Epistasis_Dataset")
    ctxm_cleaning_pipeline.save_artifacts(artifact_path)

    # open the pickle file
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(f"{artifact_csv_dir}/{artifact_name}.csv", index=False)


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
from mutcleaner import download_ctxm_source_file
from mutcleaner.cleaners import (
    CTXMCleanerConfig,
    create_ctxm_cleaner,
    clean_ctxm_dataset,
)


def main():
    # Prepare data
    download_ctxm_source_file("raw_dataset/CTXM_Epistasis_Dataset", overwrite=True, sub_dataset="CTXM_Cefotaxime_Epistasis_Dataset")

    # Set cleaning configs
    ctxm_cleaning_config = CTXMCleanerConfig()
    ctxm_cleaning_config.wt_name = "CTXM_cefotaxime"

    # File settings
    dataset_filepath = Path("raw_dataset/CTXM_Epistasis_Dataset/Doubles_C2_processed.txt")
    artifact_path = Path("logs/CTXM_Cefotaxime_Epistasis_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/CTXM_Cefotaxime_Epistasis_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    ctxm_cleaning_pipeline = create_ctxm_cleaner(dataset_filepath, ctxm_cleaning_config)
    ctxm_cleaning_pipeline, ctxm_dataset = clean_ctxm_dataset(ctxm_cleaning_pipeline)

    # Save data
    ctxm_dataset.save("cleaned_dataset/cleaned_CTXM_Epistasis_Dataset")
    ctxm_cleaning_pipeline.save_artifacts(artifact_path)

    # open the pickle file
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(f"{artifact_csv_dir}/{artifact_name}.csv", index=False)


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()

```

### Advanced Settings

See {py:class}`mutcleaner.cleaners.CTXMCleanerConfig` for details.


## RBD-ACE2 Dataset

### Basic Usage

The following example shows the complete workflow for downloading,
cleaning, saving the cleaned RBD-ACE2 dataset, and exporting the
cleaning artifacts:
```python
import pickle
from pathlib import Path

from mutcleaner import download_rbd_ace2_source_file
from mutcleaner.cleaners import (
    create_rbd_ace2_cleaner,
    clean_rbd_ace2_dataset,
)

def main():
    # Prepare data
    download_rbd_ace2_source_file("raw_dataset/RBD_ACE2_Dataset", overwrite=True)

    # File settings
    raw_data_dir = Path("raw_dataset/RBD_ACE2_Dataset")
    dataset_file_paths = sorted(raw_data_dir.glob("*.csv"))

    for dataset_file_path in dataset_file_paths:
        dataset_name = dataset_file_path.stem
        artifact_csv_dir = Path("logs/RBD_ACE2_Dataset") / dataset_name
        artifact_path = artifact_csv_dir / "artifacts.pkl"
        cleaned_dataset_dir = (
            Path("cleaned_dataset/cleaned_RBD_ACE2_Dataset") / dataset_name
        )

        artifact_csv_dir.mkdir(parents=True, exist_ok=True)

        # Clean data
        rbd_ace2_cleaning_pipeline = create_rbd_ace2_cleaner(dataset_file_path)
        rbd_ace2_cleaning_pipeline, rbd_ace2_dataset = clean_rbd_ace2_dataset(
            rbd_ace2_cleaning_pipeline
        )

        # Save data
        rbd_ace2_dataset.save(str(cleaned_dataset_dir))
        rbd_ace2_cleaning_pipeline.save_artifacts(artifact_path)

        # Read artifacts from the pickle file
        with open(artifact_path, "rb") as file:
            artifacts = pickle.load(file)

        for artifact_name, artifact_df in artifacts.items():
            artifact_df.to_csv(f"{artifact_csv_dir}/{artifact_name}.csv", index=False)

if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()

```

You can also download a specific sub-dataset directly (see
{py:func}`mutcleaner.utils.download_rbd_ace2_source_file` for details):

```python
from mutcleaner import download_rbd_ace2_source_file

download_rbd_ace2_source_file(
    "path/to/target/folder",
    sub_dataset="Omicron_EG5_FLip_BA286",
)
```

Supported sub-datasets:
- `Omicron_EG5_FLip_BA286`
- `Omicron_XBB_BQ`
- `Omicron`
- `DMS_variants`
- `Delta`

Alternatively, you can download it from
[Hugging Face](https://huggingface.co/datasets/xulab-research/MutCleaner/tree/main/RBD_ACE2_Dataset).


### Advanced Settings

See {py:class}`mutcleaner.cleaners.RBDACE2CleanerConfig` for details.

## RBD-Antibody Dataset

### Basic Usage

The following example shows the complete workflow for downloading,
cleaning, saving the cleaned RBD-antibody dataset, and exporting the
cleaning artifacts:
```python
import pickle
from pathlib import Path

from mutcleaner import download_rbd_antibody_source_file
from mutcleaner.cleaners import (
    create_rbd_antibody_cleaner,
    clean_rbd_antibody_dataset,
)

def main():
    # Prepare data
    download_rbd_antibody_source_file("raw_dataset/RBD_Antibody_Dataset", overwrite=True)

    # File settings
    raw_data_dir = Path("raw_dataset/RBD_Antibody_Dataset")
    dataset_file_paths = sorted(raw_data_dir.glob("*.csv"))

    for dataset_file_path in dataset_file_paths:
        dataset_name = dataset_file_path.stem
        artifact_csv_dir = Path("logs/RBD_Antibody_Dataset") / dataset_name
        artifact_path = artifact_csv_dir / "artifacts.pkl"
        cleaned_dataset_dir = (
            Path("cleaned_dataset/cleaned_RBD_Antibody_Dataset") / dataset_name
        )

        artifact_csv_dir.mkdir(parents=True, exist_ok=True)

        # Clean data
        rbd_antibody_cleaning_pipeline = create_rbd_antibody_cleaner(dataset_file_path)
        rbd_antibody_cleaning_pipeline, rbd_antibody_dataset = clean_rbd_antibody_dataset(
            rbd_antibody_cleaning_pipeline
        )

        # Save data
        rbd_antibody_dataset.save(str(cleaned_dataset_dir))
        rbd_antibody_cleaning_pipeline.save_artifacts(artifact_path)

        # Read artifacts from the pickle file
        with open(artifact_path, "rb") as file:
            artifacts = pickle.load(file)

        for artifact_name, artifact_df in artifacts.items():
            artifact_df.to_csv(f"{artifact_csv_dir}/{artifact_name}.csv", index=False)

if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()
```

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
    sub_dataset="Moderna",
)
```

Supported sub-datasets:
- `Moderna`
- `Rockefeller`
- `Vir_mAbs`

Alternatively, you can download it from
[Hugging Face](https://huggingface.co/datasets/xulab-research/MutCleaner/tree/main/RBD_Antibody_Dataset).

### Advanced Settings

See {py:class}`mutcleaner.cleaners.RBDAntibodyCleanerConfig` for details.


## Chitosanase dTm Dataset

### Basic Usage

You can download the source file directly by running (see {py:func}`mutcleaner.utils.download_chitosanase_dtm_source_file` for details):
```python
import pickle
from pathlib import Path
from mutcleaner import download_chitosanase_dtm_source_file
from mutcleaner.cleaners import (
    create_chitosanase_dtm_cleaner,
    clean_chitosanase_dtm_dataset,
)


def main():
    # Prepare data
    raw_data_dir = Path("raw_dataset/Chitosanase_dTm_Dataset")
    download_chitosanase_dtm_source_file(raw_data_dir, overwrite=True)

    # File settings
    artifact_path = Path(f"logs/Chitosanase_dTm_Dataset/artifacts.pkl")
    artifact_csv_dir = Path(f"logs/Chitosanase_dTm_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    chitosanase_cleaning_pipeline = create_chitosanase_dtm_cleaner("raw_dataset/Chitosanase_dTm_Dataset/Chitosanase_dTm_Dataset.csv")
    chitosanase_cleaning_pipeline, chitosanase_dataset = clean_chitosanase_dtm_dataset(chitosanase_cleaning_pipeline)

    # Save data
    chitosanase_dataset.save(f"cleaned_dataset/cleaned_Chitosanase_dTm_Dataset")
    chitosanase_cleaning_pipeline.save_artifacts(artifact_path)

    # open the pickle file
    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(f"{artifact_csv_dir}/{artifact_name}.csv", index=False)


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()

```

### Advanced Settings

See {py:class}`mutcleaner.cleaners.ChitosanasedTmCleanerConfig` for details.


## MGnify ddG Dataset

### Basic Usage

You can download the source file directly by running (see {py:func}`mutcleaner.utils.download_mgnify_ddg_source_file` for details):
```python
import pickle
from pathlib import Path
from mutcleaner import download_mgnify_ddg_source_file
from mutcleaner.cleaners import (
    create_mgnify_ddg_cleaner,
    clean_mgnify_ddg_dataset,
)


def main():
    # Prepare data
    download_mgnify_ddg_source_file("raw_dataset/MGnify_ddG_Dataset", overwrite=True)

    artifact_path = Path("logs/MGnify_ddG_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/MGnify_ddG_Dataset")
    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    mgnify_cleaning_pipeline = create_mgnify_ddg_cleaner(Path("raw_dataset/MGnify_ddG_Dataset/MGnify_ddG_Dataset.csv"))
    mgnify_cleaning_pipeline, mgnify_dataset = clean_mgnify_ddg_dataset(mgnify_cleaning_pipeline)
    
    # Save data
    mgnify_dataset.save("cleaned_dataset/cleaned_MGnify_ddG_Dataset")
    mgnify_cleaning_pipeline.save_artifacts(artifact_path)

    with open(artifact_path, "rb") as file:
        artifacts = pickle.load(file)

    for artifact_name, artifact_df in artifacts.items():
        artifact_df.to_csv(f"{artifact_csv_dir}/{artifact_name}.csv", index=False)


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.freeze_support()
    main()

```

### Advanced Settings

See {py:class}`mutcleaner.cleaners.MGnifyddGCleanerConfig` for details.
