# Supported Datasets

## Human Domainome Dataset

### Human Domainome Sup2 Dataset
The following example shows the complete workflow for downloading, cleaning,
saving the cleaned Human Domainome sup2 dataset, and exporting the
cleaning artifacts: 
```python
import pickle
from pathlib import Path

import mutcleaner


def main():
    # Prepare data
    mutcleaner.download_human_domainome_source_file("raw_dataset/Human_Domainome_Dataset", overwrite=True, sub_dataset="Human Domainome Sup2 Dataset")
    # File settings
    dataset_filepath = Path("raw_dataset/Human_Domainome_Dataset/SupplementaryTable2.txt")
    artifact_path = Path("logs/Human_Domainome_Dataset/Human_Domainome_Sup2_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/Human_Domainome_Dataset/Human_Domainome_Sup2_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    hd_cleaning_pipeline = mutcleaner.cleaners.create_human_domainome_sup2_cleaner(dataset_filepath, config=None)
    hd_cleaning_pipeline, hd_dataset = mutcleaner.cleaners.clean_human_domainome_sup2_dataset(hd_cleaning_pipeline)

    # Save data
    hd_dataset.metadata.update({
        "source_doi": "10.1038/s41586-024-08370-4",
        "source_url": "https://www.nature.com/articles/s41586-024-08370-4",
        "measured_property": "Expression",
        "measurement_type": "Normalized aPCA fitness",
        "label_unit": "dimensionless",
        "label_direction": "higher_is_more",
        "assay": "Abundance protein fragment complementation assay (aPCA)",
        "mutcleaner_version": mutcleaner.__version__
    })
    hd_dataset.save("cleaned_dataset/cleaned_Human_Domainome_Dataset/Human_Domainome_Sup2_Dataset", save_type="mutcleaner")
    hd_cleaning_pipeline.save_artifacts(artifact_path, format="pickle")

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

import mutcleaner


def main():
    # File settings
    mutcleaner.download_human_domainome_source_file("raw_dataset/Human_Domainome_Dataset", overwrite=True, sub_dataset="Human Domainome Sup4 Dataset")
    artifact_path = Path("logs/Human_Domainome_Dataset/Human_Domainome_Sup4_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/Human_Domainome_Dataset/Human_Domainome_Sup4_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    hd_cleaning_pipeline = mutcleaner.cleaners.create_human_domainome_sup4_cleaner(dataset_or_path="raw_dataset/Human_Domainome_Dataset/SupplementaryTable4.txt", sequence_source="raw_dataset/Human_Domainome_Dataset/wild_type.fasta", config=None)
    hd_cleaning_pipeline, hd_dataset = mutcleaner.cleaners.clean_human_domainome_sup4_dataset(hd_cleaning_pipeline)

    # Save data
    hd_dataset.metadata.update({
        "source_doi": "10.1038/s41586-024-08370-4",
        "source_url": "https://www.nature.com/articles/s41586-024-08370-4",
        "measured_property": "Stability",
        "measurement_type": "Sign-reversed scaled homologue-averaged folding ΔΔG",
        "label_unit": "dimensionless",
        "label_direction": "higher_is_more",
        "assay": "Abundance protein fragment complementation assay (aPCA) with MoCHI thermodynamic modelling",
        "mutcleaner_version": mutcleaner.__version__
    })
    hd_dataset.save("cleaned_dataset/cleaned_Human_Domainome_Dataset/Human_Domainome_Sup4_Dataset", save_type="mutcleaner")
    hd_cleaning_pipeline.save_artifacts(artifact_path, format="pickle")

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

## Protein cDNA Proteolysis Dataset


### ΔΔG as Label (Default Pipeline)

The following example shows the complete workflow for downloading, cleaning,
saving the cleaned cDNA Proteolysis dataset, and exporting the
cleaning artifacts, using ΔΔG as the label:
```python
import pickle
from pathlib import Path

import mutcleaner


def main():
    # Prepare dataset
    mutcleaner.download_protein_cdna_proteolysis_source_file("raw_dataset/Protein_cDNA_Proteolysis_Dataset", overwrite=True)

    # File settings
    dataset_filepath = Path("raw_dataset/Protein_cDNA_Proteolysis_Dataset/Tsuboyama2023_Dataset2_Dataset3_20230416.csv")
    artifact_path = Path("logs/Protein_cDNA_Proteolysis_ddG_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/Protein_cDNA_Proteolysis_ddG_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    cdna_cleaning_pipeline = mutcleaner.cleaners.create_cdna_proteolysis_cleaner(dataset_filepath)
    cdna_cleaning_pipeline, cdna_dataset = mutcleaner.cleaners.clean_cdna_proteolysis_dataset(cdna_cleaning_pipeline)

    # Save data
    cdna_dataset.metadata.update({
        "source_doi": "10.1038/s41586-023-06328-6",
        "source_url": "https://www.nature.com/articles/s41586-023-06328-6",
        "measured_property": "Stability",
        "measurement_type": "ΔΔG",
        "label_unit": "kcal/mol",
        "label_direction": "higher_is_more",
        "assay": "cDNA display proteolysis",
        "mutcleaner_version": mutcleaner.__version__
    })
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

import mutcleaner


def main():
    # Set cleaning configs
    cdna_cleaning_config = mutcleaner.cleaners.CDNAProteolysisCleanerConfig()
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
    cdna_cleaning_pipeline = mutcleaner.cleaners.create_cdna_proteolysis_cleaner(dataset_filepath, cdna_cleaning_config)
    cdna_cleaning_pipeline, cdna_dataset = mutcleaner.cleaners.clean_cdna_proteolysis_dataset(cdna_cleaning_pipeline)

    # Save data
    cdna_dataset.metadata.update({
        "source_doi": "10.1038/s41586-023-06328-6",
        "source_url": "https://www.nature.com/articles/s41586-023-06328-6",
        "measured_property": "Stability",
        "measurement_type": "ΔG",
        "label_unit": "kcal/mol",
        "label_direction": "higher_is_more",
        "assay": "cDNA display proteolysis",
        "mutcleaner_version": mutcleaner.__version__
    })
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

import mutcleaner


def main():
    # Prepare data
    downloaded_files = mutcleaner.download_ddg_dtm_source_file("raw_dataset/ddG_Dataset", dataset_type="ddg", overwrite=True, sub_dataset=None)

    # File settings
    for dataset_filepath in sorted(downloaded_files.values()):
        data_file = Path(dataset_filepath).stem
        artifact_path = Path(f"logs/ddG_Dataset/{data_file}/artifacts.pkl")
        artifact_csv_dir = Path(f"logs/ddG_Dataset/{data_file}")

        artifact_csv_dir.mkdir(parents=True, exist_ok=True)

        # Clean data
        ddgdtm_cleaning_pipeline = mutcleaner.cleaners.create_ddg_dtm_cleaner(dataset_filepath, config=None)
        ddgdtm_cleaning_pipeline, ddgdtm_dataset = mutcleaner.cleaners.clean_ddg_dtm_dataset(ddgdtm_cleaning_pipeline)

        # Save data
        ddgdtm_dataset.metadata.update({
            "source_doi": "10.1038/s43588-024-00716-2",
            "source_url": "https://www.nature.com/articles/s43588-024-00716-2",
            "measured_property": "Stability",
            "measurement_type": "ΔΔG",
            "label_unit": "kcal/mol",
            "label_direction": "higher_is_more",
            "assay": "Protein unfolding measurements (multiple methods)",
            "mutcleaner_version": mutcleaner.__version__
        })
        ddgdtm_dataset.save(f"cleaned_dataset/cleaned_ddG_Dataset/{data_file}", save_type="mutcleaner")
        ddgdtm_cleaning_pipeline.save_artifacts(artifact_path, format="pickle")

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

import mutcleaner


def main():
    # Prepare data
    downloaded_files = mutcleaner.download_ddg_dtm_source_file("raw_dataset/dTm_Dataset", dataset_type="dtm", overwrite=True, sub_dataset=None)

    # File settings
    for dataset_filepath in sorted(downloaded_files.values()):
        data_file = Path(dataset_filepath).stem
        artifact_path = Path(f"logs/dTm_Dataset/{data_file}/artifacts.pkl")
        artifact_csv_dir = Path(f"logs/dTm_Dataset/{data_file}")

        artifact_csv_dir.mkdir(parents=True, exist_ok=True)

        # Clean data
        ddgdtm_cleaning_pipeline = mutcleaner.cleaners.create_ddg_dtm_cleaner(dataset_filepath, config=None)
        ddgdtm_cleaning_pipeline, ddgdtm_dataset = mutcleaner.cleaners.clean_ddg_dtm_dataset(ddgdtm_cleaning_pipeline)

        # Save data
        ddgdtm_dataset.metadata.update({
            "source_doi": "10.1038/s43588-024-00716-2",
            "source_url": "https://www.nature.com/articles/s43588-024-00716-2",
            "measured_property": "Stability",
            "measurement_type": "ΔTm",
            "label_unit": "°C",
            "label_direction": "higher_is_more",
            "assay": "Thermal denaturation measurements (multiple methods)",
            "mutcleaner_version": mutcleaner.__version__
        })
        ddgdtm_dataset.save(f"cleaned_dataset/cleaned_dTm_Dataset/{data_file}", save_type="mutcleaner")
        ddgdtm_cleaning_pipeline.save_artifacts(artifact_path, format="pickle")

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

Supplementary Table 4 reports DiMSum fitness scores for intracellular protein
abundance (AbundancePCA) and ligand binding (BindingPCA). Higher scores
indicate greater abundance or binding in the corresponding assay.
[Faure et al., Supplementary Table 4](https://www.nature.com/articles/s41586-024-07966-0#Sec33)

The following example shows the complete workflow for downloading, cleaning,
saving the cleaned ArchStabMS1E10 Epistasis Sup4 dataset, and exporting the
cleaning artifacts:
```python
import pickle
from pathlib import Path

import mutcleaner


def main():
    # Prepare dataset
    mutcleaner.download_archstabms1e10_source_file("raw_dataset/ArchStabMS1E10_Epistasis_Dataset", overwrite=True, sub_dataset="ArchStabMS1E10_Epistasis_Sup4_Dataset")

    # File settings
    dataset_filepath = Path("raw_dataset/ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Sup4_Dataset.csv")
    artifact_path = Path("logs/ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Sup4_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Sup4_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    archstabms1e10_sup4_cleaning_pipeline = mutcleaner.cleaners.create_archstabms_1e10_sup4_cleaner(dataset_filepath, config=None)
    archstabms1e10_sup4_cleaning_pipeline, archstabms1e10_sup4_dataset = mutcleaner.cleaners.clean_archstabms_1e10_sup4_dataset(archstabms1e10_sup4_cleaning_pipeline)

    # Save data
    assay_metadata = {
        "Abundance": ("Expression", "Abundance fitness", "AbundancePCA"),
        "Binding": ("Binding", "Binding fitness", "BindingPCA")
    }
    for reference_id in archstabms1e10_sup4_dataset.list_reference_sequences():
        phenotype = reference_id.rsplit("_", 1)[1]
        measured_property, measurement_type, assay = assay_metadata[phenotype]
        reference_dataset = archstabms1e10_sup4_dataset.filter_by_reference(reference_id)
        reference_dataset.name = archstabms1e10_sup4_dataset.name
        reference_dataset.metadata.update({
            "source_doi": "10.1038/s41586-024-07966-0",
            "source_url": "https://www.nature.com/articles/s41586-024-07966-0",
            "measured_property": measured_property,
            "measurement_type": measurement_type,
            "label_unit": "dimensionless",
            "label_direction": "higher_is_more",
            "assay": assay,
            "mutcleaner_version": mutcleaner.__version__
        })
        reference_dataset.save("cleaned_dataset/cleaned_ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Sup4_Dataset", save_type="mutcleaner")
    archstabms1e10_sup4_cleaning_pipeline.save_artifacts(artifact_path, format="pickle")

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

Supplementary Table 5 contains MoCHI estimates of folding and binding free
energy changes. For double mutants, MutCleaner adds the coupling term to the
two single-mutant effects to obtain the total ΔΔG relative to wild type.
Lower ΔΔG indicates greater stability or stronger binding. Library 3 uses
AbundancePCA and BindingPCA jointly to infer both traits.
[Faure et al., Fig. 4](https://www.nature.com/articles/s41586-024-07966-0#Fig4)

The following example shows the complete workflow for downloading, cleaning,
saving the cleaned ArchStabMS1E10 Epistasis Sup5 dataset, and exporting the
cleaning artifacts:
```python
import pickle
from pathlib import Path

import mutcleaner


def main():
    # Prepare dataset
    mutcleaner.download_archstabms1e10_source_file("raw_dataset/ArchStabMS1E10_Epistasis_Dataset", overwrite=True, sub_dataset="ArchStabMS1E10_Epistasis_Sup5_Dataset")

    # File settings
    dataset_filepath = Path("raw_dataset/ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Sup5_Dataset.csv")
    artifact_path = Path("logs/ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Sup5_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Sup5_Dataset/")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    archstabms1e10_sup5_cleaning_pipeline = mutcleaner.cleaners.create_archstabms_1e10_sup5_cleaner(dataset_filepath, config=None)
    archstabms1e10_sup5_cleaning_pipeline, archstabms1e10_sup5_dataset = mutcleaner.cleaners.clean_archstabms_1e10_sup5_dataset(archstabms1e10_sup5_cleaning_pipeline)

    # Save data
    trait_metadata = {
        "Folding": ("Stability", "Folding free energy change (ΔΔG_f)"),
        "Binding": ("Binding", "Binding free energy change (ΔΔG_b)")
    }
    for reference_id in archstabms1e10_sup5_dataset.list_reference_sequences():
        library, trait = reference_id.split("_")
        measured_property, measurement_type = trait_metadata[trait]
        reference_dataset = archstabms1e10_sup5_dataset.filter_by_reference(reference_id)
        reference_dataset.name = archstabms1e10_sup5_dataset.name
        reference_dataset.metadata.update({
            "source_doi": "10.1038/s41586-024-07966-0",
            "source_url": "https://www.nature.com/articles/s41586-024-07966-0",
            "measured_property": measured_property,
            "measurement_type": measurement_type,
            "label_unit": "kcal/mol",
            "label_direction": "lower_is_more",
            "assay": "AbundancePCA and BindingPCA" if library == "3" else "AbundancePCA",
            "mutcleaner_version": mutcleaner.__version__
        })
        reference_dataset.save("cleaned_dataset/cleaned_ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Sup5_Dataset", save_type="mutcleaner")
    archstabms1e10_sup5_cleaning_pipeline.save_artifacts(artifact_path, format="pickle")

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

Lite et al. constructed the three-site ParD3 library. Ding et al. (2022)
measured how well its variants neutralized wild-type ParE3 in a pooled
competitive growth assay. [Library paper](https://doi.org/10.7554/eLife.60924) · [Measurement paper, Fig. 6](https://www.nature.com/articles/s41559-022-01688-0#Fig6)

Because the assay measures cell growth, we use `Organismal Fitness` for
`measured_property`, following the ProteinGym benchmark website.
[ProteinGym benchmarks](https://proteingym.org/benchmarks)

The input `label` comes from the `E3` column of the three-site table
distributed with CoVES. MutCleaner subtracts the wild-type ParD3 score from
each variant's normalized fitness. The difference is dimensionless.
A value of zero matches wild-type fitness,
and higher values indicate stronger toxin neutralization in this assay.
[CoVES data table](https://github.com/ddingding/CoVES/blob/main/data/DMS_data/df_mut_all_norm.csv)

### Basic Usage

The following example shows the complete workflow for downloading, cleaning,
saving the cleaned Antitoxin ParD3 Epistasis dataset, and exporting the
cleaning artifacts:
```python
import pickle
from pathlib import Path

import mutcleaner


def main():
    # Prepare dataset
    mutcleaner.download_antitoxin_pard3_source_file("raw_dataset/Antitoxin_ParD3_Epistasis_Dataset", overwrite=True)

    # File settings
    dataset_file_path = Path("raw_dataset/Antitoxin_ParD3_Epistasis_Dataset/Antitoxin_ParD3_Epistasis_Dataset.csv")
    artifact_path = Path("logs/Antitoxin_ParD3_Epistasis_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/Antitoxin_ParD3_Epistasis_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    antitoxin_pard3_cleaning_pipeline = mutcleaner.cleaners.create_antitoxin_pard3_cleaner(dataset_file_path, config=None)
    antitoxin_pard3_cleaning_pipeline, antitoxin_pard3_dataset = mutcleaner.cleaners.clean_antitoxin_pard3_dataset(antitoxin_pard3_cleaning_pipeline)

    # Save data
    antitoxin_pard3_dataset.metadata.update({
        "source_doi": "10.1038/s41559-022-01688-0",
        "source_url": "https://www.nature.com/articles/s41559-022-01688-0",
        "measured_property": "Organismal Fitness",
        "measurement_type": "Normalized fitness difference (variant minus wild type)",
        "label_unit": "dimensionless",
        "label_direction": "higher_is_more",
        "assay": "Pooled competitive growth assay",
        "mutcleaner_version": mutcleaner.__version__
    })
    antitoxin_pard3_dataset.save("cleaned_dataset/cleaned_Antitoxin_ParD3_Epistasis_Dataset", save_type="mutcleaner")
    antitoxin_pard3_cleaning_pipeline.save_artifacts(artifact_path, format="pickle")

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

Of the 160,000 amino acid combinations without a stop codon, the authors
experimentally quantified 159,129 and imputed fitness scores for the
remaining 871 using k-nearest neighbors.
The cleaned label is the variant's normalized growth fitness minus that of
the experimental parent (Tm9D8*, abbreviated VFVS). Positive values indicate
higher growth fitness than the parent.

The following example shows the complete workflow for downloading, cleaning,
saving the cleaned TrpB Epistasis dataset, and exporting the
cleaning artifacts:
```python
import pickle
from pathlib import Path

import mutcleaner


def main():
    # Prepare data
    dataset_filepath = Path("raw_dataset/TrpB_Epistasis_Dataset/TrpB_Epistasis_Dataset.csv")
    if not dataset_filepath.is_file():
        mutcleaner.download_trpb_source_file(str(dataset_filepath.parent), overwrite=False)

    # File settings
    artifact_path = Path("logs/TrpB_Epistasis_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/TrpB_Epistasis_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    trpb_cleaning_pipeline = mutcleaner.cleaners.create_trpb_cleaner(dataset_filepath, config=None)
    trpb_cleaning_pipeline, trpb_dataset = mutcleaner.cleaners.clean_trpb_dataset(trpb_cleaning_pipeline)

    # Save data
    trpb_dataset.metadata.update({
        "source_doi": "10.1073/pnas.2400439121",
        "source_url": "https://www.pnas.org/doi/10.1073/pnas.2400439121",
        "measured_property": "Organismal Fitness",
        "measurement_type": "Normalized growth fitness difference (variant minus parent)",
        "label_unit": "dimensionless",
        "label_direction": "higher_is_more",
        "assay": "Pooled-culture enrichment assay",
        "mutcleaner_version": mutcleaner.__version__
    })
    trpb_dataset.save("cleaned_dataset/cleaned_TrpB_Epistasis_Dataset", save_type="mutcleaner")
    trpb_cleaning_pipeline.save_artifacts(artifact_path, format="pickle")

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

## Protein Human Myoglobin Epistasis Dataset

### Basic Usage

The following example shows the complete workflow for downloading, cleaning,
saving the cleaned Protein Human Myoglobin Epistasis dataset, and exporting the
cleaning artifacts:
```python
import pickle
from pathlib import Path

import mutcleaner


def main():
    # Prepare data
    mutcleaner.download_protein_human_myoglobin_source_file("raw_dataset/Protein_Human_Myoglobin_Epistasis_Dataset", overwrite=True)
    # File settings
    dataset_filepath = Path("raw_dataset/Protein_Human_Myoglobin_Epistasis_Dataset/Protein_Human_Myoglobin_Epistasis_Dataset.csv")
    artifact_path = Path("logs/Protein_Human_Myoglobin_Epistasis_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/Protein_Human_Myoglobin_Epistasis_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    pipeline = mutcleaner.cleaners.create_protein_human_myoglobin_cleaner(dataset_filepath, config=None)
    pipeline, dataset = mutcleaner.cleaners.clean_protein_human_myoglobin_dataset(pipeline)

    # Save data
    dataset.metadata.update({
        "source_doi": "10.1002/pro.70112",
        "source_url": "https://onlinelibrary.wiley.com/doi/full/10.1002/pro.70112",
        "measured_property": "Expression",
        "measurement_type": "Expression fitness score",
        "label_unit": "dimensionless",
        "label_direction": "higher_is_more",
        "assay": "Yeast surface display, fluorescence-activated cell sorting (FACS) and deep sequencing",
        "mutcleaner_version": mutcleaner.__version__
    })
    dataset.save("cleaned_dataset/cleaned_Protein_Human_Myoglobin_Epistasis_Dataset", save_type="mutcleaner")
    pipeline.save_artifacts(artifact_path, format="pickle")

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

See {py:class}`mutcleaner.cleaners.ProteinHumanMyoglobinCleanerConfig` for details.

## CTXM Epistasis Dataset


### CTXM Ampicillin Epistasis Dataset

The following example shows the complete workflow for downloading, cleaning,
saving the cleaned CTXM Ampicillin Epistasis Dataset, and exporting the
cleaning artifacts:
```python
import pickle
from pathlib import Path

import mutcleaner


def main():
    # Prepare data
    mutcleaner.download_ctxm_source_file("raw_dataset/CTXM_Epistasis_Dataset", overwrite=True, sub_dataset="CTXM_Ampicillin_Epistasis_Dataset")

    # File settings
    dataset_filepath = Path("raw_dataset/CTXM_Epistasis_Dataset/Doubles_A3_processed.txt")
    artifact_path = Path("logs/CTXM_Ampicillin_Epistasis_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/CTXM_Ampicillin_Epistasis_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    ctxm_cleaning_pipeline = mutcleaner.cleaners.create_ctxm_cleaner(dataset_filepath)
    ctxm_cleaning_pipeline, ctxm_dataset = mutcleaner.cleaners.clean_ctxm_dataset(ctxm_cleaning_pipeline)

    # Save data
    ctxm_dataset.metadata.update({
        "source_doi": "10.1073/pnas.2313513121",
        "source_url": "https://www.pnas.org/doi/10.1073/pnas.2313513121",
        "measured_property": "Organismal Fitness",
        "measurement_type": "Relative fitness (log10 enrichment relative to wild type)",
        "label_unit": "dimensionless",
        "label_direction": "higher_is_more",
        "assay": "Ampicillin selection and deep sequencing",
        "mutcleaner_version": mutcleaner.__version__
    })
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

import mutcleaner


def main():
    # Prepare data
    mutcleaner.download_ctxm_source_file("raw_dataset/CTXM_Epistasis_Dataset", overwrite=True, sub_dataset="CTXM_Cefotaxime_Epistasis_Dataset")

    # Set cleaning configs
    ctxm_cleaning_config = mutcleaner.cleaners.CTXMCleanerConfig()
    ctxm_cleaning_config.wt_name = "CTXM_cefotaxime"

    # File settings
    dataset_filepath = Path("raw_dataset/CTXM_Epistasis_Dataset/Doubles_C2_processed.txt")
    artifact_path = Path("logs/CTXM_Cefotaxime_Epistasis_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/CTXM_Cefotaxime_Epistasis_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    ctxm_cleaning_pipeline = mutcleaner.cleaners.create_ctxm_cleaner(dataset_filepath, ctxm_cleaning_config)
    ctxm_cleaning_pipeline, ctxm_dataset = mutcleaner.cleaners.clean_ctxm_dataset(ctxm_cleaning_pipeline)

    # Save data
    ctxm_dataset.metadata.update({
        "source_doi": "10.1073/pnas.2313513121",
        "source_url": "https://www.pnas.org/doi/10.1073/pnas.2313513121",
        "measured_property": "Organismal Fitness",
        "measurement_type": "Relative fitness (log10 enrichment relative to wild type)",
        "label_unit": "dimensionless",
        "label_direction": "higher_is_more",
        "assay": "Cefotaxime selection and deep sequencing",
        "mutcleaner_version": mutcleaner.__version__
    })
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

## Codon cDNA Proteolysis Dataset

### Basic Usage

You can download the source file directly by running (see {py:func}`mutcleaner.utils.download_codon_cdna_proteolysis_source_file` for details):
```python
import pickle
from pathlib import Path

import mutcleaner


def main():
    # download the Codon cDNA Proteolysis dataset
    mutcleaner.download_codon_cdna_proteolysis_source_file("raw_dataset/Codon_cDNA_Proteolysis_Dataset", overwrite=True)

    # File settings
    dataset_filepath = Path("raw_dataset/Codon_cDNA_Proteolysis_Dataset/Codon_cDNA_Proteolysis_Dataset.csv")
    artifact_path = Path("logs/Codon_cDNA_Proteolysis_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/Codon_cDNA_Proteolysis_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    codon_cdna_cleaning_pipeline = mutcleaner.cleaners.create_codon_cdna_proteolysis_cleaner(dataset_filepath, sequence_source="raw_dataset/Codon_cDNA_Proteolysis_Dataset/wt.fasta")
    codon_cdna_cleaning_pipeline, codon_cdna_dataset = mutcleaner.cleaners.clean_codon_cdna_proteolysis_dataset(codon_cdna_cleaning_pipeline)

    # Save data
    # The provider of the preprocessed data identifies fitness as the original ddG_ML column.
    codon_cdna_dataset.metadata.update({
        "source_doi": "10.1038/s41586-023-06328-6",
        "source_url": "https://www.nature.com/articles/s41586-023-06328-6",
        "measured_property": "Stability",
        "measurement_type": "ΔΔG",
        "label_unit": "kcal/mol",
        "label_direction": "higher_is_more",
        "assay": "cDNA display proteolysis",
        "mutcleaner_version": mutcleaner.__version__
    })
    codon_cdna_dataset.save("cleaned_dataset/cleaned_Codon_cDNA_Proteolysis_Dataset")
    codon_cdna_cleaning_pipeline.save_artifacts(artifact_path)

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

See {py:class}`mutcleaner.cleaners.CodoncDNAProteolysisCleanerConfig` for details.

## Codon DMS Substitutions Dataset

### Basic Usage

You can download the source file directly by running (see {py:func}`mutcleaner.utils.download_codon_dms_substitutions_source_file` for details):
```python
import pickle
from pathlib import Path
from mutcleaner import download_codon_dms_substitutions_source_file
from mutcleaner.cleaners import (
    create_codon_dms_substitutions_cleaner,
    clean_codon_dms_substitutions_dataset,
)


def main():
    # download the Codon DMS Substitutions Dataset
    download_codon_dms_substitutions_source_file("raw_dataset/Codon_DMS_Substitutions_Dataset", overwrite=True)

    # File settings
    dataset_filepath = Path("raw_dataset/Codon_DMS_Substitutions_Dataset/Codon_DMS_Substitutions_Dataset.zip")
    artifact_path = Path("logs/Codon_DMS_Substitutions_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/Codon_DMS_Substitutions_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    codon_dms_substitutions_cleaning_pipeline = create_codon_dms_substitutions_cleaner(dataset_filepath)
    codon_dms_substitutions_cleaning_pipeline, codon_dms_substitutions_dataset = clean_codon_dms_substitutions_dataset(
        codon_dms_substitutions_cleaning_pipeline
    )

    # Save data
    codon_dms_substitutions_dataset.save("cleaned_dataset/cleaned_Codon_DMS_Substitutions_Dataset")
    codon_dms_substitutions_cleaning_pipeline.save_artifacts(artifact_path)

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

See {py:class}`mutcleaner.cleaners.CodonDMSSubstitutionsCleanerConfig` for details.

## Codon Human Myoglobin Epistasis Dataset

### Basic Usage

The following example shows the complete workflow for downloading, cleaning,
saving the cleaned Codon Human Myoglobin Epistasis dataset, and exporting the
cleaning artifacts:
```python
import pickle
from pathlib import Path

import mutcleaner


def main():
    # Prepare data
    mutcleaner.download_codon_human_myoglobin_source_file("raw_dataset/Codon_Human_Myoglobin_Epistasis_Dataset", overwrite=True)
    # File settings
    dataset_filepath = Path("raw_dataset/Codon_Human_Myoglobin_Epistasis_Dataset/Codon_Human_Myoglobin_Epistasis_Dataset.csv")
    artifact_path = Path("logs/Codon_Human_Myoglobin_Epistasis_Dataset/artifacts.pkl")
    artifact_csv_dir = Path("logs/Codon_Human_Myoglobin_Epistasis_Dataset")

    artifact_csv_dir.mkdir(parents=True, exist_ok=True)

    # Clean data
    pipeline = mutcleaner.cleaners.create_codon_human_myoglobin_cleaner(dataset_filepath, config=None)
    pipeline, dataset = mutcleaner.cleaners.clean_codon_human_myoglobin_dataset(pipeline)

    # Save data
    dataset.metadata.update({
        "source_doi": "10.1002/pro.70112",
        "source_url": "https://onlinelibrary.wiley.com/doi/full/10.1002/pro.70112",
        "measured_property": "Expression",
        "measurement_type": "Expression fitness score",
        "label_unit": "dimensionless",
        "label_direction": "higher_is_more",
        "assay": "Yeast surface display, fluorescence-activated cell sorting (FACS) and deep sequencing",
        "mutcleaner_version": mutcleaner.__version__
    })
    dataset.save("cleaned_dataset/cleaned_Codon_Human_Myoglobin_Epistasis_Dataset", save_type="mutcleaner")
    pipeline.save_artifacts(artifact_path, format="pickle")

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

See {py:class}`mutcleaner.cleaners.CodonHumanMyoglobinCleanerConfig` for details.
