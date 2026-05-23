# MutCleaner

[![License badge](https://img.shields.io/badge/License-BSD_3--Clause-yellow?logo=opensourceinitiative&logoColor=white)](https://github.com/xulab-research/MutCleaner/blob/main/LICENSE)
[![PyPI version badge](https://img.shields.io/pypi/v/mutcleaner?logo=python&logoColor=white&color=orange)](https://pypi.org/project/mutcleaner/)
[![Docs](https://github.com/xulab-research/MutCleaner/actions/workflows/docs.yml/badge.svg)](https://xulab-research.github.io/MutCleaner/)

MutCleaner is an extensible Python toolkit for cleaning and standardizing biological mutation datasets, integrating dataset-specific cleaning pipelines with core abstractions for protein, nucleotide, and codon-level mutation representations.

* **Documentation**: https://xulab-research.github.io/MutCleaner/
* **Cleaning Examples**: https://xulab-research.github.io/MutCleaner/user_guide/cleaners.html

## Overview

MutCleaner is an extensible Python toolkit for cleaning, standardizing, and analyzing biological mutation datasets. It currently focuses on protein variant data while providing core abstractions for DNA, RNA, protein sequences, and codon-level mutation representations.

The package combines dataset-specific cleaning pipelines with reusable sequence and mutation utilities, enabling reproducible preprocessing of large-scale mutational datasets for downstream bioinformatics and machine learning analyses.

### Key Capabilities

- **Mutation dataset cleaning and standardization**: Harmonize mutation annotations, sequences, labels, and metadata across heterogeneous biological mutation datasets.
- **Sequence representation and validation**: Utilities for DNA, RNA, and protein sequences, including validation, transcription, reverse transcription, translation, and mutation application.
- **Mutation parsing and transformation**: Tools for parsing amino-acid and codon-level mutations, inferring mutations from sequences, applying mutations to reference sequences, and converting codon mutations into amino-acid changes.
- **Modular pipeline architecture**: A composable pipeline interface for building reproducible dataset-cleaning workflows.
- **Parallel and scalable dataset processing**: Multi-core utilities for mutation validation, mutation application, and sequence-based mutation inference, supporting efficient processing of large tabular mutation datasets.

## Installation

### Requirements

- Python 3.13+
- Dependencies are automatically installed via pip.

### Install via pip

```bash
pip install mutcleaner
```

### Development Installation

```bash
git clone https://github.com/xulab-research/MutCleaner.git MutCleaner
cd MutCleaner
pip install -e .
```

To install development dependencies for testing and documentation:

```bash
pip install -e ".[dev]"
```

## Quick Start

See the [Data Cleaners Usage Guide](https://xulab-research.github.io/MutCleaner/user_guide/cleaners.html) for more examples.

### Supported Datasets

<table>
  <thead>
    <tr>
      <th>Dataset Name</th>
      <th>Reference</th>
      <th>File</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Human Domainome Dataset</td>
      <td><a href="https://doi.org/10.1038/s41586-024-08370-4">Site-saturation mutagenesis of 500 human protein domains</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/Human_Domainome_Dataset/SupplementaryTable2.txt">SupplementaryTable2.txt</a>, <a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/Human_Domainome_Dataset/wild_type.fasta">wild_type.fasta</a></td>
    </tr>
    <tr>
      <td>ProteinGym DMS Substitutions Dataset</td>
      <td><a href="https://doi.org/10.1101/2023.12.07.570727">ProteinGym: Large-Scale Benchmarks for Protein Design and Fitness Prediction</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/ProteinGym_DMS_Substitutions_Dataset/DMS_ProteinGym_substitutions.zip">DMS_ProteinGym_substitutions.zip</a></td>
    </tr>
    <tr>
      <td>cDNA Proteolysis Dataset</td>
      <td><a href="https://doi.org/10.1038/s41586-023-06328-6">Mega-scale experimental analysis of protein folding stability in biology and design</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/cDNA_Proteolysis_Dataset/Tsuboyama2023_Dataset2_Dataset3_20230416.csv">Tsuboyama2023_Dataset2_Dataset3_20230416.csv</a></td>
    </tr>
    <tr>
      <td>ΔΔG Dataset</td>
      <td><a href="https://doi.org/10.1038/s43588-024-00716-2">Improving the prediction of protein stability changes upon mutations by geometric learning and a pre-training strategy</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/%CE%94%CE%94G_Dataset/M1261.csv">M1261.csv</a>, <a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/%CE%94%CE%94G_Dataset/S461.csv">S461.csv</a>, <a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/%CE%94%CE%94G_Dataset/S669.csv">S669.csv</a>, <a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/%CE%94%CE%94G_Dataset/S783.csv">S783.csv</a>, <a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/%CE%94%CE%94G_Dataset/S8754.csv">S8754.csv</a></td>
    </tr>
    <tr>
      <td>ΔTm Dataset</td>
      <td><a href="https://doi.org/10.1038/s43588-024-00716-2">Improving the prediction of protein stability changes upon mutations by geometric learning and a pre-training strategy</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/%CE%94Tm_Dataset/S4346.csv">S4346.csv</a>, <a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/%CE%94Tm_Dataset/S571.csv">S571.csv</a></td>
    </tr>
    <tr>
      <td>ArchStabMS1E10 Epistasis Dataset</td>
      <td><a href="https://doi.org/10.1038/s41586-024-07966-0">The genetic architecture of protein stability</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Dataset.csv">ArchStabMS1E10_Epistasis_Dataset.csv</a></td>
    </tr>
    <tr>
      <td>Antitoxin ParD3 Epistasis Dataset</td>
      <td><a href="https://doi.org/10.1038/s41467-024-45621-4">Protein design using structure-based residue preferences</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/Antitoxin_ParD3_Epistasis_Dataset/Antitoxin_ParD3_Epistasis_Dataset.csv">Antitoxin_ParD3_Epistasis_Dataset.csv</a></td>
    </tr>
    <tr>
      <td>TrpB Epistasis Dataset</td>
      <td><a href="https://doi.org/10.1073/pnas.2400439121">A combinatorially complete epistatic fitness landscape in an enzyme active site</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/TrpB_Epistasis_Dataset/TrpB_Epistasis_Dataset.csv">TrpB_Epistasis_Dataset.csv</a></td>
    </tr>
    <tr>
      <td>Human Myoglobin Epistasis Dataset</td>
      <td><a href="https://doi.org/10.1101/2024.02.24.581358">Decoding Stability and Epistasis in Human Myoglobin by Deep Mutational Scanning and Codon-level Machine Learning</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/Human_Myoglobin_Epistasis_Dataset/Human_Myoglobin_Epistasis_Dataset.csv">Human_Myoglobin_Epistasis_Dataset.csv</a></td>
    </tr>
    <tr>
      <td>CTXM Epistasis Dataset</td>
      <td><a href="https://doi.org/10.1073/pnas.2313513121">Network of epistatic interactions in an enzyme active site revealed by DMS</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/CTXM_Epistasis_Dataset/CTXM_Cefotaxime_Epistasis_Dataset.csv">CTXM_Cefotaxime_Epistasis_Dataset.csv</a>, <a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/CTXM_Epistasis_Dataset/CTXM_Ampicillin_Epistasis_Dataset.csv">CTXM_Ampicillin_Epistasis_Dataset.csv</a></td>
    </tr>
    <tr>
      <td rowspan="4" valign="middle">RBD ACE2 Dataset</td>
      <td><a href="https://doi.org/10.1126/science.abo7896">Shifting mutational constraints in the SARS-CoV-2 receptor-binding domain during viral evolution</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_ACE2_Dataset/SARS-CoV-2-RBD_DMS_variants_bc_binding.csv">SARS-CoV-2-RBD_DMS_variants_bc_binding.csv</a><br><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_ACE2_Dataset/SARS-CoV-2-RBD_Delta_bc_binding.csv">SARS-CoV-2-RBD_Delta_bc_binding.csv</a></td>
    </tr>
    <tr>
      <td><a href="https://doi.org/10.1371/journal.ppat.1010951">Deep mutational scans for ACE2 binding, RBD expression, and antibody escape in the SARS-CoV-2 Omicron BA.1 and BA.2 receptor-binding domains</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_ACE2_Dataset/SARS-CoV-2-RBD_DMS_Omicron_bc_binding.csv">SARS-CoV-2-RBD_DMS_Omicron_bc_binding.csv</a></td>
    </tr>
    <tr>
      <td><a href="https://doi.org/10.1371/journal.ppat.1011901">Deep mutational scans of XBB.1.5 and BQ.1.1 reveal ongoing epistatic drift during SARS-CoV-2 evolution</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_ACE2_Dataset/SARS-CoV-2-RBD_DMS_Omicron-XBB-BQ_bc_binding.csv">SARS-CoV-2-RBD_DMS_Omicron-XBB-BQ_bc_binding.csv</a></td>
    </tr>
    <tr>
      <td><a href="https://doi.org/10.1093/ve/veae067">Deep mutational scanning of SARS-CoV-2 Omicron BA.2.86 and epistatic emergence of the KP.3 variant</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_ACE2_Dataset/SARS-CoV-2-RBD_DMS_Omicron-EG5-FLip-BA286_bc_binding.csv">SARS-CoV-2-RBD_DMS_Omicron-EG5-FLip-BA286_bc_binding.csv</a></td>
    </tr>
    <tr>
      <td rowspan="7" valign="middle">RBD Antibody Dataset</td>
      <td><a href="https://doi.org/10.1126/scitranslmed.abi9915">The SARS-CoV-2 mRNA-1273 vaccine elicits more RBD-focused neutralization, but with broader antibody binding within the RBD</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_Antibody_Dataset/SARS-CoV-2-RBD_MAP_Moderna.csv">SARS-CoV-2-RBD_MAP_Moderna.csv</a></td>
    </tr>
    <tr>
      <td><a href="https://doi.org/10.1016/j.chom.2021.02.003">Comprehensive mapping of mutations in the SARS-CoV-2 receptor-binding domain that affect recognition by polyclonal human plasma antibodies</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_Antibody_Dataset/SARS-CoV-2-RBD_MAP_HAARVI_sera.csv">SARS-CoV-2-RBD_MAP_HAARVI_sera.csv</a></td>
    </tr>
    <tr>
      <td><a href="https://doi.org/10.1038/s41467-021-24435-8">Mapping mutations to the SARS-CoV-2 RBD that escape binding by different classes of antibodies</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_Antibody_Dataset/SARS-CoV-2-RBD_MAP_Rockefeller.csv">SARS-CoV-2-RBD_MAP_Rockefeller.csv</a></td>
    </tr>
    <tr>
      <td><a href="https://doi.org/10.1038/s41564-021-00972-2">Genetic and structural basis for SARS-CoV-2 variant neutralization by a two-antibody cocktail</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_Antibody_Dataset/SARS-CoV-2-RBD_MAP_AZ_Abs.csv">SARS-CoV-2-RBD_MAP_AZ_Abs.csv</a></td>
    </tr>
    <tr>
      <td><a href="https://doi.org/10.1038/s41586-021-03807-6">SARS-CoV-2 RBD antibodies that maximize breadth and resistance to escape</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_Antibody_Dataset/SARS-CoV-2-RBD_MAP_Vir_mAbs.csv">SARS-CoV-2-RBD_MAP_Vir_mAbs.csv</a></td>
    </tr>
    <tr>
      <td><a href="https://doi.org/10.1126/science.abf9302">Prospective mapping of viral mutations that escape antibodies used to treat COVID-19</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_Antibody_Dataset/SARS-CoV-2-RBD_MAP_clinical_Abs.csv">SARS-CoV-2-RBD_MAP_clinical_Abs.csv</a></td>
    </tr>
  </tbody>
</table>

### Processing cDNA Proteolysis Dataset

Here's a complete example demonstrating MutCleaner's capabilities with the cDNA Proteolysis mutation dataset:

```python
from mutcleaner import cdna_proteolysis_cleaner
from mutcleaner import download_cdna_proteolysis_source_file

# Create a cDNA Proteolysis cleaning pipeline using MutCleaner's default pipeline.
cdna_proteolysis_filepath = download_cdna_proteolysis_source_file(
    "dir_path",
    "file_name",
)["filename"]

cdna_proteolysis_cleaning_pipeline = cdna_proteolysis_cleaner.create_cdna_proteolysis_cleaner(
    cdna_proteolysis_filepath,
)

# Clean and process the dataset.
cdna_proteolysis_cleaning_pipeline, cdna_proteolysis_dataset = (
    cdna_proteolysis_cleaner.clean_cdna_proteolysis_dataset(
        cdna_proteolysis_cleaning_pipeline,
    )
)

# Save the processed dataset.
cdna_proteolysis_dataset.save("output/cleaned_cdna_proteolysis_data")
```

### Basic Sequence Operations

```python
from mutcleaner.core.sequence import DNASequence

# DNA sequence analysis.
dna = DNASequence("ATGCGATCGTAA")

print(f"Reverse complement: {dna.reverse_complement()}")
print(f"Transcription: {dna.transcribe()}")
print(f"Translation: {dna.translate()}")
```

## Core Features

### Sequence Data Manipulation

- **Sequence validation**: Validate DNA, RNA, and protein sequences against predefined alphabets.
- **Sequence transformation**: Support transcription, reverse transcription, translation, and reverse-complement operations.
- **Batch processing**: Process large tabular mutation datasets through reusable cleaning utilities.

### Mutation Analysis

- **Mutation parsing**: Parse amino-acid and codon-level mutation annotations.
- **Mutation inference**: Infer mutation annotations by comparing reference and mutated sequences.
- **Mutation transformation**: Apply mutation annotations to reference sequences and convert codon-level mutations into amino-acid changes.

### Data Cleaning and Preprocessing

- **Standardization**: Harmonize mutation names, sequences, labels, and metadata across heterogeneous datasets.
- **Duplicate handling**: Remove or aggregate redundant mutation records according to dataset-specific rules.
- **Dataset-specific cleaners**: Provide reusable cleaning pipelines for commonly used mutation datasets.

### Pipeline Architecture

- **Modular design**: Compose cleaning workflows from reusable processing components.
- **Parallel processing**: Use multi-core processing for mutation validation, mutation application, and sequence-based mutation inference.
- **Progress tracking**: Monitor long-running cleaning tasks with progress bars and structured execution summaries.

## Examples and Use Cases

### Custom Processing Pipeline

```python
import pandas as pd

from mutcleaner.cleaners.basic_cleaners import (
    extract_and_rename_columns,
    filter_and_clean_data,
    convert_data_types,
    validate_mutations,
    convert_to_mutation_dataset_format,
)
from mutcleaner.cleaners.cdna_proteolysis_custom_cleaners import (
    validate_wt_sequence,
    average_labels_by_name,
    subtract_labels_by_wt,
)
from mutcleaner.core.dataset import MutationDataset
from mutcleaner.core.pipeline import create_pipeline

dataset = pd.read_csv("path/to/Tsuboyama2023_Dataset2_Dataset3_20230416.csv")

pipeline = create_pipeline(dataset, "cdna_proteolysis_cleaner")
clean_result = (
    pipeline.then(
        extract_and_rename_columns,
        column_mapping={
            "WT_name": "name",
            "aa_seq": "mut_seq",
            "mut_type": "mut_info",
            "ddG_ML": "ddG",
        },
    )
    .then(filter_and_clean_data, filters={"ddG": lambda x: x != "-"})
    .then(convert_data_types, type_conversions={"ddG": "float"})
    .then(
        validate_mutations,
        mutation_column="mut_info",
        mutation_sep="_",
        is_zero_based=False,
        num_workers=16,
    )
    .then(
        average_labels_by_name,
        name_columns=("name", "mut_info"),
        label_columns="ddG",
    )
    .then(
        validate_wt_sequence,
        name_column="name",
        mutation_column="mut_info",
        sequence_column="mut_seq",
        wt_identifier="wt",
        num_workers=16,
    )
    .then(
        subtract_labels_by_wt,
        name_column="name",
        label_columns="ddG",
        mutation_column="mut_info",
        in_place=True,
    )
    .then(
        convert_to_mutation_dataset_format,
        name_column="name",
        mutation_column="mut_info",
        mutated_sequence_column="mut_seq",
        score_column="ddG",
        is_zero_based=True,
    )
)

cdna_proteolysis_dataset_df, cdna_proteolysis_ref_seq = clean_result.data
cdna_proteolysis_dataset = MutationDataset.from_dataframe(
    cdna_proteolysis_dataset_df,
    cdna_proteolysis_ref_seq,
)

execution_info = pipeline.get_execution_summary()
artifacts = pipeline.artifacts
pipeline.save_structured_data("cdna_proteolysis_cleaner_pipeline.pkl")
```

## Citation

If you use MutCleaner in your research, please cite:

```bibtex
@software{mutcleaner,
  title={
    MutCleaner: An efficient framework for cleaning, standardizing, and processing biological mutation data.
  },
  author={Yuxiang Tang and Ziyu Shi},
  year={2026},
  url={https://github.com/xulab-research/MutCleaner}
}
```

## License

This project is licensed under the BSD 3-Clause License. See the [LICENSE](LICENSE) file for details.

## Support

- **Issues**: [GitHub Issues](https://github.com/xulab-research/MutCleaner/issues)
- **Discussions**: [GitHub Discussions](https://github.com/xulab-research/MutCleaner/discussions)
