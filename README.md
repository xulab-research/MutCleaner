# MutCleaner

[![License badge](https://img.shields.io/badge/License-BSD_3--Clause-yellow?logo=opensourceinitiative&logoColor=white)](https://github.com/xulab-research/MutCleaner/blob/main/LICENSE)
[![PyPI version badge](https://img.shields.io/pypi/v/mutcleaner.svg?logo=python&logoColor=white&color=orange)](https://pypi.org/project/mutcleaner/)
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

## Package Structure
```text
mutcleaner/
├── cleaners/                # Reusable cleaners and dataset-specific pipelines
│   ├── basic_cleaners.py    # Reusable basic cleaning functions
│   ├── base_config.py       # Shared cleaner configuration
│   └── *_cleaner.py         # Dataset-specific cleaning pipelines
├── core/                    # Core data structures and processing logic
│   ├── alphabet.py          # Biological alphabet definitions and validation
│   ├── codon.py             # Codon table and sequence translation utilities
│   ├── constants.py         # Shared biological constants
│   ├── dataset.py           # Standard dataset abstraction and export logic
│   ├── mutation.py          # Mutation parsing, representation and validation
│   ├── pipeline.py          # Reusable data cleaning pipeline framework
│   ├── sequence.py          # DNA, RNA and protein sequence abstractions
│   └── types.py             # Shared type aliases and annotations
└── utils/                   # General helper functions for downloading and I/O
```
### Module Overview
#### mutcleaner.cleaners

The `cleaners` module provides both reusable cleaning functions and predefined dataset-specific cleaning pipelines.  

`basic_cleaners.py` contains general-purpose cleaning functions that can be reused across different datasets, such as column checking, missing value handling, sequence validation, mutation validation, and common formatting operations.  
`base_config.py` defines shared configuration used by dataset cleaners, such as common column names, required fields, and reusable cleaner settings.  
`*_cleaner.py` files define dataset-specific cleaning pipelines. Each file is designed for a particular dataset or benchmark and combines reusable cleaning functions into a complete workflow.  
#### mutcleaner.core

The `core` module contains the fundamental data structures and processing logic used throughout MutCleaner.  

`alphabet.py` defines biological alphabets and validation rules for DNA, RNA, and protein sequences.  
`codon.py` provides codon table utilities and sequence translation functionality.
`constants.py` stores shared biological constants, such as amino acid symbols, nucleotide symbols, complements, and stop codon definitions.  
`dataset.py` defines the standard dataset abstraction used to store, manage, validate, and export cleaned mutation datasets.  
`mutation.py` provides mutation parsing, mutation representation, and mutation validation logic.  
`pipeline.py` defines the reusable data cleaning pipeline framework, including pipeline steps, execution order, and artifact tracking.  
`sequence.py` defines biological sequence abstractions, including DNA, RNA, and protein sequence classes.  
`types.py` stores shared type aliases and annotations used across the package.  

#### mutcleaner.utils

The `utils` module contains helper functions that support common operations outside the core cleaning logic.  
It includes utilities for downloading source files, handling paths, extracting files, checking file existence, and managing common input/output operations.  
## Quick Start

See the [Data Cleaners Usage Guide](https://xulab-research.github.io/MutCleaner/user_guide/cleaners.html) for more examples.

### Supported Datasets

| Dataset Name | Reference | File |
| --- | --- | --- |
| Human Domainome Dataset | [Site-saturation mutagenesis of 500 human protein domains](https://doi.org/10.1038/s41586-024-08370-4) | [SupplementaryTable2.txt](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/Human_Domainome_Dataset/SupplementaryTable2.txt), [wild_type.fasta](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/Human_Domainome_Dataset/wild_type.fasta) |
| ProteinGym DMS Substitutions Dataset | [ProteinGym: Large-Scale Benchmarks for Protein Design and Fitness Prediction](https://doi.org/10.1101/2023.12.07.570727) | [DMS_ProteinGym_substitutions.zip](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/ProteinGym_DMS_Substitutions_Dataset/DMS_ProteinGym_substitutions.zip) |
| cDNA Proteolysis Dataset | [Mega-scale experimental analysis of protein folding stability in biology and design](https://doi.org/10.1038/s41586-023-06328-6) | [Tsuboyama2023_Dataset2_Dataset3_20230416.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/cDNA_Proteolysis_Dataset/Tsuboyama2023_Dataset2_Dataset3_20230416.csv) |
| ΔΔG Dataset | [Improving the prediction of protein stability changes upon mutations by geometric learning and a pre-training strategy](https://doi.org/10.1038/s43588-024-00716-2) | [M1261.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/%CE%94%CE%94G_Dataset/M1261.csv), [S461.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/%CE%94%CE%94G_Dataset/S461.csv), [S669.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/%CE%94%CE%94G_Dataset/S669.csv), [S783.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/%CE%94%CE%94G_Dataset/S783.csv), [S8754.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/%CE%94%CE%94G_Dataset/S8754.csv) |
| ΔTm Dataset | [Improving the prediction of protein stability changes upon mutations by geometric learning and a pre-training strategy](https://doi.org/10.1038/s43588-024-00716-2) | [S4346.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/%CE%94Tm_Dataset/S4346.csv), [S571.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/%CE%94Tm_Dataset/S571.csv) |
| ArchStabMS1E10 Epistasis Dataset | [The genetic architecture of protein stability](https://doi.org/10.1038/s41586-024-07966-0) | [ArchStabMS1E10_Epistasis_Dataset.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Dataset.csv) |
| Antitoxin ParD3 Epistasis Dataset | [Protein design using structure-based residue preferences](https://doi.org/10.1038/s41467-024-45621-4) | [Antitoxin_ParD3_Epistasis_Dataset.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/Antitoxin_ParD3_Epistasis_Dataset/Antitoxin_ParD3_Epistasis_Dataset.csv) |
| TrpB Epistasis Dataset | [A combinatorially complete epistatic fitness landscape in an enzyme active site](https://doi.org/10.1073/pnas.2400439121) | [TrpB_Epistasis_Dataset.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/TrpB_Epistasis_Dataset/TrpB_Epistasis_Dataset.csv) |
| Human Myoglobin Epistasis Dataset | [Decoding Stability and Epistasis in Human Myoglobin by Deep Mutational Scanning and Codon-level Machine Learning](https://doi.org/10.1101/2024.02.24.581358) | [Human_Myoglobin_Epistasis_Dataset.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/Human_Myoglobin_Epistasis_Dataset/Human_Myoglobin_Epistasis_Dataset.csv) |
| CTXM Epistasis Dataset | [Network of epistatic interactions in an enzyme active site revealed by DMS](https://doi.org/10.1073/pnas.2313513121) | [CTXM_Cefotaxime_Epistasis_Dataset.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/CTXM_Epistasis_Dataset/CTXM_Cefotaxime_Epistasis_Dataset.csv), [CTXM_Ampicillin_Epistasis_Dataset.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/CTXM_Epistasis_Dataset/CTXM_Ampicillin_Epistasis_Dataset.csv) |
| RBD ACE2 Dataset | [Shifting mutational constraints in the SARS-CoV-2 receptor-binding domain during viral evolution](https://doi.org/10.1126/science.abo7896), [Deep mutational scans for ACE2 binding, RBD expression, and antibody escape in the SARS-CoV-2 Omicron BA.1 and BA.2 receptor-binding domains](https://doi.org/10.1371/journal.ppat.1010951), [Deep mutational scans of XBB.1.5 and BQ.1.1 reveal ongoing epistatic drift during SARS-CoV-2 evolution](https://doi.org/10.1371/journal.ppat.1011901), [Deep mutational scanning of SARS-CoV-2 Omicron BA.2.86 and epistatic emergence of the KP.3 variant](https://doi.org/10.1093/ve/veae067) | [SARS-CoV-2-RBD_DMS_variants_bc_binding.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_ACE2_Dataset/SARS-CoV-2-RBD_DMS_variants_bc_binding.csv), [SARS-CoV-2-RBD_Delta_bc_binding.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_ACE2_Dataset/SARS-CoV-2-RBD_Delta_bc_binding.csv), [SARS-CoV-2-RBD_DMS_Omicron_bc_binding.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_ACE2_Dataset/SARS-CoV-2-RBD_DMS_Omicron_bc_binding.csv), [SARS-CoV-2-RBD_DMS_Omicron-XBB-BQ_bc_binding.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_ACE2_Dataset/SARS-CoV-2-RBD_DMS_Omicron-XBB-BQ_bc_binding.csv), [SARS-CoV-2-RBD_DMS_Omicron-EG5-FLip-BA286_bc_binding.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_ACE2_Dataset/SARS-CoV-2-RBD_DMS_Omicron-EG5-FLip-BA286_bc_binding.csv) |
| RBD Antibody Dataset | [Antibodies elicited by mRNA-1273 vaccination bind more broadly to the receptor binding domain than do those from SARS-CoV-2 infection](https://doi.org/10.1126/scitranslmed.abi9915), [Comprehensive mapping of mutations in the SARS-CoV-2 receptor-binding domain that affect recognition by polyclonal human plasma antibodies](https://doi.org/10.1016/j.chom.2021.02.003), [Mapping mutations to the SARS-CoV-2 RBD that escape binding by different classes of antibodies](https://doi.org/10.1038/s41467-021-24435-8), [Genetic and structural basis for SARS-CoV-2 variant neutralization by a two-antibody cocktail](https://doi.org/10.1038/s41564-021-00972-2), [SARS-CoV-2 RBD antibodies that maximize breadth and resistance to escape](https://doi.org/10.1038/s41586-021-03807-6), [Prospective mapping of viral mutations that escape antibodies used to treat COVID-19](https://doi.org/10.1126/science.abf9302) | [SARS-CoV-2-RBD_MAP_Moderna_scores.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_Antibody_Dataset/SARS-CoV-2-RBD_MAP_Moderna_scores.csv), [SARS-CoV-2-RBD_MAP_HAARVI_sera_scores.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_Antibody_Dataset/SARS-CoV-2-RBD_MAP_HAARVI_sera_scores.csv), [SARS-CoV-2-RBD_MAP_Rockefeller_scores.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_Antibody_Dataset/SARS-CoV-2-RBD_MAP_Rockefeller_scores.csv), [SARS-CoV-2-RBD_MAP_AZ_Abs_scores.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_Antibody_Dataset/SARS-CoV-2-RBD_MAP_AZ_Abs_scores.csv), [SARS-CoV-2-RBD_MAP_Vir_mAbs_scores.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_Antibody_Dataset/SARS-CoV-2-RBD_MAP_Vir_mAbs_scores.csv), [SARS-CoV-2-RBD_MAP_clinical_Abs_scores.csv](https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_Antibody_Dataset/SARS-CoV-2-RBD_MAP_clinical_Abs_scores.csv) |
| Chitosanase Dataset |  |  |

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


## Citation

If you use MutCleaner in your research, please cite:

```bibtex
@software{mutcleaner,
  title={
    MutCleaner: An efficient framework for cleaning, standardizing, and processing biological mutation data.
  },
  author={Ziyu Shi and Yuxiang Tang},
  year={2026},
  url={https://github.com/xulab-research/MutCleaner}
}
```

## License

This project is licensed under the BSD 3-Clause License. See the [LICENSE](LICENSE) file for details.

## Support

- **Issues**: [GitHub Issues](https://github.com/xulab-research/MutCleaner/issues)
- **Discussions**: [GitHub Discussions](https://github.com/xulab-research/MutCleaner/discussions)
