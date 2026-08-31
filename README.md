# MutCleaner

[![PyPI version badge](https://img.shields.io/pypi/v/mutcleaner.svg?logo=python&logoColor=white&color=orange)](https://pypi.org/project/mutcleaner/)
[![Python version badge](https://img.shields.io/python/required-version-toml.svg?logo=python&logoColor=white&color=orange&tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Fxulab-research%2FMutCleaner%2Frefs%2Fheads%2Fmain%2Fpyproject.toml)](https://pypi.org/project/mutcleaner/)
[![License badge](https://img.shields.io/badge/License-Apache_2.0-blue?logo=apache&logoColor=white)](https://github.com/xulab-research/MutCleaner/blob/main/LICENSE)
[![Hugging Face](https://img.shields.io/badge/HuggingFace-datasets-yellow)](https://huggingface.co/datasets/xulab-research/MutCleaner)
[![Docs](https://github.com/xulab-research/MutCleaner/actions/workflows/docs.yml/badge.svg)](https://xulab-research.github.io/MutCleaner/)

MutCleaner is an extensible Python framework that cleans, validates, and standardizes protein- and codon-level mutation datasets through composable cleaning pipelines, unified sequence and mutation data structures, and dataset-specific cleaners.

* **Documentation**: https://xulab-research.github.io/MutCleaner
* **Cleaning Examples**: https://xulab-research.github.io/MutCleaner/user_guide/supported_datasets.html

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
│   ├── basic_cleaners.py    # Reusable data-cleaning and standardization steps
│   ├── base_config.py       # Base configuration for dataset-specific cleaners
│   └── ...                  # Dataset-specific cleaner modules
├── core/                    # Core data structures and processing logic
│   ├── alphabet.py          # Biological alphabets and sequence validation
│   ├── codon.py             # Codon-table definitions and translation
│   ├── constants.py         # Biological alphabets, mappings and genetic-code constants
│   ├── dataset.py           # MutationDataset representation, validation and export
│   ├── mutation.py          # Mutation parsing, representation and validation
│   ├── pipeline.py          # Composable data-cleaning pipeline framework
│   ├── sequence.py          # DNA, RNA and protein sequence representations and operations
│   └── types.py             # Shared type variables and type annotations
└── utils/                   # Supporting utilities for conversion, processing and I/O
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

See the [Data Cleaners Usage Guide](https://xulab-research.github.io/MutCleaner/user_guide) for more examples.

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
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/Human_Domainome_Dataset/SupplementaryTable2.txt">SupplementaryTable2.txt</a>,<a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/Human_Domainome_Dataset/SupplementaryTable4.txt">SupplementaryTable4.txt</a>,<a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/Human_Domainome_Dataset/wild_type.fasta">wild_type.fasta</a></td>
    </tr>
    <tr>
      <td>ProteinGym DMS Substitutions Dataset</td>
      <td><a href="https://doi.org/10.1101/2023.12.07.570727">ProteinGym: Large-Scale Benchmarks for Protein Design and Fitness Prediction</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/ProteinGym_DMS_Substitutions_Dataset/DMS_ProteinGym_substitutions.zip">DMS_ProteinGym_substitutions.zip</a></td>
    </tr>
    <tr>
      <td>Protein cDNA Proteolysis Dataset</td>
      <td><a href="https://doi.org/10.1038/s41586-023-06328-6">Mega-scale experimental analysis of protein folding stability in biology and design</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/cDNA_Proteolysis_Dataset/Tsuboyama2023_Dataset2_Dataset3_20230416.csv">Tsuboyama2023_Dataset2_Dataset3_20230416.csv</a></td>
    </tr>
    <tr>
      <td>ddG Dataset</td>
      <td><a href="https://doi.org/10.1038/s43588-024-00716-2">Improving the prediction of protein stability changes upon mutations by geometric learning and a pre-training strategy</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/resolve/main/ddG_Dataset/M1261.csv">M1261.csv</a>, <a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/ddG_Dataset/S461.csv">S461.csv</a>, <a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/ddG_Dataset/S669.csv">S669.csv</a>, <a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/ddG_Dataset/S783.csv">S783.csv</a>, <a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/ddG_Dataset/S8754.csv">S8754.csv</a></td>
    </tr>
    <tr>
      <td>dTm Dataset</td>
      <td><a href="https://doi.org/10.1038/s43588-024-00716-2">Improving the prediction of protein stability changes upon mutations by geometric learning and a pre-training strategy</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/dTm_Dataset/S4346.csv">S4346.csv</a>, <a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/dTm_Dataset/S557.csv">S557.csv</a></td>
    </tr>
    <tr>
      <td>ArchStabMS1E10 Epistasis Dataset</td>
      <td><a href="https://doi.org/10.1038/s41586-024-07966-0">The genetic architecture of protein stability</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Sup4_Dataset.csv">ArchStabMS1E10_Epistasis_Sup4_Dataset.csv</a>,<a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Sup5_Dataset.csv">ArchStabMS1E10_Epistasis_Sup5_Dataset.csv</a></td>
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
      <td>Protein Human Myoglobin Epistasis Dataset</td>
      <td><a href="https://doi.org/10.1101/2024.02.24.581358">Decoding Stability and Epistasis in Human Myoglobin by Deep Mutational Scanning and Codon-level Machine Learning</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/Protein_Human_Myoglobin_Epistasis_Dataset/Protein_Human_Myoglobin_Epistasis_Dataset.csv">Protein_Human_Myoglobin_Epistasis_Dataset.csv</a></td>
    </tr>
    <tr>
      <td>CTXM Epistasis Dataset</td>
      <td><a href="https://doi.org/10.1073/pnas.2313513121">Network of epistatic interactions in an enzyme active site revealed by DMS</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/CTXM_Epistasis_Dataset/CTXM_Cefotaxime_Epistasis_Dataset.csv">CTXM_Cefotaxime_Epistasis_Dataset.csv</a>, <a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/CTXM_Epistasis_Dataset/CTXM_Ampicillin_Epistasis_Dataset.csv">CTXM_Ampicillin_Epistasis_Dataset.csv</a></td>
    </tr>
    <tr>
      <td rowspan="4" valign="middle">RBD ACE2 Dataset</td>
      <td><a href="https://doi.org/10.1126/science.abo7896">Shifting mutational constraints in the SARS-CoV-2 receptor-binding domain during viral evolution</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_ACE2_Dataset/SARS-CoV-2-RBD_DMS_variants_bc_binding.csv">SARS-CoV-2-RBD_DMS_variants_bc_binding.csv</a>,<a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_ACE2_Dataset/SARS-CoV-2-RBD_Delta_bc_binding.csv">SARS-CoV-2-RBD_Delta_bc_binding.csv</a></td>
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
      <td rowspan="3" valign="middle">RBD Antibody Dataset</td>
      <td><a href="https://doi.org/10.1126/scitranslmed.abi9915">The SARS-CoV-2 mRNA-1273 vaccine elicits more RBD-focused neutralization, but with broader antibody binding within the RBD</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_Antibody_Dataset/SARS-CoV-2-RBD_MAP_Moderna.csv">SARS-CoV-2-RBD_MAP_Moderna.csv</a></td>
    </tr>
    <tr>
      <td><a href="https://doi.org/10.1038/s41467-021-24435-8">Mapping mutations to the SARS-CoV-2 RBD that escape binding by different classes of antibodies</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_Antibody_Dataset/SARS-CoV-2-RBD_MAP_Rockefeller.csv">SARS-CoV-2-RBD_MAP_Rockefeller.csv</a></td>
    </tr>
    <tr>
      <td><a href="https://doi.org/10.1038/s41586-021-03807-6">SARS-CoV-2 RBD antibodies that maximize breadth and resistance to escape</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/RBD_Antibody_Dataset/SARS-CoV-2-RBD_MAP_Vir_mAbs.csv">SARS-CoV-2-RBD_MAP_Vir_mAbs.csv</a></td>
    </tr>
      <td>Chitosanase dTm Dataset</td>
      <td>In-house wet-lab data, no reference available yet</td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/Chitosanase_dTm_Dataset/Chitosanase_dTm_Dataset.csv">Chitosanase_dTm_Dataset.csv</a></td>
    </tr>
    <tr>
      <td>MGnify ddG Dataset</td>
      <td><a href="https://doi.org/10.64898/2026.05.19.726285">Accurate protein stability prediction for small domains using mega-scale experiments</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/MGnify_ddG_Dataset/MGnify_ddG_Dataset.csv">MGnify_ddG_Dataset.csv</a></td>
    </tr>
    </tr>
      <td>Codon cDNA Proteolysis Dataset</td>
      <td><a href="https://doi.org/10.1038/s41586-023-06328-6">Mega-scale experimental analysis of protein folding stability in biology and design</td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/Codon_cDNA_Proteolysis_Dataset/Codon_cDNA_Proteolysis_Dataset.csv">Codon_cDNA_Proteolysis_Dataset.csv</a>,<a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/Codon_cDNA_Proteolysis_Dataset/wt.fasta">wt.fasta</a>
      </td>
    </tr>
    </tr>
      <td>Codon DMS Substitutions Dataset</td>
      <td><a href="https://doi.org/10.1186/s13059-025-03476-y">MaveDB 2024: a curated community database with over seven million variant effects from multiplexed functional assays</td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/
      Codon_DMS_Substitutions_Dataset/Codon_DMS_Substitutions_Datasetzip">Codon_DMS_Substitutions_Dataset.zip</a>
      </td>
    </tr>
    <tr>
      <td>Codon Human Myoglobin Epistasis Dataset</td>
      <td><a href="https://doi.org/10.1101/2024.02.24.581358">Decoding Stability and Epistasis in Human Myoglobin by Deep Mutational Scanning and Codon-level Machine Learning</a></td>
      <td><a href="https://huggingface.co/datasets/xulab-research/MutCleaner/blob/main/Codon_Human_Myoglobin_Epistasis_Dataset/Codon_Human_Myoglobin_Epistasis_Dataset.csv">Codon_Human_Myoglobin_Epistasis_Dataset.csv</a></td>
    </tr>
  </tbody>
</table>


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


## Citation

If you use MutCleaner in your research, please cite:

```bibtex
@misc{mutcleaner,
  title = {MutCleaner: Cleaning and Standardizing Biological Mutation Datasets for Variant Effect Prediction},
  author = {Shi, Ziyu and Tang, Yuxiang and Yang, Mengxin and Shi, Yancheng and Yu, Shize and Xu, Yunxin},
  year = {2026},
  url = {https://github.com/xulab-research/MutCleaner}
}
```

## License

This project is licensed under the Apache License 2.0.

Unless otherwise stated, the source code, model architecture, training scripts,
inference scripts, and released model weights/checkpoints are licensed under
Apache-2.0.

Datasets used in this project may be subject to their original licenses and
terms of use. Please refer to the corresponding dataset sources for details.

This software is provided for research purposes and is not intended for clinical
diagnosis, medical decision-making, or direct therapeutic use.
