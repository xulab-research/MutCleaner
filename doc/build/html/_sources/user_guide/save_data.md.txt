# Data Saving Usage Guide

## Overview

This guide provides usage examples for saving cleaned datasets and cleaning artifacts, organized by database:

- [**HumanDomainome**](#human-domainome-database): Site-saturation mutagenesis of 500 human protein domains
    - [**Supplementary Table 2**](#supplementarytable2-cleaner): Fitness scores and errors.
    - [**Supplementary Table 4**](#supplementarytable4-cleaner): Homolog-averaged ∆∆G predictions across families mapped to homologous domains proteome-wide.
- [**ProteinGym**](#protein-gym-database): ProteinGym: Large-Scale Benchmarks for Protein Design and Fitness Prediction
- [**cDNAProteolysis**](#cdna-proteolysis-database): Mega-scale experimental analysis of protein folding stability in biology and design
- [**ddG-dTm Datasets**](#ddg-dtm-datasets): A collection of datasets providing single- and multiple-mutant measurements, labeled by thermodynamic parameters (ΔΔG, ΔTm)
- [**ArchStabMS1E10 Datasets**](#archstabms1e10-datasets): High-order multi-mutant libraries (“1e10”) measuring protein stability for GRB2-SH3 and SRC.
- [**Antitoxin ParD3 Datasets**](#antitoxin-pard3): The antitoxin ParD3 3-position library is a combinatorially exhaustive dataset of 8,000 variants demonstrating that simple, independent per-residue mutation preferences are sufficient to almost perfectly predict combinatorial protein fitness.
- [**TrpB Datasets**](#trpb-datasets): a combinatorially complete sequence-fitness landscape comprising 160,000 variants across four active-site residues of the enzyme tryptophan synthase, capturing significant epistatic interactions to serve as a benchmark for model-guided enzyme engineering.
- [**Human Myoglobin Datasets**](#human-myoglobin-datasets): a deep mutational scanning library detailing the expression fitness scores for near-comprehensive single-codon and small-fraction double-codon mutations in yeast surface-displayed human myoglobin, which was used to train machine learning models for predicting epistatic effects and discovering stability-enhancing variants.
- [**CTXM Datasets**](#ctxm-database): a comprehensive deep mutational scanning library of 49,096 pairwise double mutations across 17 active site residues of the CTX-M-14 $\beta$-lactamase enzyme, constructed to systematically map the epistatic interaction network driving antibiotic resistance.
    - [**CTXM ampicillin**]: A subset of the CTX-M library quantifying the functional fitness and epistatic interactions of the enzyme variants under ampicillin selection, revealing a broader mutational tolerance and distinct compensatory pathways.
    - [**CTXM cefotaxime**]: A subset of the CTX-M library quantifying the functional fitness and epistatic interactions of the enzyme variants under cefotaxime selection, characterized by highly stringent sequence requirements and substrate-specific epistasis.
- [**RBD ACE2 Database**](#rbd-ace2-database): SARS-CoV-2 RBD sequences with ACE2 binding affinity scores, labeled by `log10Ka` where higher values indicate stronger ACE2 binding affinity.
- [**RBD Antibody Database**](#rbd-antibody-database): SARS-CoV-2 RBD antibody escape data with `score` computed as the negative logarithm of escape. Higher scores indicate weaker escape, reflecting better binding capacity.

### Main idea

- **Main dataset** → use `dataset.save(path)`  
  → stored in the **mutcleaner format** (default `save_type="mutcleaner"`), which can be reloaded directly by mutcleaner.
- **Filtered / discarded / intermediate rows** (“artifacts”) → use `cleaning_pipeline.save_artifacts(path)`  
  → saved as a single **pickle** file (`.pkl`) containing multiple DataFrames.  
  → since pickle is binary and not convenient to view in VS Code, we **recommend reading it back in Python and exporting individual tables to CSV**.

## Prerequisites

You have already cleaned the raw data using the corresponding cleaning pipeline  
(see the *Data Cleaners Usage Guide*).

---

## Common Patterns

### Save the main cleaned dataset (mutcleaner format)

```python
# Default: save_type="mutcleaner"
dataset.save("path/to/output_dir")
```

### Save artifacts (filtered / discarded data) as pickle, then export to CSV

```python
# Save all artifacts (e.g., filtered-out rows, QC tables, per-step outputs)
cleaning_pipeline.save_artifacts("path/to/artifacts.pkl")

import pickle
import pandas as pd

# Load the binary pickle (not convenient to inspect directly in VS Code)
with open("path/to/artifacts.pkl", "rb") as f:
    artifacts = pickle.load(f)  # typically a dict: {name: DataFrame}

# Export each artifact DataFrame to a separate CSV for easy viewing
for key, df in artifacts.items():
    out_path = f"path/to/target_dir/{key}.csv"
    df.to_csv(out_path, index=False)
```

## Human Domainome Database

### SupplementaryTable2 Dataset

**save the main data**

```python
hd_dataset.save("path/to/output_dir")
```
**save the artifacts data**

```python
hd_cleaning_pipeline.save_artifacts("path/to/artifacts.pkl")

with open("path/to/artifacts.pkl", "rb") as f:
    obj = pickle.load(f)
    
for key, val in obj.items():
    out_path = f"path/to/target_dir/{key}.csv"
    val.to_csv(out_path, index=False)
```


### SupplementaryTable4 Dataset

**save the main data**

```python
hd_dataset.save("path/to/output_dir")
```
**save the artifacts data**

```python
hd_cleaning_pipeline.save_artifacts("path/to/artifacts.pkl")

with open("path/to/artifacts.pkl", "rb") as f:
    obj = pickle.load(f)
    
for key, val in obj.items():
    out_path = f"path/to/target_dir/{key}.csv"
    val.to_csv(out_path, index=False)
```

## cDNA Proteolysis Database

### ΔΔG Dataset

**save the main data**

```python
hd_dataset.save("path/to/output_dir")
```
**save the artifacts data**

```python
hd_cleaning_pipeline.save_artifacts("path/to/artifacts.pkl")

with open("path/to/artifacts.pkl", "rb") as f:
    obj = pickle.load(f)
    
for key, val in obj.items():
    out_path = f"path/to/target_dir/{key}.csv"
    val.to_csv(out_path, index=False)
```

### ΔG Dataset

**save the main data**

```python
hd_dataset.save("path/to/output_dir")
```
**save the artifacts data**

```python
hd_cleaning_pipeline.save_artifacts("path/to/artifacts.pkl")

with open("path/to/artifacts.pkl", "rb") as f:
    obj = pickle.load(f)
    
for key, val in obj.items():
    out_path = f"path/to/target_dir/{key}.csv"
    val.to_csv(out_path, index=False)
```

## ddG-dTm Datasets

**save the main data**

```python
ddgdtm_dataset.save("path/to/artifacts.pkl")
```
**save the artifacts data**

```python
ddgdtm_cleaning_pipeline.save_artifacts.save_artifacts("path/to/artifacts.pkl")

with open("path/to/artifacts.pkl", "rb") as f:
    obj = pickle.load(f)
    
for key, val in obj.items():
    out_path = f"path/to/target_dir/{key}.csv"
    val.to_csv(out_path, index=False)
```

## archstabms_1e10 Datasets

**save the main data**

```python
archstabms_dataset.save("path/to/output_dir")
```
**save the artifacts data**

```python
archstabms_cleaning_pipeline.save_artifacts("path/to/artifacts.pkl")

with open("path/to/artifacts.pkl", "rb") as f:
    obj = pickle.load(f)
    
for key, val in obj.items():
    out_path = f"path/to/target_dir/{key}.csv"
    val.to_csv(out_path, index=False)
```

## Antitoxin ParD3

**save the main data**

```python
antitoxin_pard3_dataset.save("path/to/output_dir")
```
**save the artifacts data**

```python
antitoxin_pard3_cleaning_pipeline.save_artifacts("path/to/artifacts.pkl")

with open("path/to/artifacts.pkl", "rb") as f:
    obj = pickle.load(f)
    
for key, val in obj.items():
    out_path = f"path/to/target_dir/{key}.csv"
    val.to_csv(out_path, index=False)
```


## Trpb Datasets

**save the main data**

```python
trpB_dataset.save("path/to/output_dir")
```
**save the artifacts data**

```python
trpB_cleaning_pipeline.save_artifacts("path/to/artifacts.pkl")

with open("path/to/artifacts.pkl", "rb") as f:
    obj = pickle.load(f)
    
for key, val in obj.items():
    out_path = f"path/to/target_dir/{key}.csv"
    val.to_csv(out_path, index=False)
```

## Human Myoglobin Datasets
**save the main data**

```python
human_myoglobin_dataset.save("path/to/output_dir")
```
**save the artifacts data**

```python
human_myoglobin_cleaning_pipeline.save_artifacts("path/to/artifacts.pkl")

with open("path/to/artifacts.pkl", "rb") as f:
    obj = pickle.load(f)
    
for key, val in obj.items():
    out_path = f"path/to/target_dir/{key}.csv"
    val.to_csv(out_path, index=False)
```

## CTXM Database

**save the main data**

```python
ctxm_dataset.save("path/to/output_dir")
```
**save the artifacts data**

```python
ctxm_cleaning_pipeline.save_artifacts("path/to/artifacts.pkl")

with open("path/to/artifacts.pkl", "rb") as f:
    obj = pickle.load(f)
    
for key, val in obj.items():
    out_path = f"path/to/target_dir/{key}.csv"
    val.to_csv(out_path, index=False)
```


## RBD ACE2 Database

**save the main data**

```python
rbd_ace2_dataset.save("path/to/output_dir")
```

**save the artifacts data**

```python
import pickle

rbd_ace2_cleaning_pipeline.save_artifacts("path/to/artifacts.pkl")

with open("path/to/artifacts.pkl", "rb") as f:
    obj = pickle.load(f)

for key, val in obj.items():
    out_path = f"path/to/target_dir/{key}.csv"
    val.to_csv(out_path, index=False)
```

## RBD Antibody Database

**save the main data**

```python
rbd_antibody_dataset.save("path/to/output_dir")
```

**save the artifacts data**

```python
import pickle

rbd_antibody_cleaning_pipeline.save_artifacts("path/to/artifacts.pkl")

with open("path/to/artifacts.pkl", "rb") as f:
    obj = pickle.load(f)

for key, val in obj.items():
    out_path = f"path/to/target_dir/{key}.csv"
    val.to_csv(out_path, index=False)
```
