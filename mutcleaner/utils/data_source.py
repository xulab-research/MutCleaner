# mutcleaner/utils/data_source.py
DATASETS = {
    "Protein cDNA Proteolysis Dataset": {
        "paper_title": "Mega-scale experimental analysis of protein folding stability in biology and design",
        "official_doi": "https://doi.org/10.1038/s41586-023-06328-6",
        "huggingface_repos": [
            "datasets/xulab-research/MutCleaner/resolve/main/Protein_cDNA_Proteolysis_Dataset/Tsuboyama2023_Dataset2_Dataset3_20230416.csv?download=true"
        ],
        "file_name": ["Tsuboyama2023_Dataset2_Dataset3_20230416.csv"],
    },
    "ProteinGym DMS Substitutions Dataset": {
        "paper_title": "ProteinGym: Large-Scale Benchmarks for Protein Design and Fitness Prediction",
        "official_doi": "https://doi.org/10.1101/2023.12.07.570727",
        "huggingface_repos": [
            "datasets/xulab-research/MutCleaner/resolve/main/ProteinGym_DMS_Substitutions_Dataset/DMS_ProteinGym_substitutions.zip?download=true"
        ],
        "file_name": ["ProteinGym_DMS_substitutions.zip"],
    },
    "Human Domainome Dataset": {
        "paper_title": "Site-saturation mutagenesis of 500 human protein domains",
        "official_doi": "https://doi.org/10.1038/s41586-024-08370-4",
        "huggingface_repos": [
            "datasets/xulab-research/MutCleaner/resolve/main/Human_Domainome_Dataset/SupplementaryTable2.txt?download=true",
            "datasets/xulab-research/MutCleaner/resolve/main/Human_Domainome_Dataset/SupplementaryTable4.txt?download=true",
            "datasets/xulab-research/MutCleaner/resolve/main/Human_Domainome_Dataset/wild_type.fasta?download=true",
        ],
        "file_name": [
            "SupplementaryTable2.txt",
            "SupplementaryTable4.txt",
            "wild_type.fasta",
        ],
        "sub_datasets": {
            "Human Domainome Sup2 Dataset": {
                "huggingface_repos": [
                    "datasets/xulab-research/MutCleaner/resolve/main/Human_Domainome_Dataset/SupplementaryTable2.txt?download=true",
                ],
                "file_name": ["SupplementaryTable2.txt"],
            },
            "Human Domainome Sup4 Dataset": {
                "huggingface_repos": [
                    "datasets/xulab-research/MutCleaner/resolve/main/Human_Domainome_Dataset/SupplementaryTable4.txt?download=true",
                    "datasets/xulab-research/MutCleaner/resolve/main/Human_Domainome_Dataset/wild_type.fasta?download=true",
                ],
                "file_name": [
                    "SupplementaryTable4.txt",
                    "wild_type.fasta",
                ],
            },
        },
    },
    "ΔΔG Dataset": {
        "paper_title": "Improving the prediction of protein stability changes upon mutations by geometric learning and a pre-training strategy",
        "official_doi": "https://doi.org/10.1038/s43588-024-00716-2",
        "huggingface_repos": [
            "datasets/xulab-research/MutCleaner/resolve/main/ddG_Dataset/M1261.csv?download=true",
            "datasets/xulab-research/MutCleaner/resolve/main/ddG_Dataset/S461.csv?download=true",
            "datasets/xulab-research/MutCleaner/resolve/main/ddG_Dataset/S669.csv?download=true",
            "datasets/xulab-research/MutCleaner/resolve/main/ddG_Dataset/S783.csv?download=true",
            "datasets/xulab-research/MutCleaner/resolve/main/ddG_Dataset/S8754.csv?download=true",
        ],
        "file_name": [
            "M1261.csv",
            "S461.csv",
            "S669.csv",
            "S783.csv",
            "S8754.csv",
        ],
        "sub_datasets": {
            "M1261": {
                "huggingface_repos": [
                    "datasets/xulab-research/MutCleaner/resolve/main/ddG_Dataset/M1261.csv?download=true"
                ],
                "file_name": ["M1261.csv"],
            },
            "S461": {
                "huggingface_repos": [
                    "datasets/xulab-research/MutCleaner/resolve/main/ddG_Dataset/S461.csv?download=true"
                ],
                "file_name": ["S461.csv"],
            },
            "S669": {
                "huggingface_repos": [
                    "datasets/xulab-research/MutCleaner/resolve/main/ddG_Dataset/S669.csv?download=true"
                ],
                "file_name": ["S669.csv"],
            },
            "S783": {
                "huggingface_repos": [
                    "datasets/xulab-research/MutCleaner/resolve/main/ddG_Dataset/S783.csv?download=true"
                ],
                "file_name": ["S783.csv"],
            },
            "S8754": {
                "huggingface_repos": [
                    "datasets/xulab-research/MutCleaner/resolve/main/ddG_Dataset/S8754.csv?download=true"
                ],
                "file_name": ["S8754.csv"],
            },
        },
    },
    "ΔTm Dataset": {
        "paper_title": "Improving the prediction of protein stability changes upon mutations by geometric learning and a pre-training strategy",
        "official_doi": "https://doi.org/10.1038/s43588-024-00716-2",
        "huggingface_repos": [
            "datasets/xulab-research/MutCleaner/resolve/main/dTm_Dataset/S4346.csv?download=true",
            "datasets/xulab-research/MutCleaner/resolve/main/dTm_Dataset/S571.csv?download=true",
        ],
        "file_name": [
            "S4346.csv",
            "S571.csv",
        ],
        "sub_datasets": {
            "S4346": {
                "huggingface_repos": [
                    "datasets/xulab-research/MutCleaner/resolve/main/dTm_Dataset/S4346.csv?download=true"
                ],
                "file_name": ["S4346.csv"],
            },
            "S557": {
                "huggingface_repos": [
                    "datasets/xulab-research/MutCleaner/resolve/main/dTm_Dataset/S557.csv?download=true"
                ],
                "file_name": ["S557.csv"],
            },
        },
    },
    "ArchStabMS1E10 Epistasis Dataset": {
        "paper_title": "The genetic architecture of protein stability",
        "official_doi": "https://doi.org/10.1038/s41586-024-07966-0",
        "huggingface_repos": [
            "datasets/xulab-research/MutCleaner/resolve/main/ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Sup4_Dataset.csv?download=true",
            "datasets/xulab-research/MutCleaner/resolve/main/ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Sup5_Dataset.csv?download=true",
        ],
        "file_name": ["ArchStabMS1E10_Epistasis_Sup4_Dataset.csv","ArchStabMS1E10_Epistasis_Sup5_Dataset.csv"],
        "sub_datasets": {
            "ArchStabMS1E10_Epistasis_Sup4_Dataset": {
                "huggingface_repos": [
                    "datasets/xulab-research/MutCleaner/resolve/main/ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Sup4_Dataset.csv?download=true"
                ],
                "file_name": ["ArchStabMS1E10_Epistasis_Sup4_Dataset.csv"],
            },
            "ArchStabMS1E10_Epistasis_Sup5_Dataset": {
                "huggingface_repos": [
                    "datasets/xulab-research/MutCleaner/resolve/main/ArchStabMS1E10_Epistasis_Dataset/ArchStabMS1E10_Epistasis_Sup5_Dataset.csv?download=true"
                ],
                "file_name": ["ArchStabMS1E10_Epistasis_Sup5_Dataset.csv"],
            }
        }
    },
    "Human Myoglobin Epistasis Dataset":{
        "paper_title": "Decoding Stability and Epistasis in Human Myoglobin by Deep Mutational Scanning and Codon-level Machine Learning",
        "official_doi": "https://doi.org/10.1101/2024.02.24.581358",
        "huggingface_repos": [
            "datasets/xulab-research/MutCleaner/resolve/main/Human_Myoglobin_Epistasis_Dataset/Human_Myoglobin_Epistasis_Dataset.csv?download=true"
        ],
        "file_name": ["Human_Myoglobin_Epistasis_Dataset.csv"],
    },
    "CTXM Epistasis Dataset":{
        "paper_title": "Network of epistatic interactions in an enzyme active site revealed by DMS",
        "official_doi": "https://doi.org/10.1073/pnas.2313513121",
        "huggingface_repos": [
            "datasets/xulab-research/MutCleaner/resolve/main/CTXM_Epistasis_Dataset/Doubles_A3_processed.txt?download=true",
            "datasets/xulab-research/MutCleaner/resolve/main/CTXM_Epistasis_Dataset/Doubles_C2_processed.txt?download=true",
        ],
        "file_name": [
            "Doubles_A3_processed.txt",
            "Doubles_C2_processed.txt",
        ],
        "sub_datasets": {
            "CTXM_Ampicillin_Epistasis_Dataset": {
                "huggingface_repos": [
                    "datasets/xulab-research/MutCleaner/resolve/main/CTXM_Epistasis_Dataset/Doubles_A3_processed.txt?download=true"
                ],
                "file_name": ["Doubles_A3_processed.txt"],
            },
            "CTXM_Cefotaxime_Epistasis_Dataset": {
                "huggingface_repos": [
                    "datasets/xulab-research/MutCleaner/resolve/main/CTXM_Epistasis_Dataset/Doubles_C2_processed.txt?download=true"
                ],
                "file_name": ["Doubles_C2_processed.txt"],
            },
        }
    },
    "TrpB Epistasis Dataset":{
        "paper_title": "A combinatorially complete epistatic fitness landscape in an enzyme active site",
        "official_doi": "https://doi.org/10.1073/pnas.2400439121",
        "huggingface_repos": [
            "datasets/xulab-research/MutCleaner/resolve/main/TrpB_Epistasis_Dataset/TrpB_Epistasis_Dataset.csv?download=true"
        ],
        "file_name": ["TrpB_Epistasis_Dataset.csv"],
    },
    "Antitoxin ParD3 Epistasis Dataset":{
        "paper_title": "Antitoxin_ParD3_datasets",
        "official_doi": "https://doi.org/10.1038/s41467-024-45621-4",
        "huggingface_repos": [
            "datasets/xulab-research/MutCleaner/resolve/main/Antitoxin_ParD3_Epistasis_Dataset/Antitoxin_ParD3_Epistasis_Dataset.csv?download=true"
        ],
        "file_name": ["Antitoxin_ParD3_Epistasis_Dataset.csv"],
    },
    "RBD Antibody Dataset": {
        "paper_title": "RBD_Antibody_datasets",
        "official_doi": None,
        "huggingface_repos": [
            "datasets/xulab-research/MutCleaner/resolve/main/RBD_Antibody_Dataset/SARS-CoV-2-RBD_MAP_Moderna.csv?download=true",
            "datasets/xulab-research/MutCleaner/resolve/main/RBD_Antibody_Dataset/SARS-CoV-2-RBD_MAP_Rockefeller.csv?download=true",
            "datasets/xulab-research/MutCleaner/resolve/main/RBD_Antibody_Dataset/SARS-CoV-2-RBD_MAP_Vir_mAbs.csv?download=true",
        ],
        "file_name": [
            "SARS-CoV-2-RBD_MAP_Moderna.csv",
            "SARS-CoV-2-RBD_MAP_Rockefeller.csv",
            "SARS-CoV-2-RBD_MAP_Vir_mAbs.csv",
        ],
        "sub_datasets": {
            "Moderna": {
                "huggingface_repos": [
                    "datasets/xulab-research/MutCleaner/resolve/main/RBD_Antibody_Dataset/SARS-CoV-2-RBD_MAP_Moderna.csv?download=true"
                ],
                "file_name": ["SARS-CoV-2-RBD_MAP_Moderna.csv"],
            },
            "Rockefeller": {
                "huggingface_repos": [
                    "datasets/xulab-research/MutCleaner/resolve/main/RBD_Antibody_Dataset/SARS-CoV-2-RBD_MAP_Rockefeller.csv?download=true"
                ],
                "file_name": ["SARS-CoV-2-RBD_MAP_Rockefeller.csv"],
            },
            "Vir_mAbs": {
                "huggingface_repos": [
                    "datasets/xulab-research/MutCleaner/resolve/main/RBD_Antibody_Dataset/SARS-CoV-2-RBD_MAP_Vir_mAbs.csv?download=true"
                ],
                "file_name": ["SARS-CoV-2-RBD_MAP_Vir_mAbs.csv"],
            },
        },
    },
    "RBD ACE2 Dataset": {
        "paper_title": "RBD_ACE2_datasets",
        "official_doi": None,
        "huggingface_repos": [
            "datasets/xulab-research/MutCleaner/resolve/main/RBD_ACE2_Dataset/SARS-CoV-2-RBD_DMS_Omicron-EG5-FLip-BA286_bc_binding.csv?download=true",
            "datasets/xulab-research/MutCleaner/resolve/main/RBD_ACE2_Dataset/SARS-CoV-2-RBD_DMS_Omicron-XBB-BQ_bc_binding.csv?download=true",
            "datasets/xulab-research/MutCleaner/resolve/main/RBD_ACE2_Dataset/SARS-CoV-2-RBD_DMS_Omicron_bc_binding.csv?download=true",
            "datasets/xulab-research/MutCleaner/resolve/main/RBD_ACE2_Dataset/SARS-CoV-2-RBD_DMS_variants_bc_binding.csv?download=true",
            "datasets/xulab-research/MutCleaner/resolve/main/RBD_ACE2_Dataset/SARS-CoV-2-RBD_Delta_bc_binding.csv?download=true",
        ],
        "file_name": [
            "SARS-CoV-2-RBD_DMS_Omicron-EG5-FLip-BA286_bc_binding.csv",
            "SARS-CoV-2-RBD_DMS_Omicron-XBB-BQ_bc_binding.csv",
            "SARS-CoV-2-RBD_DMS_Omicron_bc_binding.csv",
            "SARS-CoV-2-RBD_DMS_variants_bc_binding.csv",
            "SARS-CoV-2-RBD_Delta_bc_binding.csv",
        ],
        "sub_datasets": {
            "Omicron_EG5_FLip_BA286": {
                "huggingface_repos": [
                    "datasets/xulab-research/MutCleaner/resolve/main/RBD_ACE2_Dataset/SARS-CoV-2-RBD_DMS_Omicron-EG5-FLip-BA286_bc_binding.csv?download=true"
                ],
                "file_name": ["SARS-CoV-2-RBD_DMS_Omicron-EG5-FLip-BA286_bc_binding.csv"],
            },
            "Omicron_XBB_BQ": {
                "huggingface_repos": [
                    "datasets/xulab-research/MutCleaner/resolve/main/RBD_ACE2_Dataset/SARS-CoV-2-RBD_DMS_Omicron-XBB-BQ_bc_binding.csv?download=true"
                ],
                "file_name": ["SARS-CoV-2-RBD_DMS_Omicron-XBB-BQ_bc_binding.csv"],
            },
            "Omicron": {
                "huggingface_repos": [
                    "datasets/xulab-research/MutCleaner/resolve/main/RBD_ACE2_Dataset/SARS-CoV-2-RBD_DMS_Omicron_bc_binding.csv?download=true"
                ],
                "file_name": ["SARS-CoV-2-RBD_DMS_Omicron_bc_binding.csv"],
            },
            "DMS_variants": {
                "huggingface_repos": [
                    "datasets/xulab-research/MutCleaner/resolve/main/RBD_ACE2_Dataset/SARS-CoV-2-RBD_DMS_variants_bc_binding.csv?download=true"
                ],
                "file_name": ["SARS-CoV-2-RBD_DMS_variants_bc_binding.csv"],
            },
            "Delta": {
                "huggingface_repos": [
                    "datasets/xulab-research/MutCleaner/resolve/main/RBD_ACE2_Dataset/SARS-CoV-2-RBD_Delta_bc_binding.csv?download=true"
                ],
                "file_name": ["SARS-CoV-2-RBD_Delta_bc_binding.csv"],
            },
        },
    },
    "Chitosanase dTm Dataset": {
        "paper_title": "Chitosanase_dTm_dataset",
        "official_doi": None,
        "huggingface_repos": [
            "datasets/xulab-research/MutCleaner/resolve/main/Chitosanase_dTm_Dataset/Chitosanase_dTm_Dataset.csv?download=true",
        ],
        "file_name": ["Chitosanase_dTm_Dataset.csv"],
    },
    "MGnify ddG Dataset": {
        "paper_title": "MGnify_ddG_dataset",
        "official_doi": None,
        "huggingface_repos": [
            "datasets/xulab-research/MutCleaner/resolve/main/MGnify_ddG_Dataset/MGnify_ddG_Dataset.csv?download=true",
        ],
        "file_name": ["MGnify_ddG_Dataset.csv"],
    },
    "Codon cDNA Proteolysis Dataset": {
        "paper_title": "Mega-scale experimental analysis of protein folding stability in biology and design",
        "official_doi": "https://doi.org/10.1038/s41586-023-06328-6",
        "huggingface_repos": [
            "datasets/xulab-research/MutCleaner/resolve/main/Codon_cDNA_Proteolysis_Dataset/Codon_cDNA_Proteolysis_Dataset.csv?download=true",
            "datasets/xulab-research/MutCleaner/resolve/main/Codon_cDNA_Proteolysis_Dataset/wt.fasta?download=true"            
        ],
        "file_name": ["Codon_cDNA_Proteolysis_Dataset.csv","wt.fasta"],
    },
}


def list_datasets_with_built_in_cleaners() -> None:
    """
    List built-in datasets with predefined processing pipelines.

    These are public datasets for which this package includes pre-defined
    data cleaning pipelines. The datasets themselves are not distributed
    with the package and must be downloaded manually.

    You can also define custom cleaner functions for your own datasets using
    the same `@pipeline_step` framework.

    Predefined datasets:

    - cDNA Proteolysis Dataset
    - Codon cDNA Proteolysis Dataset
    - ProteinGym DMS Substitutions Dataset
    - Human Domainome Dataset
    - ΔΔG Dataset
    - ΔTm Dataset
    - Antitoxin ParD3 Epistasis Dataset
    - TrpB Epistasis Dataset
    - Human Myoglobin Epistasis Dataset
    - CTXM Epistasis Dataset
    - RBD Antibody Dataset
    - RBD ACE2 Dataset
    - Chitosanase dTm Dataset
    - MGnify ddG Dataset
    """
    print("Public datasets with ready-to-use cleaning pipelines:")
    for key, info in DATASETS.items():
        print(f"- {key}: {info['paper_title']}")
        print(f"  - Official DOI: {info['official_doi']}")
    print(
        "\nUse the `show_download_instructions` function to see detailed download instructions."
    )


def show_download_instructions(dataset_key: str) -> None:
    """
    Show download instructions for a specific dataset.
    """
    info = DATASETS.get(dataset_key)
    if not info:
        raise KeyError(f"Dataset key not found: {dataset_key}")

    print(f"Dataset: {info['paper_title']}")
    for i, file in enumerate(info["file_name"]):
        print(f"  - File: {file}")
        print(f"    - Download link: {info['huggingface_repos'][i]}")
    print(f"\nSub-datasets:")
    for sub_dataset, sub_info in info.get("sub_datasets", {}).items():
        print(f"- Sub-dataset: {sub_dataset}")
        for i, file in enumerate(sub_info["file_name"]):
            print(f"  - File: {file}")
            print(f"    - Download link: {sub_info['huggingface_repos'][i]}")
