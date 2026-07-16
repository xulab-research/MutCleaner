# mutcleaner/__init__.py
"""
MutCleaner: An efficient framework for cleaning, standardizing, and processing biological mutation data.
"""

__author__ = "Ziyu Shi and Yuxiang Tang"

__version__ = "0.8.0"

from .core import (
    # Alphabet
    alphabet,
    # Codon
    codon,
    # Mutation
    mutation,
    # Sequence
    sequence,
    # Dataset
    MutationDataset,
    # Pipeline
    Pipeline,
    pipeline_step,
    multiout_step,
    create_pipeline,
)

from .cleaners import (
    basic_cleaners,
    rbd_antibody_cleaner,
    cdna_proteolysis_cleaner,
    human_domainome_sup2_cleaner,
    human_domainome_sup4_cleaner,
    proteingym_dms_substitutions_cleaner,
    ddg_dtm_cleaners,
    archstabms_1e10_sup4_cleaner,
    archstabms_1e10_sup5_cleaner,
    human_myoglobin_cleaner,
    ctxm_cleaner,
    trpb_cleaner,
    antitoxin_pard3_cleaner,
    rbd_ace2_cleaner,
    chitosanase_dtm_cleaner,
    mgnify_ddg_cleaner,
)

from .utils import (
    download,
    download_cdna_proteolysis_source_file,
    download_proteingym_source_file,
    download_human_domainome_source_file,
    download_ddg_dtm_source_file,
    list_datasets_with_built_in_cleaners,
    show_download_instructions,
    download_archstabms1e10_source_file,
    download_human_myoglobin_source_file,
    download_ctxm_source_file,
    download_trpb_source_file,
    download_antitoxin_pard3_source_file,
    download_rbd_antibody_source_file,
    download_rbd_ace2_source_file,
    download_chitosanase_dtm_source_file,
    download_mgnify_ddg_source_file,
)

__all__ = [
    "alphabet",
    "codon",
    "mutation",
    "sequence",
    "MutationDataset",
    "Pipeline",
    "pipeline_step",
    "multiout_step",
    "create_pipeline",
    "basic_cleaners",
    "rbd_antibody_cleaner",
    "cdna_proteolysis_cleaner",
    "human_domainome_sup2_cleaner",
    "human_domainome_sup4_cleaner",
    "proteingym_dms_substitutions_cleaner",
    "ddg_dtm_cleaners",
    "archstabms_1e10_sup4_cleaner",
    "archstabms_1e10_sup5_cleaner",
    "human_myoglobin_cleaner",
    "ctxm_cleaner",
    "trpb_cleaner",
    "antitoxin_pard3_cleaner",
    "rbd_ace2_cleaner",
    "chitosanase_dtm_cleaner",
    "mgnify_ddg_cleaner",
    "download",
    "download_cdna_proteolysis_source_file",
    "download_proteingym_source_file",
    "download_human_domainome_source_file",
    "download_ddg_dtm_source_file",
    "list_datasets_with_built_in_cleaners",
    "show_download_instructions",
    "download_archstabms1e10_source_file",
    "download_human_myoglobin_source_file",
    "download_ctxm_source_file",
    "download_trpb_source_file",
    "download_antitoxin_pard3_source_file",
    "download_rbd_antibody_source_file",
    "download_rbd_ace2_source_file",
    "download_chitosanase_dtm_source_file",
    "download_mgnify_ddg_source_file",
]