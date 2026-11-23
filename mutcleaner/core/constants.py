# mutcleaner/core/constants.py
"""Alphabet and genetic-code constants used across mutcleaner.

This module collects IUPAC alphabets for DNA/RNA/amino acids, base-complement
mappings, 1↔3 letter amino-acid code conversions, and the standard genetic code
(codon → amino acid) for DNA and RNA.

Examples
--------
Validate a DNA sequence contains only standard bases::

    all(b in STANDARD_DNA_BASES for b in "ATGCGT")

Translate a DNA codon::

    AA = STANDARD_GENETIC_CODE_DNA["ATG"]   # 'M'

Get the 1-letter code from a PDB-style residue name::

    AA1 = AA3_TO_1["ASP"]    # 'D'
    AA3 = AA1_TO_3["D"]      # 'Asp'

Attributes
----------
STANDARD_DNA_BASES : set of str
    Canonical DNA bases (``{'A','T','C','G'}``).
AMBIGUOUS_DNA_BASES : set of str
    IUPAC ambiguous DNA symbols (e.g. ``'R','Y','S','W','K','M','B','D','H','V','N'``).
STANDARD_RNA_BASES : set of str
    Canonical RNA bases (``{'A','U','C','G'}``).
AMBIGUOUS_RNA_BASES : set of str
    IUPAC ambiguous RNA symbols.
STANDARD_AMINO_ACIDS : set of str
    20 standard amino-acid one-letter codes.
AMBIGUOUS_AMINO_ACIDS : set of str
    Symbols for Asp or Asn (B), Glu or Gln (Z), Leu or Ile (J), and an unknown residue (X).
SEC_PYL_AMINO_ACIDS : set of str
    One-letter codes for selenocysteine (U) and pyrrolysine (O).

DNA_BASE_COMPLEMENTS : dict[str, str]
    DNA Watson-Crick complements.
RNA_BASE_COMPLEMENTS : dict[str, str]
    RNA Watson-Crick complements.

AA3_TO_1 : dict[str, str]
    Amino-acid 3-letter → 1-letter code map (case-tolerant; includes ``'TER' → '*'``).
AA1_TO_3 : dict[str, str]
    Amino-acid 1-letter → 3-letter code map (includes stop ``'*' → 'TER'``).

STANDARD_GENETIC_CODE_DNA : dict[str, str]
    DNA codon table (triplet of ``ATCG`` → ``ACDEFGHIKLMNPQRSTVWY*``).
STANDARD_START_CODONS_DNA : set[str]
    DNA start codons (default ``{'ATG'}``).

STANDARD_GENETIC_CODE_RNA : dict[str, str]
    RNA codon table (triplet of ``AUCG`` → ``ACDEFGHIKLMNPQRSTVWY*``).
STANDARD_START_CODONS_RNA : set[str]
    RNA start codon set (default ``{'AUG'}``).
"""

# fmt: off
# ==== Alphabet Constants ====
STANDARD_DNA_BASES = {"A", "T", "C", "G"}
AMBIGUOUS_DNA_BASES = {"R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"}
STANDARD_RNA_BASES = {"A", "U", "C", "G"}
AMBIGUOUS_RNA_BASES = {"R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"}
STANDARD_AMINO_ACIDS = {
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
}
AMBIGUOUS_AMINO_ACIDS = {"B", "Z", "J", "X"}
SEC_PYL_AMINO_ACIDS = {"U", "O"}

# ==== Base Conversion ====
DNA_BASE_COMPLEMENTS = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
RNA_BASE_COMPLEMENTS = {"A": "U", "U": "A", "C": "G", "G": "C", "N": "N"}

# ==== Amino Acid Conversion ====
AA3_TO_1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", 
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    "TER": "*",
    "Ala": "A", "Cys": "C", "Asp": "D", "Glu": "E", "Phe": "F",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", 
    "Met": "M", "Asn": "N", "Pro": "P", "Gln": "Q", "Arg": "R",
    "Ser": "S", "Thr": "T", "Val": "V", "Trp": "W", "Tyr": "Y",
    "Ter": "*"
}
AA1_TO_3 = {
    "A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe",
    "G": "Gly", "H": "His", "I": "Ile", "K": "Lys", "L": "Leu",
    "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln", "R": "Arg",
    "S": "Ser", "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr",
    "*": "Ter"
}

# ==== Codon AA Conversion ====
STANDARD_GENETIC_CODE_DNA = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
STANDARD_START_CODONS_DNA = {"ATG"}
STANDARD_GENETIC_CODE_RNA = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
STANDARD_START_CODONS_RNA = {"AUG"}
