# mutcleaner/utils/sequence_io.py
"""Utilities for reading and writing sequence files without BioPython dependency."""

from __future__ import annotations

import json
import pandas as pd
import warnings
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any, Callable, Dict, Optional, Tuple, Union


def parse_uniprot_header(header: str) -> Tuple[str, Dict[str, str]]:
    """Parse UniProt FASTA header to extract ID and metadata

    Parameters
    ----------
    header : str
        FASTA header line (without '>')

    Returns
    -------
    Tuple[str, Dict[str, str]]
        (sequence_id, metadata_dict)

    Examples
    --------
    >>> parse_uniprot_header("sp|P12345|PROT_HUMAN Protein description OS=Homo sapiens")
    ('P12345', {'db': 'sp', 'entry_name': 'PROT_HUMAN', 'description': 'Protein description OS=Homo sapiens'})
    >>> parse_uniprot_header("P12345|PROT_HUMAN Description")
    ('P12345', {'entry_name': 'PROT_HUMAN', 'description': 'Description'})
    >>> parse_uniprot_header("P12345")
    ('P12345', {})
    """
    metadata = {}

    if "|" in header:
        parts = header.split("|")
        if len(parts) >= 3 and parts[0] in ["sp", "tr"]:
            # Format: sp|P12345|PROT_HUMAN Description
            metadata["db"] = parts[0]
            seq_id = parts[1]
            remaining = "|".join(parts[2:])
            if " " in remaining:
                entry_name, description = remaining.split(" ", 1)
                metadata["entry_name"] = entry_name
                metadata["description"] = description
            else:
                metadata["entry_name"] = remaining
        elif len(parts) >= 2:
            # Format: P12345|PROT_HUMAN Description
            seq_id = parts[0]
            remaining = "|".join(parts[1:])
            if " " in remaining:
                entry_name, description = remaining.split(" ", 1)
                metadata["entry_name"] = entry_name
                metadata["description"] = description
            else:
                metadata["entry_name"] = remaining
        else:
            seq_id = parts[0]
    else:
        # No pipe, take everything before first space
        if " " in header:
            seq_id, description = header.split(" ", 1)
            metadata["description"] = description
        else:
            seq_id = header

    return seq_id, metadata


def parse_ncbi_header(header: str) -> Tuple[str, Dict[str, str]]:
    """Parse NCBI FASTA header to extract ID and metadata

    Parameters
    ----------
    header : str
        FASTA header line (without '>')

    Returns
    -------
    Tuple[str, Dict[str, str]]
        (sequence_id, metadata_dict)

    Examples
    --------
    >>> parse_ncbi_header("gi|123456|ref|NP_000001.1| protein description [Homo sapiens]")
    ('NP_000001.1', {'gi': '123456', 'db': 'ref', 'description': 'protein description [Homo sapiens]'})
    >>> parse_ncbi_header("NP_000001.1 protein description")
    ('NP_000001.1', {'description': 'protein description'})
    """
    metadata = {}

    if header.startswith("gi|"):
        # NCBI format: gi|123456|db|accession| description
        parts = header.split("|")
        if len(parts) >= 4:
            metadata["gi"] = parts[1]
            metadata["db"] = parts[2]
            seq_id = parts[3]
            if len(parts) > 4:
                description = "|".join(parts[4:]).strip()
                if description:
                    metadata["description"] = description
        else:
            seq_id = header
    else:
        # Simple format: accession description
        if " " in header:
            seq_id, description = header.split(" ", 1)
            metadata["description"] = description
        else:
            seq_id = header

    return seq_id, metadata


def parse_simple_header(header: str) -> Tuple[str, Dict[str, str]]:
    """Simple header parser that uses the first word as ID

    Parameters
    ----------
    header : str
        FASTA header line (without '>')

    Returns
    -------
    Tuple[str, Dict[str, str]]
        (sequence_id, metadata_dict)

    Examples
    --------
    >>> parse_simple_header("GENE1 some description text")
    ('GENE1', {'description': 'some description text'})
    >>> parse_simple_header("GENE1")
    ('GENE1', {})
    """
    if " " in header:
        seq_id, description = header.split(" ", 1)
        return seq_id, {"description": description}
    else:
        return header, {}


def parse_custom_delimiter_header(delimiter: str = "|", id_position: int = 0) -> Callable:
    """Create a header parser for custom delimiter-based formats

    Parameters
    ----------
    delimiter : str, default='|'
        Delimiter character to split the header
    id_position : int, default=0
        Position of the ID in the split parts (0-based)

    Returns
    -------
    Callable
        Header parser function

    Examples
    --------
    >>> parser = parse_custom_delimiter_header('|', 1)
    >>> parser("db|GENE1|other|info")
    ('GENE1', {'parts': ['db', 'other', 'info']})
    """

    def parser(header: str) -> Tuple[str, Dict[str, str]]:
        parts = header.split(delimiter)
        if len(parts) > id_position:
            seq_id = parts[id_position]
            # Store other parts in metadata
            other_parts_idx = [i for i in range(len(parts)) if i != id_position]
            return seq_id, dict(zip([f"parts{i}" for i in other_parts_idx], parts[other_parts_idx:]))
        else:
            return header, {}

    return parser


def parse_fasta(
    file_path: Union[str, Path],
    header_parser: Optional[Callable[[str], Tuple[str, Dict[str, str]]]] = None,
    clean_sequence: bool = True,
) -> Dict[str, Dict[str, Any]]:
    """Parse FASTA file with custom header parsing

    Parameters
    ----------
    file_path : Union[str, Path]
        Path to FASTA file
    header_parser : Optional[Callable], default=None
        Function to parse headers. Should take header string and return (id, metadata).
        If None, uses parse_uniprot_header as default.
    clean_sequence : bool, default=True
        Whether to clean sequences (remove whitespace, numbers, etc.)

    Returns
    -------
    Dict[str, Dict[str, Any]]
        Dictionary mapping sequence IDs to {'sequence': str, 'metadata': dict}

    Examples
    --------
    >>> # Use default UniProt parser
    >>> sequences = parse_fasta("proteins.fasta")

    >>> # Use NCBI parser
    >>> sequences = parse_fasta("ncbi_proteins.fasta", header_parser=parse_ncbi_header)

    >>> # Use simple parser
    >>> sequences = parse_fasta("genes.fasta", header_parser=parse_simple_header)

    >>> # Custom parser
    >>> def my_parser(header):
    ...     return header.split('_')[0], {'full_header': header}
    >>> sequences = parse_fasta("custom.fasta", header_parser=my_parser)
    """
    if header_parser is None:
        header_parser = parse_uniprot_header

    sequences = {}
    current_id = None
    current_seq = []
    current_metadata = {}

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                # Save previous sequence if exists
                if current_id is not None:
                    seq = "".join(current_seq)
                    sequences[current_id] = {
                        "sequence": seq,
                        "metadata": current_metadata,
                    }

                # Parse new header
                header = line[1:].strip()
                try:
                    current_id, current_metadata = header_parser(header)
                except Exception as e:
                    warnings.warn(f"Failed to parse header '{header}': {e}. Using full header as ID.")
                    current_id = header
                    current_metadata = {}

                current_seq = []
            else:
                # Sequence line
                if clean_sequence:
                    # Remove any whitespace and numbers
                    cleaned_line = "".join(c for c in line if c.isalpha())
                else:
                    cleaned_line = line.strip()

                if cleaned_line:
                    current_seq.append(cleaned_line)

    # Save the final parsed sequence.
    if current_id is not None:
        seq = "".join(current_seq)
        sequences[current_id] = {"sequence": seq, "metadata": current_metadata}

    return sequences


def load_sequences(
    file_path: Union[str, Path],
    header_parser: Optional[Callable[[str], Tuple[str, Dict[str, str]]]] = None,
    format: Optional[str] = None,
    id_column: Optional[str] = None,
    sequence_column: Optional[str] = None,
) -> Dict[str, str]:
    """Load sequences from various file formats

    Parameters
    ----------
    file_path : Union[str, Path]
        Path to sequence file
    header_parser : Optional[Callable], default=None
        Function to parse FASTA headers (only used for FASTA format)
    format : Optional[str], default=None
        File format. If None, inferred from extension.
        Supported: 'fasta', 'csv', 'tsv', 'json'
    id_column : Optional[str], default=None
        Column name for sequence IDs (CSV/TSV only)
    sequence_column : Optional[str], default=None
        Column name for sequences (CSV/TSV only)

    Returns
    -------
    Dict[str, str]
        Dictionary mapping sequence IDs to sequences

    Examples
    --------
    >>> # Load UniProt FASTA
    >>> seqs = load_sequences("uniprot.fasta")

    >>> # Load FASTA with custom parser
    >>> seqs = load_sequences("genes.fasta", header_parser=parse_simple_header)

    >>> # Load CSV with specified columns
    >>> seqs = load_sequences("sequences.csv", id_column="protein_id", sequence_column="aa_sequence")

    >>> # Load with automatic column detection
    >>> seqs = load_sequences("sequences.csv")
    """
    path = Path(file_path)

    # Infer format from extension if not specified
    if format is None:
        format = path.suffix.lower().lstrip(".")

    # Normalize format names
    format_map = {
        "fa": "fasta",
        "faa": "fasta",
        "fas": "fasta",
        "txt": "fasta",
    }
    format = format_map.get(format, format)

    if format == "fasta":
        # Load FASTA file
        fasta_data = parse_fasta(file_path, header_parser=header_parser)
        seq_dict = {}
        for seq_id, data in fasta_data.items():
            if isinstance(data, dict) and "sequence" in data:
                seq_dict[seq_id] = data["sequence"]
        return seq_dict

    elif format in ["csv", "tsv"]:
        # Load CSV/TSV
        sep = "\t" if format == "tsv" else ","
        df = pd.read_csv(path, sep=sep)

        # Auto-detect columns if not specified
        if id_column is None or sequence_column is None:
            id_col, seq_col = _detect_sequence_columns(df.columns)
            if id_column is None:
                id_column = id_col
            if sequence_column is None:
                sequence_column = seq_col

        if id_column is None or sequence_column is None:
            raise ValueError(f"Could not detect ID and sequence columns. " f"Please specify id_column and sequence_column. " f"Available columns: {list(df.columns)}")

        # Create dictionary
        seq_dict = {}
        for _, row in df.iterrows():
            id_val = row[id_column]
            seq_val = row[sequence_column]

            if pd.notna(id_val) and pd.notna(seq_val):
                seq_dict[str(id_val).strip()] = str(seq_val).strip()

        return seq_dict

    elif format == "json":
        # Load JSON
        with open(path, "r") as f:
            data = json.load(f)
            # Ensure all values are strings
            return {str(k): str(v) for k, v in data.items()}

    else:
        raise ValueError(f"Unsupported format: {format}. " f"Supported formats: fasta, csv, tsv, json")


def _detect_sequence_columns(columns: pd.Index) -> Tuple[Optional[str], Optional[str]]:
    """Auto-detect ID and sequence columns from column names

    Returns
    -------
    Tuple[Optional[str], Optional[str]]
        (id_column, sequence_column)
    """
    id_col = None
    seq_col = None

    # Patterns for detecting columns
    id_patterns = [
        "uniprot",
        "accession",
        "protein_id",
        "entry",
        "gene_id",
        "gene_name",
        "id",
        "name",
        "identifier",
    ]
    seq_patterns = [
        "sequence",
        "seq",
        "aa_seq",
        "aa_sequence",
        "protein_seq",
        "protein_sequence",
        "peptide",
    ]

    for col in columns:
        col_lower = col.lower()

        # Check sequence patterns
        if seq_col is None:
            for pattern in seq_patterns:
                if pattern in col_lower:
                    seq_col = col
                    break

        # Check ID patterns
        if id_col is None:
            for pattern in id_patterns:
                if pattern in col_lower:
                    id_col = col
                    break

    # Fallback: if only 2 columns, assume first is ID, second is sequence
    if (id_col is None or seq_col is None) and len(columns) == 2:
        warnings.warn(f"Could not identify columns by name. " f"Assuming {columns[0]} is ID and {columns[1]} is sequence.")
        id_col = columns[0]
        seq_col = columns[1]

    return id_col, seq_col


def write_fasta(
    sequences: Union[Dict[str, str], Dict[str, Dict[str, Any]]],
    file_path: Union[str, Path],
    wrap_length: int = 60,
    header_formatter: Optional[Callable[[str, Dict], str]] = None,
) -> None:
    """Write sequences to FASTA file

    Parameters
    ----------
    sequences : Union[Dict[str, str], Dict[str, Dict[str, Any]]]
        Dictionary mapping IDs to sequences or {'sequence': str, 'metadata': dict}
    file_path : Union[str, Path]
        Output file path
    wrap_length : int, default=60
        Line length for sequence wrapping (0 for no wrapping)
    header_formatter : Optional[Callable], default=None
        Function to format headers. Takes (id, metadata) and returns header string.

    Examples
    --------
    >>> # Simple sequences
    >>> seqs = {'GENE1': 'ACDEF', 'GENE2': 'KLMNO'}
    >>> write_fasta(seqs, 'output.fasta')

    >>> # With metadata
    >>> seqs = {
    ...     'P12345': {
    ...         'sequence': 'ACDEF',
    ...         'metadata': {'description': 'Protein 1', 'organism': 'Human'}
    ...     }
    ... }
    >>> write_fasta(seqs, 'output.fasta')
    """

    def default_formatter(seq_id: str, metadata: Dict) -> str:
        if not metadata:
            return seq_id
        if "description" in metadata:
            return f"{seq_id} {metadata['description']}"
        return seq_id

    if header_formatter is None:
        header_formatter = default_formatter

    with open(file_path, "w") as f:
        for seq_id, data in sequences.items():
            # Handle both simple sequences and dict format
            if isinstance(data, str):
                sequence = data
                metadata = {}
            else:
                sequence = data["sequence"]
                metadata = data.get("metadata", {})

            # Write header
            header = header_formatter(seq_id, metadata)
            f.write(f">{header}\n")

            # Write sequence
            if wrap_length > 0:
                for i in range(0, len(sequence), wrap_length):
                    f.write(sequence[i : i + wrap_length] + "\n")
            else:
                f.write(sequence + "\n")
