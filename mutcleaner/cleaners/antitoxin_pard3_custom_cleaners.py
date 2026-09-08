import re
import pandas as pd
import numpy as np
from typing import List
from tqdm import tqdm

from ..core.pipeline import pipeline_step
from ..core.sequence import ProteinSequence

__all__ = ["infer_wt_sequence", "add_wild_type_sequence", "simplify_mutations"]


@pipeline_step
def infer_wt_sequence(
    dataset: pd.DataFrame,
    mutation_col: str,
    sequence_col: str = "mut_seq",
    wt_column: str = "wt_seq",
) -> pd.DataFrame:
    """
    infer wild type seuquence from mutation_sequence and mutation information

    Argument
    --------
        dataset : pd.DataFrame
            the dataset which can be inferred wild-type sequence with mutation sequence and mutation information
        mutation_col : str
            the column contatins mutation information
        sequence_col : str, default="mut_seq"
            the column contains mutation sequence
        wt_column : str, default="wt_seq"
            the column that can be inferred by mutation sequence and mutation information

    Returns
    -------
        dataset: pd.DataFrame
            the dataset contatins added wild-type sequence
    """
    if dataset.empty:
        raise ValueError("Dataset is empty. Cannot infer WT sequence.")

    dataset = dataset.copy()

    seq_len = len(dataset[sequence_col].iat[0])
    sequences: List[str] = dataset[sequence_col].tolist()
    mutation_seqs: List[str] = dataset[mutation_col].tolist()

    wt_sequence_list = []

    try:
        mut_matrix = np.array(mutation_seqs, dtype="U").view("U1").reshape(len(dataset), -1)
        seq_matrix = np.array(sequences, dtype="U").view("U1").reshape(len(dataset), -1)
    except ValueError:
        raise ValueError("Sequences in the DataFrame are not of equal length.")

    for i in range(seq_len):
        match_indices = np.where(mut_matrix[:, i] == "_")[0]

        if match_indices.size > 0:
            idx = match_indices[0]
            original_char = seq_matrix[idx, i]
            wt_sequence_list.append(original_char)
        else:
            print(f"Warning: Position {i} is mutated in All sequences. Cannot recover wild-type sequence.")
            wt_sequence_list.append("X")

    wt_seq = "".join(wt_sequence_list)

    dataset[wt_column] = wt_seq

    return dataset


@pipeline_step
def simplify_mutations(
    dataset: pd.DataFrame,
    mutation_column: str = "mut_info",
    mutation_sep: str = ":",
) -> pd.DataFrame:
    """
    Simplify mutation strings by removing no-change tokens like L47L.
    Keeps unknown tokens unchanged (does not filter rows).

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    >>> 'name':["protein1", "protein1", "protein1"]
    >>> 'mut_info':["L47L:D51D:I52I","L47A:D51D:I52A","L47_:D51D:I52I"]
    >>> })
    >>> df = simplify_mutations(df)
    >>> df[mut_info]
        pd.Series("WT","L47A,I52A","L47_") # only simplify the mutations but not fliter
    """
    if mutation_column not in dataset.columns:
        raise ValueError(f"the {mutation_column} is not in the dataset")

    dataset = dataset.copy()

    # only for simplifying
    # validate the legality of mutations
    mutation_regex = re.compile(r"^([A-Z\*_])(\d+)([A-Z\*_])$")

    def simplify_single_mutations(mutations):
        if pd.isna(mutations):
            return mutations

        mutations = str(mutations).strip()
        if not mutations:
            return mutations

        keep = []
        for single_mutation_tok in mutations.split(mutation_sep):
            single_mutation_tok = single_mutation_tok.strip()
            if not single_mutation_tok:
                continue
            single_mutation = mutation_regex.match(single_mutation_tok)
            if single_mutation is None:
                keep.append(single_mutation_tok)
                continue
            wt, pos, mut = (
                single_mutation.group(1),
                single_mutation.group(2),
                single_mutation.group(3),
            )
            if wt != mut:
                keep.append(f"{wt}{pos}{mut}")

        return "WT" if len(keep) == 0 else ",".join(keep)

    dataset[mutation_column] = dataset[mutation_column].apply(simplify_single_mutations)

    return dataset
