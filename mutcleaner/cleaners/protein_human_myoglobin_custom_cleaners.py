from __future__ import annotations

import re
from typing import TYPE_CHECKING

import pandas as pd
from tqdm import tqdm

from ..core.codon import CodonTable
from ..core.pipeline import multiout_step

if TYPE_CHECKING:
    from typing import Literal, Tuple

__all__ = ["convert_codon_to_amino_acid"]


@multiout_step(main="success", failed="failed")
def convert_codon_to_amino_acid(
    dataset: pd.DataFrame,
    codon_column: str = "codon_mutations",
    amino_acid_column: str = "mut_info",
    seq_type: Literal["DNA", "RNA"] = "DNA",
    strict: bool = True,
    drop_codon_column: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Convert codon mutation annotations to amino acid mutations.

    Codon substitutions are translated using the standard genetic code while
    preserving the original mutation positions. Records containing mutations
    that introduce stop codons are separated from successfully converted
    records.

    Parameters
    ----------
    dataset : pandas.DataFrame
        Input dataset containing codon mutation annotations.
    codon_column : str, default="codon_mutations"
        Column containing codon mutation annotations. Multiple mutations within
        one record must be separated by commas, for example
        ``"ATG0GTG,GAA5GAC"``.
    amino_acid_column : str, default="aa"
        Name of the output column containing converted amino acid mutation
        annotations.
    seq_type : {"DNA", "RNA"}, default="DNA"
        Type of nucleotide sequence used to interpret codons.
    strict : bool, default=True
        Whether invalid codon mutation annotations should raise a
        ``ValueError``. If False, invalid records are returned in the failed
        DataFrame with an appropriate ``error_message``.
    drop_codon_column : bool, default=False
        Whether to remove ``codon_column`` from successfully converted records.
        The original codon column is retained in failed records for
        traceability.

    Returns
    -------
    successful : pandas.DataFrame
        Records successfully converted to amino acid mutation annotations and
        not introducing stop codons.
    failed : pandas.DataFrame
        Records excluded during conversion. With ``strict=True``, this normally
        contains records introducing stop codons. With ``strict=False``,
        malformed or unrecognized codon annotations are also included.
        An ``error_message`` column describes the reason for exclusion.

    Raises
    ------
    ValueError
        If ``codon_column`` is absent, ``seq_type`` is unsupported, or an
        invalid codon mutation annotation is encountered while ``strict=True``.

    Notes
    -----
    Mutation positions are copied directly from the codon mutation annotation.
    This function does not convert between zero-based and one-based indexing.

    Synonymous codon substitutions are retained during conversion. For example,
    ``AAA0AAG`` is converted to ``K0K``. Such redundant amino acid mutations
    may be removed by subsequent mutation validation steps.

    Examples
    --------
    >>> df = pd.DataFrame({
    ...     "codon_mutations": ["ATG0GTG", "TGG1TAA"],
    ...     "label": [1.2, -0.5],
    ... })
    >>> successful, failed = convert_codon_to_amino_acid(
    ...     df,
    ...     codon_column="codon_mutations",
    ...     amino_acid_column="mut_info",
    ... )
    >>> successful["mut_info"].tolist()
    ['M0V']
    >>> len(failed)
    1
    """
    if codon_column not in dataset.columns:
        raise ValueError(f"Column '{codon_column}' not found in dataset")

    seq_type = seq_type.upper()
    if seq_type not in {"DNA", "RNA"}:
        raise ValueError(f"seq_type must be 'DNA' or 'RNA', got '{seq_type}'")

    tqdm.write("Converting codon mutations to amino acid mutations...")

    alphabet_pattern = r"[ACGT]" if seq_type == "DNA" else r"[ACGU]"
    token_re = re.compile(
        rf"^({alphabet_pattern}{{3}})(\d+)({alphabet_pattern}{{3}})$",
        re.IGNORECASE,
    )
    table = CodonTable.get_standard_table(seq_type=seq_type)

    def _convert_field(value) -> tuple[object, str | None]:
        if pd.isna(value) or not str(value).strip():
            message = "Missing codon mutation information"
            if strict:
                raise ValueError(message)
            return pd.NA, message

        tokens = [token.strip() for token in str(value).split(",") if token.strip()]
        amino_acid_mutations = []

        for token in tokens:
            match = token_re.fullmatch(token)
            if match is None:
                message = f"Invalid codon mutation format: {token}"
                if strict:
                    raise ValueError(message)
                return pd.NA, message

            wt_codon, position, mut_codon = match.groups()
            wt_codon = wt_codon.upper()
            mut_codon = mut_codon.upper()

            wt_aa = table.translate_codon(wt_codon)
            mut_aa = table.translate_codon(mut_codon)

            if "X" in (wt_aa, mut_aa):
                message = (
                    f"Unknown codon translation: {token} "
                    f"(wild_type={wt_aa}, mutant={mut_aa})"
                )
                if strict:
                    raise ValueError(message)
                return pd.NA, message

            if table.is_stop_codon(mut_codon):
                return pd.NA, f"Mutation introduces a stop codon: {token}"

            amino_acid_mutations.append(f"{wt_aa}{int(position)}{mut_aa}")

        return ",".join(amino_acid_mutations), None

    converted = dataset[codon_column].apply(_convert_field)
    failed_mask = converted.map(lambda result: result[1] is not None)

    successful = dataset.loc[~failed_mask].copy()
    failed = dataset.loc[failed_mask].copy()

    successful[amino_acid_column] = converted.loc[~failed_mask].map(
        lambda result: result[0]
    )

    if not failed.empty:
        failed["error_message"] = converted.loc[failed_mask].map(
            lambda result: result[1]
        )

    if drop_codon_column:
        successful = successful.drop(columns=[codon_column])

    tqdm.write(
        f"Codon-to-amino-acid conversion completed: "
        f"{len(successful)} records retained, {len(failed)} records filtered."
    )

    return successful, failed