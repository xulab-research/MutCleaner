import re
import pandas as pd
from typing import TYPE_CHECKING
from tqdm import tqdm

from ..core.pipeline import pipeline_step

from ..core.codon import CodonTable

__all__ = ["convert_codon_to_amino_acid"]


@pipeline_step
def convert_codon_to_amino_acid(
    dataset: pd.DataFrame,
    codon_column: str = "codon_mutations",
    amino_acid_column: str = "aa",
    seq_type: str = "DNA",
    drop_stop: bool = True,
    strict: bool = True,
    drop_codon_column: bool = False,
) -> pd.DataFrame:
    if codon_column not in dataset.columns:
        raise ValueError(f"codon_column '{codon_column}' not in dataset columns")

    tqdm.write("Coverting codons to amino acids...")

    token_re = re.compile(r"^([ACGT]{3})(\d+)([ACGT]{3})$")
    split_re = re.compile(r"[,]+")
    table = CodonTable.get_standard_table(seq_type=seq_type)

    def _convert_field(x):
        if pd.isna(x):
            return pd.NA

        tokens = [t for t in split_re.split(str(x).strip()) if t]
        aa_tokens = []

        for tok in tokens:
            m = token_re.match(tok)
            if not m:
                if strict:
                    raise ValueError(f"Bad COD token: {tok}")
                return pd.NA

            wt_codon, pos_str, mut_codon = m.group(1), m.group(2), m.group(3)

            wt_aa = table.translate_codon(wt_codon)
            mut_aa = table.translate_codon(mut_codon)

            # unknown codon -> 'X'
            if strict and ("X" in (wt_aa, mut_aa)):
                raise ValueError(
                    f"Unknown codon in token: {tok} (wt_aa={wt_aa}, mut_aa={mut_aa})"
                )
            if not strict and ("X" in (wt_aa, mut_aa)):
                return pd.NA

            if drop_stop and (
                table.is_stop_codon(wt_codon) or table.is_stop_codon(mut_codon)
            ):
                return pd.NA

            aa_tokens.append(f"{wt_aa}{int(pos_str)}{mut_aa}")

        return ",".join(aa_tokens)

    out_df = dataset.copy()
    out_df[amino_acid_column] = out_df[codon_column].apply(_convert_field)

    if drop_stop:
        out_df = out_df.dropna(subset=[amino_acid_column]).copy()

    if drop_codon_column:
        out_df = out_df.drop(columns=[codon_column], errors="ignore")

    tqdm.write("Successfully converting codons to amino acids!")

    return out_df
