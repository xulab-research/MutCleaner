import pandas as pd


def expand_mutations(
    dataset: pd.DataFrame,
    mutation_column: str = "mut_info",
    label_column: str = "fitness",
    sep: str = ",",
) -> pd.DataFrame:
    """Expand single and double mutation records into a long-format table.

    Each input row contains two single mutations (described by separate
    wild-type amino acid, position, and mutant amino acid columns) along
    with their individual fitness values and a combined double-mutation
    fitness value. This function unpacks each row into three records: the
    first single mutation, the second single mutation, and the double
    mutation, each paired with its corresponding label value.

    Parameters
    ----------
    dataset : pandas.DataFrame
        Input DataFrame. Must contain the columns ``WT_AA1``, ``Pos1``,
        ``Mut1``, ``WT_AA2``, ``Pos2``, ``Mut2``, ``fitness1``,
        ``fitness2``, and ``fitness``.
    mutation_column : str, optional
        Name of the output column holding the mutation descriptor.
        Defaults to ``"mut_info"``.
    label_column : str, optional
        Name of the output column holding the label (fitness) value.
        Defaults to ``"fitness"``.
    sep : str, optional
        Separator used to join the two single mutations into a double
        mutation descriptor. Defaults to ``","``.

    Returns
    -------
    pandas.DataFrame
        Long-format DataFrame with two columns named according to
        ``mutation_column`` and ``label_column``.
    """
    single1 = dataset["WT_AA1"].astype(str) + dataset["Pos1"].astype(str) + dataset["Mut1"].astype(str)
    single2 = dataset["WT_AA2"].astype(str) + dataset["Pos2"].astype(str) + dataset["Mut2"].astype(str)
    double = single1 + sep + single2

    part1 = pd.DataFrame({mutation_column: single1, label_column: dataset["fitness1"]})
    part2 = pd.DataFrame({mutation_column: single2, label_column: dataset["fitness2"]})
    part3 = pd.DataFrame({mutation_column: double, label_column: dataset["fitness"]})

    return pd.concat([part1, part2, part3], ignore_index=True)


def map_ambler_positions(
    dataset: pd.DataFrame,
    mapping: dict[str, str],
    mutation_column: str = "mut_info",
    sep: str = ",",
) -> pd.DataFrame:
    """Remap mutation positions from Ambler numbering to sequence numbering.

    The mutation descriptors in the dataset use Ambler numbering, which is
    not aligned with the true amino acid sequence numbering. This function
    rewrites each mutation descriptor by replacing its wild-type-plus-position
    prefix (e.g. ``"S70"``) according to the provided mapping, while keeping
    the trailing mutant amino acid unchanged. Both single and multi-mutation
    records (joined by ``sep``) are handled.

    Parameters
    ----------
    dataset : pandas.DataFrame
        Input DataFrame containing the mutation descriptor column.
    mapping : dict of str to str
        Mapping from Ambler-numbered prefixes to sequence-numbered prefixes
        (e.g. ``{"S70": "S66"}``). Prefixes consist of the wild-type amino
        acid followed by the position. Prefixes absent from the mapping are
        left unchanged.
    mutation_column : str, optional
        Name of the column holding the mutation descriptor.
        Defaults to ``"mut_info"``.
    sep : str, optional
        Separator splitting multiple mutations within a single descriptor.
        Defaults to ``","``.

    Returns
    -------
    pandas.DataFrame
        The input DataFrame with ``mutation_column`` updated in place to use
        the remapped positions.
    """
    dataset[mutation_column] = dataset[mutation_column].str.split(sep).apply(lambda muts: sep.join(mapping.get(m[:-1], m[:-1]) + m[-1] for m in muts))
    return dataset
