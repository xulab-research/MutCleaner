import re
import subprocess
import pandas as pd
import mutcleaner.cli.measurements
import mutcleaner.cli.output
import mutcleaner.cli.tables
import mutcleaner.utils.dataset_builders

from datetime import datetime, timezone
from pathlib import Path


def sequence_settings(sequence_type: str, include_ambiguous: bool, include_sec_pyl: bool) -> tuple:
    if sequence_type == "protein":
        sequence_class = mutcleaner.sequence.ProteinSequence
        mutation_type = mutcleaner.mutation.AminoAcidMutation
        alphabet = mutcleaner.alphabet.ProteinAlphabet(include_ambiguous=include_ambiguous, include_sec_pyl=include_sec_pyl, include_stop=False)
    else:
        if include_ambiguous or include_sec_pyl:
            raise ValueError("The special residue flags apply only to protein input")
        mutation_type = mutcleaner.mutation.CodonMutation
        if sequence_type == "dna":
            sequence_class, alphabet = mutcleaner.sequence.DNASequence, mutcleaner.alphabet.DNAAlphabet()
        else:
            sequence_class, alphabet = mutcleaner.sequence.RNASequence, mutcleaner.alphabet.RNAAlphabet()
    return sequence_class, mutation_type, alphabet


def reference_map(reference_rows: pd.DataFrame, id_pattern) -> pd.Series:
    tools = mutcleaner.cli.tables
    identifiers = reference_rows["reference_id"].str.strip()
    if id_pattern is not None:
        identifiers = tools.extract_id(identifiers, id_pattern)
    reference_rows = reference_rows.assign(reference_id=identifiers, reference_sequence=reference_rows["reference_sequence"].str.strip().str.upper())
    reference_rows = reference_rows.loc[reference_rows["reference_id"].ne("") & reference_rows["reference_sequence"].ne("")].drop_duplicates()
    conflicts = reference_rows.loc[reference_rows["reference_id"].duplicated(), "reference_id"].unique()
    if len(conflicts):
        raise ValueError(f"Reference IDs with conflicting reference sequences: {', '.join(conflicts[:10])}")
    return reference_rows.set_index("reference_id")["reference_sequence"]


def reference_slice_bounds(bounds: tuple, sequence_type: str) -> tuple:
    start, stop = bounds
    if not 0 <= start < stop:
        raise ValueError("--reference-slice requires nonempty zero-based START STOP coordinates")
    unit = 1 if sequence_type == "protein" else 3
    return start * unit, stop * unit


def prepare_rows(table: pd.DataFrame, arguments, source: str, reference_column, variant_column, mutation_column: str | None, sequence_class, alphabet, measurement_settings: dict) -> pd.DataFrame:
    tools = mutcleaner.cli.tables
    tools.require_columns(table, [arguments.reference_id_column, reference_column, variant_column, mutation_column])
    empty_cells = pd.Series("", index=table.index, dtype="string")
    identifiers = table[arguments.reference_id_column].str.strip() if arguments.reference_id_column is not None else pd.Series(arguments.reference_id.strip(), index=table.index, dtype="string")
    if arguments.reference_id_pattern is not None:
        if arguments.reference_id_column is None:
            raise ValueError("--reference-id-pattern requires --reference-id-column")
        identifiers = tools.extract_id(identifiers, arguments.reference_id_pattern)

    rows = pd.DataFrame(
        {
            "source_row": range(1, len(table) + 1),
            "reference_id": identifiers,
            "reference_sequence": empty_cells if reference_column is None else table[reference_column].str.strip().str.upper(),
            "variant_sequence": empty_cells if variant_column is None else table[variant_column].str.strip().str.upper(),
            "input_mutation": empty_cells if mutation_column is None else table[mutation_column].str.strip(),
            **{alias: table[column] for column, alias in measurement_settings["labels"].items()},
            **{alias: table[column].str.strip() for column, alias in measurement_settings["conditions"].items()},
        }
    )
    rows["filter_reason"] = tools.filter_reasons(table, arguments)
    rows["error_message"] = ""
    invalid_id = rows["reference_id"].eq("") | rows["reference_id"].str.contains(r"[\x00-\x1f\x7f]")
    marked_wild_type = tools.wild_type_mask(rows["input_mutation"], arguments.wild_type_marker)
    if source in {"column", "wild_type_rows"}:
        sequence_column = "reference_sequence" if source == "column" else "variant_sequence"
        selected = ~invalid_id & rows["filter_reason"].eq("")
        if source == "wild_type_rows":
            selected &= marked_wild_type
            if not (selected & rows["variant_sequence"].ne("")).any():
                raise ValueError("--reference-from-wt found no usable reference sequences with the specified WT marker after filtering")
        reference_rows = rows.loc[selected, ["reference_id", sequence_column]].rename(columns={sequence_column: "reference_sequence"})
        rows["reference_sequence"] = rows["reference_id"].map(reference_map(reference_rows, None)).fillna("")
    elif source == "fasta":
        references = mutcleaner.sequence.load_sequences_from_fasta(arguments.reference_fasta, sequence_class, alphabet=alphabet)
        reference_rows = pd.DataFrame([(identifier, str(sequence)) for identifier, sequence in references.items()], columns=["reference_id", "reference_sequence"], dtype="string")
        rows["reference_sequence"] = rows["reference_id"].map(reference_map(reference_rows, arguments.reference_file_id_pattern)).fillna("")
    elif source == "reference_table":
        reference_rows = tools.read_reference_table(arguments.reference_table, arguments.reference_table_id_column, arguments.reference_table_sequence_column)
        rows["reference_sequence"] = rows["reference_id"].map(reference_map(reference_rows, arguments.reference_file_id_pattern)).fillna("")
    elif source == "literal":
        rows["reference_sequence"] = str(sequence_class(arguments.reference_sequence.strip().upper(), alphabet=alphabet))

    if arguments.reference_slice is not None:
        if source == "inferred":
            raise ValueError("--reference-slice requires a supplied reference")
        start, stop = reference_slice_bounds(arguments.reference_slice, arguments.sequence_type)
        too_short = rows["reference_sequence"].ne("") & rows["reference_sequence"].str.len().lt(stop)
        rows.loc[too_short, "error_message"] = "Reference is shorter than the requested slice"
        rows["reference_sequence"] = rows["reference_sequence"].str.slice(start, stop)
        if source == "wild_type_rows":
            rows.loc[marked_wild_type, "variant_sequence"] = rows.loc[marked_wild_type, "variant_sequence"].str.slice(start, stop)

    if variant_column is not None:
        rows.loc[~marked_wild_type & rows["variant_sequence"].eq(""), "error_message"] = "Variant sequence is missing"
    if mutation_column is not None:
        rows.loc[rows["input_mutation"].eq(""), "error_message"] = "Mutation annotation is missing"
    if source != "inferred":
        mismatched_wild_type = marked_wild_type & rows["variant_sequence"].ne("") & rows["variant_sequence"].ne(rows["reference_sequence"])
        rows.loc[mismatched_wild_type, "error_message"] = "Marked wild-type sequence differs from the reference"
        rows.loc[rows["reference_sequence"].eq("") & rows["error_message"].eq(""), "error_message"] = "No reference sequence found for this reference ID"
    if arguments.sequence_type != "protein":
        wrong_length = rows["reference_sequence"].str.len().mod(3).ne(0) | rows["variant_sequence"].str.len().mod(3).ne(0)
        rows.loc[wrong_length & rows["error_message"].eq(""), "error_message"] = "Codon sequences must contain a whole number of codons"
    for column in measurement_settings["baseline_columns"]:
        alias = measurement_settings["conditions"][column]
        rows.loc[rows[alias].eq("") & rows["error_message"].eq(""), "error_message"] = f"Required grouping condition is missing: {column}"
    rows.loc[invalid_id, "error_message"] = "Reference ID is missing or invalid"
    rows["is_marked_wild_type"] = marked_wild_type
    rows["is_unchanged"] = marked_wild_type if mutation_column is not None else rows["variant_sequence"].eq(rows["reference_sequence"])
    return rows


def read_position_rules(arguments) -> tuple[dict, dict]:
    tools = mutcleaner.cli.tables
    offsets = {identifier: int(offset) for identifier, offset in tools.parse_assignments(arguments.offset_for, "--offset-for").items()}
    position_maps = {}
    if arguments.position_map is not None:
        rows, _ = tools.read_table(arguments.position_map, None, None, None, None, list_sheets=False)
        columns = ["reference_id", "input_position", "reference_index"]
        tools.require_columns(rows, columns)
        for reference_id, source_position, target_position in rows[columns].itertuples(index=False, name=None):
            mapping = position_maps.setdefault(reference_id, {})
            source_position, target_position = int(source_position), int(target_position)
            if not reference_id or source_position < 0 or target_position < 0 or source_position in mapping:
                raise ValueError("Position maps require nonnegative positions and unique reference/source pairs")
            mapping[source_position] = target_position
    return offsets, position_maps


def remap_annotation(annotation: str, separator: str, offset: int | None, position_map: dict) -> str:
    remapped = []
    for token in annotation.split(separator):
        match = re.fullmatch(r"([A-Za-z*]{1,3})(\d+)([A-Za-z*]{1,3})", token.strip())
        if match is None:
            raise ValueError("Annotation contains a token that could not be mapped")
        wild_type, position, mutant = match.groups()
        position = int(position)
        target = position_map.get(position, position - offset if offset is not None else -1)
        if target < 0:
            raise ValueError("Annotation contains a position that could not be mapped")
        remapped.append(f"{wild_type}{target}{mutant}")
    return separator.join(remapped)


def normalize_annotations(candidates: pd.DataFrame, arguments, mutation_type, alphabet, offsets: dict, position_maps: dict) -> tuple[pd.DataFrame, list]:
    cleaners = mutcleaner.basic_cleaners
    separator = arguments.mutation_separator
    if arguments.offset is not None and not offsets and not position_maps:
        result, artifacts = cleaners.validate_mutations(candidates, mutation_column="mutation_name", mutation_sep=separator, mutation_type=mutation_type, alphabet=alphabet, offset=arguments.offset, num_workers=arguments.workers)
        return result, [artifacts["failed"]]

    candidates = candidates.copy()
    candidates["normalization_error"] = ""
    candidates["annotation_offset"] = pd.Series(pd.NA, index=candidates.index, dtype="Int64")
    missing_rules = []
    for reference_id, group in candidates.groupby("reference_id", sort=False):
        offset = offsets.get(reference_id, arguments.offset)
        position_map = position_maps.get(reference_id, {})
        if not position_map and offset is not None:
            candidates.loc[group.index, "annotation_offset"] = offset
            continue

        replacements = {}
        errors = {}
        missing_positions = set()
        affected_rows = 0
        for annotation, count in group["mutation_name"].value_counts(sort=False).items():
            try:
                if offset is None:
                    mutations = mutcleaner.mutation.MutationSet.from_string(annotation, sep=separator, offset=0, mutation_type=mutation_type, alphabet=alphabet)
                    uncovered = {mutation.position for mutation in mutations} - position_map.keys()
                    if uncovered:
                        missing_positions.update(uncovered)
                        affected_rows += count
                        continue
                replacements[annotation] = remap_annotation(annotation, separator, offset, position_map)
            except ValueError as error:
                errors[annotation] = str(error)
        if affected_rows:
            examples = ", ".join(map(str, sorted(missing_positions)[:8]))
            if len(missing_positions) > 8:
                examples += ", ..."
            missing_rules.append(f"{reference_id!r}: {affected_rows} row(s), input positions {examples}")
        candidates.loc[group.index, "mutation_name"] = group["mutation_name"].map(replacements).fillna(group["mutation_name"])
        candidates.loc[group.index, "normalization_error"] = group["mutation_name"].map(errors).fillna("")
        candidates.loc[group.index, "annotation_offset"] = 0

    if missing_rules:
        raise ValueError("Missing numbering rules. Specify --offset, --offset-for or complete --position-map entries for " + "; ".join(missing_rules))

    failed = candidates.loc[candidates["normalization_error"].ne("")]
    failures = [failed.assign(error_message=failed["normalization_error"]).drop(columns=["annotation_offset", "normalization_error"])]
    successful = [candidates.iloc[:0].drop(columns=["annotation_offset", "normalization_error"])]
    valid = candidates.loc[candidates["normalization_error"].eq("")]
    for offset, group in valid.groupby("annotation_offset", sort=False):
        normalized, artifacts = cleaners.validate_mutations(group.drop(columns=["annotation_offset", "normalization_error"]), mutation_column="mutation_name", mutation_sep=separator, mutation_type=mutation_type, alphabet=alphabet, offset=int(offset), num_workers=arguments.workers)
        successful.append(normalized)
        failures.append(artifacts["failed"])
    return pd.concat(successful, ignore_index=True).sort_values("source_row"), failures


def resolve_mutations(candidates: pd.DataFrame, arguments, workflow: str, mutation_type, alphabet, offsets, position_maps) -> tuple[pd.DataFrame, pd.DataFrame]:
    cleaners = mutcleaner.basic_cleaners
    accepted = candidates.assign(mutation_name=candidates["input_mutation"].str.upper())
    empty_failures = accepted.iloc[:0].assign(error_message="")
    if accepted.empty:
        return accepted, empty_failures
    if workflow == "compare_sequences":
        accepted, artifacts = cleaners.infer_mutations_from_sequences(accepted, wt_sequence_column="reference_sequence", mut_sequence_column="variant_sequence", sequence_type=arguments.sequence_type, alphabet=alphabet, num_workers=arguments.workers)
        accepted["mutation_name"] = accepted.pop("inferred_mutations")
        return accepted, artifacts["failed"]

    accepted, failures = normalize_annotations(accepted, arguments, mutation_type, alphabet, offsets, position_maps)
    if accepted.empty:
        return accepted, pd.concat(failures, ignore_index=True)

    if workflow == "infer_reference":
        inferred, artifacts = cleaners.infer_wildtype_sequences(accepted, name_column="reference_id", mutation_column="mutation_name", sequence_column="variant_sequence", label_columns=[], sequence_type=arguments.sequence_type, alphabet=alphabet, offset=0, num_workers=arguments.workers)
        group_failures = artifacts["failed"]
        if not group_failures.empty:
            errors_by_reference = group_failures.set_index("reference_id")["error_message"]
            failed = accepted["reference_id"].isin(errors_by_reference.index)
            failed_rows = accepted.loc[failed]
            failures.append(failed_rows.assign(error_message=failed_rows["reference_id"].map(errors_by_reference)))
            accepted = accepted.loc[~failed]
        if not accepted.empty:
            references = inferred.loc[inferred["mutation_name"].eq("WT")].set_index("reference_id")["variant_sequence"]
            accepted = accepted.assign(reference_sequence=accepted["reference_id"].map(references))
    else:
        accepted, artifacts = cleaners.apply_mutations_to_sequences(accepted, sequence_column="reference_sequence", mutated_sequence_column="variant_sequence", name_column="reference_id", mutation_column="mutation_name", mutation_type=mutation_type, sequence_type=arguments.sequence_type, alphabet=alphabet, offset=0, num_workers=arguments.workers)
        failures.append(artifacts["failed"])

    return accepted, pd.concat(failures, ignore_index=True)


def filter_mutations(accepted: pd.DataFrame, arguments, mutation_type, alphabet) -> tuple[pd.DataFrame, pd.DataFrame]:
    reasons = pd.Series("", index=accepted.index, dtype="string")
    depth = accepted["mutation_name"].astype("string").str.count(",") + 1
    if arguments.min_substitutions is not None:
        reasons.loc[depth.lt(arguments.min_substitutions)] = "Below minimum substitution count"
    if arguments.max_substitutions is not None:
        reasons.loc[depth.gt(arguments.max_substitutions)] = "Above maximum substitution count"
    if arguments.filter_stop or arguments.mutation_effect:
        codon_table = mutcleaner.codon.CodonTable.get_standard_table(arguments.sequence_type.upper()) if arguments.sequence_type != "protein" else None
        effect_reasons = {}
        for annotation in accepted["mutation_name"].unique():
            mutations = mutcleaner.mutation.MutationSet.from_string(annotation, sep=",", offset=0, mutation_type=mutation_type, alphabet=alphabet)
            if arguments.filter_stop and any(codon_table.is_stop_codon(mutation.mutant_codon) for mutation in mutations):
                effect_reasons[annotation] = "Introduces a stop codon"
            elif arguments.mutation_effect:
                amino_acid_mutations = mutations if codon_table is None else (mutation.to_amino_acid_mutation(codon_table) for mutation in mutations)
                if any(mutation.get_mutation_category() not in arguments.mutation_effect for mutation in amino_acid_mutations):
                    effect_reasons[annotation] = "Substitution effect is outside the requested set"
        reasons = reasons.mask(reasons.eq(""), accepted["mutation_name"].map(effect_reasons).fillna(""))
    return accepted.loc[reasons.eq("")], accepted.loc[reasons.ne("")].assign(filter_reason=reasons.loc[reasons.ne("")])


def translate_rows(accepted: pd.DataFrame, sequence_class, workers: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    unique_sequences = pd.concat([accepted["reference_sequence"], accepted["variant_sequence"]]).unique()
    translations = {sequence: str(sequence_class(sequence).translate(require_mod3=True)) for sequence in unique_sequences}
    translated = accepted.assign(reference_sequence=accepted["reference_sequence"].map(translations), variant_sequence=accepted["variant_sequence"].map(translations))
    synonymous = translated["reference_sequence"].eq(translated["variant_sequence"])
    filtered = translated.loc[synonymous].assign(filter_reason="No amino acid change after translation")
    variants = translated.loc[~synonymous]
    if variants.empty:
        return variants, filtered, variants.assign(error_message="")
    protein_rows, artifacts = mutcleaner.basic_cleaners.infer_mutations_from_sequences(variants, wt_sequence_column="reference_sequence", mut_sequence_column="variant_sequence", num_workers=workers)
    protein_rows["mutation_name"] = protein_rows.pop("inferred_mutations")
    return protein_rows, filtered, artifacts["failed"]


def clean_table(table: pd.DataFrame, arguments) -> dict:

    #######################################################################
    # validate options and determine the workflow
    #######################################################################

    started_at = datetime.now(timezone.utc).isoformat()
    if arguments.output_path.exists() or arguments.output_path.is_symlink():
        raise FileExistsError(f"Output path already exists: {arguments.output_path}")
    tools = mutcleaner.cli.tables
    measurements = mutcleaner.cli.measurements
    steps = ["read"]
    custom_statistics = custom_source = None
    if custom_paths := getattr(arguments, "custom_cleaner", ()):
        table, custom_statistics, custom_source = tools.apply_custom_cleaner(table, custom_paths[0])
        steps.append("custom")
    table, mutation_column = tools.preprocess(table, arguments)
    if any((arguments.rename, arguments.constant, arguments.join_column, arguments.split_column, arguments.replace, arguments.extract, arguments.mutation_parts)):
        steps.append("prepare_fields")
    reference_column = arguments.reference_sequence_column
    variant_column = arguments.variant_sequence_column
    reference_sources = {
        "column": reference_column is not None,
        "fasta": arguments.reference_fasta is not None,
        "reference_table": arguments.reference_table is not None,
        "literal": arguments.reference_sequence is not None,
        "wild_type_rows": arguments.reference_from_wt,
        "inferred": arguments.infer_reference,
    }
    selected_sources = [source for source, selected in reference_sources.items() if selected]
    if arguments.reference_table is None and (arguments.reference_table_id_column is not None or arguments.reference_table_sequence_column is not None):
        raise ValueError("Reference table columns require --reference-table")
    if arguments.reference_file_id_pattern is not None and arguments.reference_fasta is None and arguments.reference_table is None:
        raise ValueError("--reference-file-id-pattern requires --reference-fasta or --reference-table")
    if mutation_column is None:
        if arguments.offset is not None or arguments.offset_for or arguments.position_map is not None or arguments.wild_type_marker is not None or arguments.no_wild_type_marker:
            raise ValueError("Numbering and wild-type marker options require mutation annotations")
        if variant_column is None:
            raise ValueError("Sequence comparison requires --variant-sequence-column")
    else:
        if not arguments.mutation_separator:
            raise ValueError("Mutation annotations require a nonempty --mutation-separator")
        if (arguments.wild_type_marker is not None) == arguments.no_wild_type_marker:
            raise ValueError("Mutation annotations require exactly one of --wild-type-marker or --no-wild-type-marker")
        if arguments.wild_type_marker is not None and not arguments.wild_type_marker.strip():
            raise ValueError("--wild-type-marker must be nonempty")

    if len(selected_sources) != 1:
        raise ValueError("Specify exactly one reference source: --reference-sequence-column, --reference-fasta, --reference-table, --reference-sequence, --reference-from-wt, or --infer-reference")
    reference_source = selected_sources[0]
    if reference_source == "literal":
        if arguments.reference_id is None or not arguments.reference_sequence.strip():
            raise ValueError("--reference-sequence requires --reference-id and a nonempty sequence")
    if reference_source == "inferred" and (mutation_column is None or variant_column is None):
        raise ValueError("--infer-reference requires mutation annotations and --variant-sequence-column")
    if reference_source == "wild_type_rows" and (arguments.wild_type_marker is None or variant_column is None):
        raise ValueError("--reference-from-wt requires --wild-type-marker and --variant-sequence-column")
    if arguments.to_protein and arguments.sequence_type == "protein":
        raise ValueError("--to-protein requires --sequence-type dna or rna")
    if "synonymous" in arguments.mutation_effect and arguments.sequence_type == "protein":
        raise ValueError("Synonymous substitutions can only be filtered from codon inputs")
    if arguments.filter_stop and arguments.sequence_type == "protein":
        raise ValueError("--filter-stop applies to DNA or RNA codon inputs")
    if arguments.min_substitutions is not None and arguments.max_substitutions is not None and arguments.min_substitutions > arguments.max_substitutions:
        raise ValueError("Minimum substitutions must not exceed maximum substitutions")

    if mutation_column is None:
        workflow = "compare_sequences"
    elif reference_source == "inferred":
        workflow = "infer_reference"
    elif variant_column is None:
        workflow = "apply_annotations"
    else:
        workflow = "cross_check"

    settings = measurements.measurement_settings(table, arguments, mutation_column)
    sequence_options = {key: getattr(arguments, key) for key in ("sequence_type", "include_ambiguous", "include_sec_pyl")}
    sequence_class, mutation_type, alphabet = sequence_settings(**sequence_options)
    offsets, position_maps = read_position_rules(arguments)

    #######################################################################
    # process records and prepare measurement tasks
    #######################################################################

    rows = prepare_rows(table, arguments, reference_source, reference_column, variant_column, mutation_column, sequence_class, alphabet, settings)
    if arguments.filter or arguments.exclude or arguments.drop_missing:
        steps.append("filter_input")
    filtered = rows.loc[rows["filter_reason"].ne("")]
    eligible = rows["filter_reason"].eq("") & rows["error_message"].eq("")
    unchanged = rows.loc[eligible & rows["is_unchanged"]]
    candidates = rows.loc[eligible & ~rows["is_unchanged"]].drop(columns="error_message")
    accepted, mutation_failures = resolve_mutations(candidates, arguments, workflow, mutation_type, alphabet, offsets, position_maps)
    if not candidates.empty:
        steps.append("resolve_references_and_variants")
    elif reference_source != "inferred":
        steps.append("resolve_references")
    failures = pd.concat([rows.loc[rows["filter_reason"].eq("") & rows["error_message"].ne("")], mutation_failures], ignore_index=True)
    if workflow == "infer_reference" and not unchanged.empty:
        inferred = accepted.drop_duplicates("reference_id").set_index("reference_id")["reference_sequence"]
        expected = unchanged["reference_id"].map(inferred)
        inconsistent_wt = expected.notna() & unchanged["variant_sequence"].ne("") & unchanged["variant_sequence"].ne(expected)
        failures = pd.concat([failures, unchanged.loc[inconsistent_wt].assign(error_message="Marked wild-type sequence differs from the inferred reference")], ignore_index=True)
        unchanged = unchanged.loc[~inconsistent_wt]
    if not accepted.empty and any((arguments.min_substitutions, arguments.max_substitutions, arguments.mutation_effect, arguments.filter_stop)):
        accepted, mutation_filtered = filter_mutations(accepted, arguments, mutation_type, alphabet)
        filtered = pd.concat([filtered, mutation_filtered], ignore_index=True)
        steps.append("filter_mutations")
    output_sequence_options = sequence_options
    output_sequence_class, output_mutation_type, output_alphabet = sequence_class, mutation_type, alphabet
    if arguments.to_protein:
        output_sequence_options = {"sequence_type": "protein", "include_ambiguous": False, "include_sec_pyl": False}
        output_sequence_class, output_mutation_type, output_alphabet = sequence_settings(**output_sequence_options)
        if not accepted.empty:
            accepted, translation_filtered, translation_failures = translate_rows(accepted, sequence_class, arguments.workers)
            filtered = pd.concat([filtered, translation_filtered], ignore_index=True)
            failures = pd.concat([failures, translation_failures], ignore_index=True)
            steps.append("translate_codons")
    tasks, task_summaries, measurement_steps = measurements.measurement_tasks(accepted, unchanged, arguments, settings)
    steps.extend(measurement_steps)

    for task in tasks:
        converted, references = mutcleaner.utils.dataset_builders.convert_format_2(task["rows"], name_column="reference_id", mutation_column="mutation_name", sequence_column="reference_sequence", score_column="label", mutation_set_prefix="set", offset=0, additional_metadata=None, sequence_class=output_sequence_class, mutation_type=output_mutation_type, alphabet=output_alphabet)
        dataset = mutcleaner.MutationDataset.from_dataframe(converted, references, name=arguments.input_path.stem)
        validation = dataset.validate_against_references()
        if validation["invalid_mutation_sets"] or validation["position_mismatches"]:
            raise RuntimeError("Prepared mutations failed reference validation")
        task["dataset"] = dataset

    #######################################################################
    # record effective parameters and publish verified datasets
    #######################################################################

    output = mutcleaner.cli.output
    unassigned = failures["reference_id"].eq("") | failures["reference_id"].str.contains(r"[\x00-\x1f\x7f]")
    parameters = {key: str(value.resolve()) if isinstance(value, Path) else value for key, value in vars(arguments).items() if key not in {"command_options", "read_settings", "custom_cleaner"}}
    parameters.update(reference_sequence_column=reference_column, variant_sequence_column=variant_column, mutation_column=mutation_column)
    if arguments.wild_type_marker is not None:
        parameters["wild_type_marker"] = arguments.wild_type_marker.strip()
    command_options = dict(arguments.command_options)
    if reference_column is not None:
        command_options["--reference-sequence-column"] = [(reference_column,)]
    if variant_column is not None:
        command_options["--variant-sequence-column"] = [(variant_column,)]
    if arguments.wild_type_marker is not None:
        command_options["--wild-type-marker"] = [(parameters["wild_type_marker"],)]
    command_arguments = ["clean"]
    for flag, occurrences in command_options.items():
        for values in occurrences:
            command_arguments.extend((flag, *(str(value) for value in values)))
    command_arguments.extend(("--", str(arguments.input_path), str(arguments.output_path)))
    source_root_path = Path(mutcleaner.__file__).resolve().parents[1]
    git_commit = None
    if (source_root_path / ".git").exists():
        try:
            result = subprocess.run(["git", "-C", str(source_root_path), "rev-parse", "HEAD"], capture_output=True, text=True, check=False, timeout=2)
            if result.returncode == 0:
                git_commit = result.stdout.strip()
        except (OSError, subprocess.SubprocessError):
            pass

    report = {
        "started_at": started_at,
        "mutcleaner_version": mutcleaner.__version__,
        "git_commit": git_commit,
        "workflow": workflow,
        "command_arguments": command_arguments,
        "read_settings": arguments.read_settings,
        "effective_parameters": parameters,
        "reference_source": reference_source,
        "input_position_unit": "residue" if arguments.sequence_type == "protein" else "codon",
        "output_sequence_type": output_sequence_options["sequence_type"],
        "position_unit": "residue" if output_sequence_options["sequence_type"] == "protein" else "codon",
        "output_offset": 0,
        "steps": steps,
        "wild_type_policy": "subtract measured WT before variant aggregation" if arguments.subtract_wt else "count and omit WT from the dataset",
        "label_formula": "measurement - matched WT baseline" if arguments.subtract_wt else "measurement",
        "nearest_distance": "lexicographic" if settings["targets"] else None,
        "global_statistics": {**output.summarize_rows(rows, failures, unchanged, filtered, accepted), "unassigned_invalid_rows": int(unassigned.sum()), "finite_wild_type_measurements": {column: int(measurements.finite_numbers(unchanged[alias]).notna().sum()) for column, alias in settings["labels"].items()}},
        "measurement_tasks": task_summaries,
    }
    if custom_statistics is not None:
        report["custom"] = custom_statistics
    return output.write_results(arguments.output_path, tasks, report, (rows, failures, unchanged, filtered, accepted), output_sequence_options, arguments.save_format, command_options, custom_source)
