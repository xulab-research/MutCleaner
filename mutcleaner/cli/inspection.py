import json
import math
import re
import pandas as pd
import mutcleaner.cli.clean
import mutcleaner.cli.measurements
import mutcleaner.cli.output
import mutcleaner.cli.tables
import mutcleaner.utils.mutation_converter

from collections import Counter
from functools import lru_cache
from pathlib import Path


def read_fasta_records(fasta_path: Path) -> list:
    records = []
    identifier = None
    lines = []
    with fasta_path.open(encoding="utf-8-sig") as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                if identifier is not None:
                    records.append((identifier, "".join(lines)))
                header = line[1:].split(maxsplit=1)
                if not header:
                    raise ValueError("FASTA record has an empty header")
                identifier = header[0]
                lines = []
            elif line:
                if identifier is None:
                    raise ValueError("FASTA sequence appeared before the first header")
                lines.append("".join(line.split()).upper())
    if identifier is not None:
        records.append((identifier, "".join(lines)))
    if not records:
        raise ValueError("No records found in the FASTA file")
    return records


def describe_fasta(fasta_path: Path, sample_rows: int) -> dict:
    records = read_fasta_records(fasta_path)
    counts = Counter(identifier for identifier, _ in records)
    characters = Counter(character for _, sequence in records for character in sequence)
    lengths = [len(sequence) for _, sequence in records]
    return {"kind": "fasta", "input_path": str(fasta_path.resolve()), "records": len(records), "distinct_ids": len(counts), "duplicate_ids": sum(count > 1 for count in counts.values()), "length_min": min(lengths), "length_max": max(lengths), "characters": dict(sorted(characters.items())), "samples": [{"id": identifier, "length": len(sequence)} for identifier, sequence in records[:sample_rows]]}


def dataset_directories(input_path: Path) -> list:
    required = {"data.csv", "wt.fasta", "metadata.json"}
    cleaner_header = (mutcleaner.cli.output.CLEANER_HEADER.format(save_format="mutcleaner") + "\n").encode("ascii")
    directories = [(input_path, False)]
    candidates = []
    while directories:
        current_path, is_task = directories.pop()
        files = {entry.name for entry in current_path.iterdir() if entry.is_file()}
        present = required.intersection(files)
        candidate = len(present) >= 2 or "wt.fasta" in present
        if not candidate and "metadata.json" in present:
            try:
                metadata = json.loads((current_path / "metadata.json").read_text(encoding="utf-8"))
                candidate = isinstance(metadata, dict) and "reference_id" in metadata
            except (OSError, ValueError):
                candidate = False
        if not candidate and "data.csv" in present:
            try:
                header = pd.read_csv(current_path / "data.csv", nrows=0)
                candidate = {"mutation_name", "mutated_sequence", "label"}.issubset(header.columns)
            except (OSError, ValueError):
                candidate = False
        if not candidate and "cleaner.py" in files:
            try:
                with (current_path / "cleaner.py").open("rb") as script:
                    candidate = script.read(len(cleaner_header)) == cleaner_header
            except OSError:
                candidate = False
        if candidate:
            candidates.append((current_path, sorted(required - files)))
        elif current_path == input_path or is_task:
            for child in sorted(current_path.iterdir(), reverse=True):
                if child.is_dir() and not child.is_symlink():
                    is_child_task = re.fullmatch(r"(?:label|condition)_\d{3,}", child.name) is not None
                    directories.append((child, is_child_task))
    return candidates


def verify_reference(reference_path: Path, metadata: dict, variant_rows: pd.DataFrame) -> dict:
    sequence_type = {"ProteinSequence": "protein", "DNASequence": "dna", "RNASequence": "rna"}.get(metadata.get("sequence_type"))
    if sequence_type is None:
        raise ValueError(f"Unrecognized sequence type: {metadata.get('sequence_type')!r}")
    sequence_class, mutation_type, alphabet = mutcleaner.cli.clean.sequence_settings(sequence_type, sequence_type == "protein", sequence_type == "protein")
    records = read_fasta_records(reference_path / "wt.fasta")
    if len(records) != 1:
        raise ValueError("Expected one FASTA reference")
    reference = sequence_class(records[0][1], alphabet=alphabet)
    failures = 0
    covered = set()
    for mutation_name, variant_sequence in variant_rows[["mutation_name", "mutated_sequence"]].itertuples(index=False, name=None):
        try:
            mutations = mutcleaner.mutation.MutationSet.from_string(mutation_name, sep=",", offset=0, mutation_type=mutation_type, alphabet=alphabet)
            if str(reference.apply_mutation(mutations)) == variant_sequence:
                covered.update(mutation.position for mutation in mutations)
            else:
                failures += 1
        except (TypeError, ValueError, IndexError):
            failures += 1
    unit_size = 1 if sequence_type == "protein" else 3
    position_unit = "residue" if unit_size == 1 else "codon"
    length = len(reference) // unit_size
    coverage = len(covered) / length * 100 if length else 0.0
    stored_coverage = metadata.get("coverage_percentage")
    metadata_matches = len(reference) % unit_size == 0 and len(variant_rows) == metadata.get("num_mutation_sets") and len(reference) == metadata.get("sequence_length") and len(covered) == metadata.get("covered_positions") and metadata.get("position_unit") == position_unit and isinstance(stored_coverage, (int, float)) and not isinstance(stored_coverage, bool) and math.isclose(coverage, stored_coverage, rel_tol=1e-12, abs_tol=1e-12)
    invalid_labels = int(mutcleaner.cli.measurements.finite_numbers(variant_rows["label"]).isna().sum())
    return {
        "status": "passed" if not failures and not invalid_labels and metadata_matches else "failed",
        "checked_rows": len(variant_rows),
        "sequence_failures": failures,
        "invalid_labels": invalid_labels,
        "covered_positions": len(covered),
        "coverage_percentage": coverage,
        "position_unit": position_unit,
        "metadata_matches": metadata_matches,
    }


def describe_directory(input_path: Path, verify: bool) -> dict:
    entries = sorted(input_path.iterdir())
    candidates = dataset_directories(input_path)
    report = {
        "kind": "directory",
        "input_path": str(input_path.resolve()),
        "entries": len(entries),
        "inventory": [{"name": entry.name, "kind": "directory" if entry.is_dir() else "file", "format": entry.suffix.lower().lstrip(".") if entry.is_file() else None} for entry in entries[:50]],
        "inventory_truncated": len(entries) > 50,
        "datasets": [],
    }
    for reference_path, missing in candidates:
        item = {"directory": reference_path.name, "relative_path": str(reference_path.relative_to(input_path)), "missing_files": missing, "verification": None}
        if missing:
            if verify:
                item["verification"] = {"status": "failed", "error": "Incomplete reference directory"}
        else:
            try:
                metadata = json.loads((reference_path / "metadata.json").read_text(encoding="utf-8"))
                variant_rows = pd.read_csv(reference_path / "data.csv", dtype="string", keep_default_na=False)
                if not isinstance(metadata, dict):
                    raise ValueError("Reference metadata must be a JSON object")
                coverage = metadata.get("coverage_percentage")
                coverage = coverage if isinstance(coverage, (int, float)) and math.isfinite(coverage) else None
                positions = metadata.get("covered_positions")
                positions = positions if type(positions) is int else None
                item.update(reference_id=str(metadata["reference_id"]), rows=len(variant_rows), duplicate_mutation_rows=int(variant_rows["mutation_name"].duplicated().sum()), covered_positions=positions, coverage_percentage=coverage)
                if verify:
                    item["verification"] = verify_reference(reference_path, metadata, variant_rows)
            except (KeyError, OSError, TypeError, ValueError) as error:
                item["error"] = str(error)
                if verify:
                    item["verification"] = {"status": "failed", "error": str(error)}
        report["datasets"].append(item)
    if verify:
        report["verification"] = "passed" if candidates and all(item["verification"]["status"] == "passed" for item in report["datasets"]) else "failed"
        if not candidates:
            report["error"] = "No standard reference datasets found"
    return report


def describe_table(table: pd.DataFrame, sample_rows: int, profile_rows: int) -> dict:
    profile = table if profile_rows >= len(table) else table.head(profile_rows)
    columns = {}
    for column in profile:
        text = profile[column].str.strip()
        values = text.loc[text.ne("")]
        numbers = mutcleaner.cli.measurements.finite_numbers(values).dropna()
        median = None
        if len(numbers):
            ordered = numbers.sort_values()
            low, high = float(ordered.iloc[(len(ordered) - 1) // 2]), float(ordered.iloc[len(ordered) // 2])
            median = low if low == high else low / 2 + high / 2
        observed_type = "empty" if values.empty else ("numeric" if len(values) == len(numbers) else ("mixed" if len(numbers) else "text"))
        columns[column] = {"observed_type": observed_type, "empty": int(text.eq("").sum()), "unique_nonempty": int(values.nunique()), "finite_numeric": len(numbers), "non_finite_or_non_numeric": len(values) - len(numbers), "min": float(numbers.min()) if len(numbers) else None, "max": float(numbers.max()) if len(numbers) else None, "median": median}
    sample = [{column: value[:180] for column, value in row.items()} for row in profile.head(sample_rows).to_dict(orient="records")]
    return {"kind": "table", "row_count": len(table), "column_count": len(table.columns), "columns": table.columns.tolist(), "profiled_rows": len(profile), "full_table_profiled": len(profile) == len(table), "empty_cells_by_column": {column: fields["empty"] for column, fields in columns.items()}, "field_statistics": columns, "sample_rows": sample, "sample_character_limit_per_cell": 180}


def mutation_separator_report(annotations: pd.Series, mutation_type, alphabet, requested_separator, sample_rows: int) -> dict:
    if requested_separator == "":
        raise ValueError("The mutation separator must be nonempty")

    separator_pattern = re.compile(r"([ \t]*[,;:/_|\\][ \t]*|[ \t]+)")
    evidence = {}
    substitution_positions = {}
    for annotation, count in annotations.value_counts(sort=False).items():
        parts = separator_pattern.split(annotation)
        if len(parts) < 3:
            continue

        substitutions = [part.strip() for part in parts[::2]]
        for substitution in substitutions:
            if substitution not in substitution_positions:
                try:
                    mutation = mutation_type.from_string(substitution, offset=0, alphabet=alphabet)
                except ValueError:
                    substitution_positions[substitution] = None
                else:
                    substitution_positions[substitution] = mutation.position if mutation.wild_type != mutation.mutant_type else None
        positions = [substitution_positions[substitution] for substitution in substitutions]
        if None in positions or len(set(positions)) != len(positions):
            continue

        for separator in dict.fromkeys(gap.strip() or ("\t" if "\t" in gap else " ") for gap in parts[1::2]):
            item = evidence.setdefault(separator, {"rows": 0, "examples": []})
            item["rows"] += int(count)
            if len(item["examples"]) < sample_rows:
                item["examples"].append(annotation)

    separator = requested_separator
    if requested_separator is not None:
        status = "explicit"
    elif len(evidence) == 1:
        status, separator = "detected", next(iter(evidence))
    elif evidence:
        status = "mixed"
    else:
        status, separator = "no_multi_substitution_evidence", ","
    return {"status": status, "used": separator, "evidence": evidence}


def inspect_candidates(table: pd.DataFrame, options) -> dict:
    tools = mutcleaner.cli.tables
    profile = table if options.profile_rows >= len(table) else table.head(options.profile_rows)
    columns = [options.reference_id_column, options.reference_sequence_column, options.variant_sequence_column, options.mutation_column, *options.label_column, *options.condition_column]
    tools.require_columns(profile, columns)
    empty = pd.Series("", index=profile.index, dtype="string")
    ids = profile[options.reference_id_column].str.strip() if options.reference_id_column else (empty if options.reference_id is None else pd.Series(options.reference_id.strip(), index=profile.index, dtype="string"))
    if options.reference_id_pattern is not None:
        if options.reference_id_column is None:
            raise ValueError("--reference-id-pattern requires --reference-id-column")
        ids = tools.extract_id(ids, options.reference_id_pattern)
    if options.sequence_type is None:
        if any((options.reference_sequence_column, options.reference_fasta, options.reference_table, options.variant_sequence_column, options.mutation_column, options.reference_slice)):
            raise ValueError("Specify --sequence-type to inspect sequences or substitution annotations")
        nonempty_ids = ids.loc[ids.ne("")]
        labels = {column: {"finite": int(mutcleaner.cli.measurements.finite_numbers(profile[column]).notna().sum())} for column in options.label_column}
        return {"profiled_rows": len(profile), "reference_id_nonempty": len(nonempty_ids), "distinct_reference_ids": int(nonempty_ids.nunique()), "labels": labels}
    variants = empty if options.variant_sequence_column is None else profile[options.variant_sequence_column].str.strip().str.upper()
    annotations = empty if options.mutation_column is None else profile[options.mutation_column].str.strip().str.upper()
    marker = None
    marker_mode = "none"
    if options.mutation_column and not options.no_wild_type_marker:
        marker = "WT" if options.wild_type_marker is None else options.wild_type_marker
        marker_mode = "candidate" if options.wild_type_marker is None else "explicit"
    wt = tools.wild_type_mask(annotations, marker)
    observed_references = empty if options.reference_sequence_column is None else profile[options.reference_sequence_column].str.strip().str.upper()
    if options.reference_fasta is not None:
        source = pd.DataFrame(read_fasta_records(options.reference_fasta), columns=["reference_id", "reference_sequence"], dtype="string")
    elif options.reference_table is not None:
        source = tools.read_reference_table(options.reference_table, options.reference_table_id_column, options.reference_table_sequence_column)
    elif options.reference_from_wt:
        source = pd.DataFrame({"reference_id": ids[wt], "reference_sequence": variants[wt]})
    else:
        source = pd.DataFrame({"reference_id": ids, "reference_sequence": observed_references})
    if options.reference_file_id_pattern is not None:
        source["reference_id"] = tools.extract_id(source["reference_id"], options.reference_file_id_pattern)
    source = source.assign(reference_id=source["reference_id"].str.strip(), reference_sequence=source["reference_sequence"].str.strip().str.upper())
    source = source.loc[source["reference_id"].ne("") & source["reference_sequence"].ne("")].drop_duplicates()
    conflict_ids = source.loc[source["reference_id"].duplicated(), "reference_id"].unique()
    unique = source.loc[~source["reference_id"].isin(conflict_ids)].set_index("reference_id")["reference_sequence"]
    references = ids.map(unique).fillna("")
    supplied_reference = options.reference_sequence_column is not None or options.reference_fasta is not None or options.reference_table is not None or options.reference_from_wt
    if options.reference_sequence_column is not None and options.reference_id_column is None and options.reference_id is None:
        references = observed_references

    out_of_slice = pd.Series(False, index=profile.index, dtype=bool)
    if options.reference_slice is not None:
        if not supplied_reference:
            raise ValueError("--reference-slice requires a supplied reference")
        start, stop = mutcleaner.cli.clean.reference_slice_bounds(options.reference_slice, options.sequence_type)
        out_of_slice = references.ne("") & references.str.len().lt(stop)
        references = references.str.slice(start, stop).mask(out_of_slice, "")
        if options.reference_from_wt:
            variants = variants.mask(wt, variants.str.slice(start, stop))

    sequence_class, mutation_type, alphabet = mutcleaner.cli.clean.sequence_settings(options.sequence_type, options.include_ambiguous, options.include_sec_pyl)
    default_alphabet = mutcleaner.cli.clean.sequence_settings(options.sequence_type, False, False)[2]
    ambiguous = set("BZJX" if options.sequence_type == "protein" else "RYSWKMBDHVN")
    special = set("UO") if options.sequence_type == "protein" else set()
    stop_symbols = {"*"} if options.sequence_type == "protein" else set()
    recognized = default_alphabet.letters | ambiguous | special | stop_symbols
    sequence_summary = {}
    for name, values in (("reference", references), ("variant", variants)):
        lengths = values.loc[values.ne("")].str.len()
        characters = Counter(character for value in values for character in value)
        sequence_summary[name] = {"nonempty": len(lengths), "length_min": int(lengths.min()) if len(lengths) else None, "length_max": int(lengths.max()) if len(lengths) else None, "whole_codon_lengths": int(lengths.mod(3).eq(0).sum()), "characters": dict(sorted(characters.items()))}
        sequence_summary[name]["symbol_counts"] = {
            "canonical": sum(characters[letter] for letter in default_alphabet.letters),
            "ambiguous": sum(characters[letter] for letter in ambiguous),
            "sec_pyl": sum(characters[letter] for letter in special),
            "stop": sum(characters[letter] for letter in stop_symbols),
            "invalid": sum(count for letter, count in characters.items() if letter not in recognized),
            "disallowed_by_selected_alphabet": sum(count for letter, count in characters.items() if letter not in alphabet.letters),
        }
    both_sequences = references.ne("") & variants.ne("")
    equal_sequences = both_sequences & references.eq(variants)
    labels = {column: {"finite": int(mutcleaner.cli.measurements.finite_numbers(profile[column]).notna().sum()), "marked_wild_type_finite": int(mutcleaner.cli.measurements.finite_numbers(profile.loc[wt, column]).notna().sum()), "equal_sequence_finite": int(mutcleaner.cli.measurements.finite_numbers(profile.loc[equal_sequences, column]).notna().sum()), "wild_type_or_equal_finite": int(mutcleaner.cli.measurements.finite_numbers(profile.loc[wt | equal_sequences, column]).notna().sum())} for column in options.label_column}
    result = {"profiled_rows": len(profile), "reference_id_nonempty": int(ids.ne("").sum()), "distinct_reference_ids": int(ids.loc[ids.ne("")].nunique()), "reference_matches": int(references.ne("").sum()), "reference_conflicts": len(conflict_ids), "conflicting_id_examples": list(conflict_ids[:10]), "sequences": sequence_summary, "sequence_pairs": int(both_sequences.sum()), "equal_length_pairs": int((both_sequences & references.str.len().eq(variants.str.len())).sum()), "equal_sequences": int(equal_sequences.sum()), "labels": labels}
    result["references_shorter_than_slice"] = int(out_of_slice.sum())

    if options.mutation_column is None:
        usable = both_sequences & ~equal_sequences
        names = pd.Series("", index=profile.index, dtype="string")
        if usable.any():
            pairs = pd.DataFrame({"reference": references[usable], "variant": variants[usable]})
            found, failures = mutcleaner.basic_cleaners.infer_mutations_from_sequences(pairs, wt_sequence_column="reference", mut_sequence_column="variant", sequence_type=options.sequence_type, alphabet=alphabet, num_workers=options.workers)
            names.loc[found.index] = found["inferred_mutations"]
            result["pair_inference_failures"] = len(failures["failed"])
        result["repeats"] = describe_repeats(profile, ids, names, options.label_column, options.condition_column)
        return result

    has_annotation = annotations.ne("") & ~wt
    separators = {separator: int(annotations.loc[has_annotation].str.contains(separator, regex=False).sum()) for separator in (",", ";", ":", "/", " ")}
    result["annotations"] = {"nonempty": int(annotations.ne("").sum()), "marked_wild_type": int(wt.sum()), "suspected_non_substitutions": int(annotations.loc[has_annotation].str.contains(r"DEL|INS|DUP|FS").sum()), "special_symbol_annotations": int(annotations.loc[has_annotation].str.contains(r"[BJOUXZ*]" if options.sequence_type == "protein" else r"[^ACGTU0-9\s,;:/]").sum()), "separator_evidence": separators, "samples": annotations.loc[has_annotation].drop_duplicates().head(options.sample_rows).tolist()}
    result["annotations"]["wild_type_marker"] = {"mode": marker_mode, "value": marker.strip() if marker is not None else None}
    result["offset_candidates"] = []
    separator_report = mutation_separator_report(annotations.loc[has_annotation], mutation_type, alphabet, options.mutation_separator, options.sample_rows)
    result["annotations"]["mutation_separator"] = separator_report
    mutation_separator = separator_report["used"]
    if mutation_separator is None:
        return result

    raw_positions = []
    for annotation in annotations.loc[has_annotation].unique():
        raw_positions.extend(int(position) for position in re.findall(r"\d+", annotation))
    result["annotations"].update(position_min=min(raw_positions) if raw_positions else None, position_max=max(raw_positions) if raw_positions else None)
    unit = 1 if options.sequence_type == "protein" else 3
    reference_offsets, position_maps = mutcleaner.cli.clean.read_position_rules(options)
    result["position_rules"] = {"offset_for": reference_offsets, "position_map": str(options.position_map.resolve()) if options.position_map is not None else None}
    candidate_offsets = options.offset or ((None,) if options.offset_for or options.position_map is not None else (0, 1))
    initial_counts = dict.fromkeys(("annotated_rows", "parsed_rows", "parse_failed_rows", "single_rows", "multi_rows", "reference_rows_checked", "in_range_rows", "out_of_range_rows", "reference_match_rows", "variant_rows_checked", "variant_match_rows", "reference_inference_failed_rows"), 0)

    @lru_cache(maxsize=8192)
    def parse_annotation(annotation, offset, mapped_reference_id):
        try:
            if mapped_reference_id is not None:
                annotation = mutcleaner.cli.clean.remap_annotation(annotation, mutation_separator, offset, position_maps[mapped_reference_id])
                offset = 0
            if offset is None:
                return None, "", False
            mutations = mutcleaner.mutation.MutationSet.from_string(annotation, sep=mutation_separator, offset=offset, mutation_type=mutation_type, alphabet=alphabet)
        except (TypeError, ValueError):
            below_range = mapped_reference_id is None and offset is not None and any(int(position) < offset for position in re.findall(r"\d+", annotation))
            return None, "", below_range
        return mutations, str(mutations), False

    @lru_cache(maxsize=8192)
    def rebuild_sequence(sequence, mutations, invert):
        if not invert and any(mutation.position < 0 or (mutation.position + 1) * unit > len(sequence) for mutation in mutations):
            return False, None
        try:
            if invert:
                mutations = mutcleaner.utils.mutation_converter.invert_mutation_set(mutations)
            return True, str(sequence_class(sequence, alphabet=alphabet).apply_mutation(mutations))
        except (TypeError, ValueError, IndexError):
            return True, None

    annotated_rows = pd.DataFrame({"reference_id": ids, "annotation": annotations, "reference": references, "variant": variants}).loc[has_annotation]
    for offset in candidate_offsets:
        per_reference = {}
        normalized = pd.Series("", index=profile.index, dtype="string")
        inferred_references = {}
        for index, reference_id, annotation, reference, variant in annotated_rows.itertuples(name=None):
            counts = per_reference.setdefault(reference_id, Counter(initial_counts))
            counts["annotated_rows"] += 1
            effective_offset = reference_offsets.get(reference_id, offset)
            mapped_reference_id = reference_id if reference_id in position_maps else None
            mutations, normalized_annotation, below_range = parse_annotation(annotation, effective_offset, mapped_reference_id)
            if mutations is None:
                counts["parse_failed_rows"] += 1
                if reference and below_range:
                    counts["out_of_range_rows"] += 1
                continue
            counts["parsed_rows"] += 1
            counts["single_rows" if len(mutations) == 1 else "multi_rows"] += 1
            normalized.at[index] = normalized_annotation
            if reference:
                counts["reference_rows_checked"] += 1
                counts["variant_rows_checked"] += bool(variant)
                in_range, rebuilt = rebuild_sequence(reference, mutations, False)
                if not in_range:
                    counts["out_of_range_rows"] += 1
                    continue
                counts["in_range_rows"] += 1
                if rebuilt is None:
                    continue
                counts["reference_match_rows"] += 1
                counts["variant_match_rows"] += bool(variant and rebuilt == variant)
            elif variant and reference_id and not supplied_reference:
                _, inferred = rebuild_sequence(variant, mutations, True)
                if inferred is None:
                    counts["reference_inference_failed_rows"] += 1
                else:
                    inferred_references.setdefault(reference_id, set()).add(inferred)
        totals = Counter(initial_counts)
        for counts in per_reference.values():
            totals.update(counts)
        detailed_references = [key for key in per_reference if key]
        result["offset_candidates"].append(
            {
                "offset": offset,
                **dict(totals),
                "inferred_reference_groups": len(inferred_references),
                "conflicting_inferred_reference_groups": sum(len(values) > 1 for values in inferred_references.values()),
                "reference_details": [{"reference_id": key, **per_reference[key], "conflicting_inferred_reference": len(inferred_references.get(key, ())) > 1} for key in detailed_references[:50]],
                "reference_details_truncated": len(detailed_references) > 50,
                "reference_groups": len(detailed_references),
                "repeats": describe_repeats(profile, ids, normalized, options.label_column, options.condition_column),
            }
        )
    return result


def describe_repeats(profile, ids, annotations, labels, conditions) -> dict:
    selected = ids.ne("") & annotations.ne("")
    details = profile.loc[selected].assign(_audit_reference=ids[selected], _audit_annotation=annotations[selected])
    groups = details.groupby(["_audit_reference", "_audit_annotation"], sort=False)
    size = groups.size()
    return {"groups_with_repeats": int(size.gt(1).sum()), "rows_beyond_first": int(size.sub(1).sum()), "groups_with_different_measurements": {column: int(groups[column].nunique().gt(1).sum()) for column in labels}, "groups_with_different_conditions": {column: int(groups[column].nunique().gt(1).sum()) for column in conditions}}
