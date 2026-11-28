import json
import math
import tempfile
import pandas as pd
import mutcleaner

from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

CLEANER_HEADER = '"""MutCleaner generated cleaner ({save_format}).'


def summarize_rows(rows, failures, unchanged, filtered, accepted) -> dict:
    return {
        "input_rows": len(rows),
        "valid_variant_rows": len(accepted),
        "invalid_rows": len(failures),
        "invalid_reasons": {str(reason): int(count) for reason, count in failures["error_message"].value_counts().items()},
        "filtered_rows": len(filtered),
        "filter_reasons": {str(reason): int(count) for reason, count in filtered["filter_reason"].value_counts().items()},
        "wild_type_or_unchanged_rows": len(unchanged),
        "explicit_wild_type_rows": int(unchanged["is_marked_wild_type"].sum()),
        "equal_sequence_rows": int(len(unchanged) - unchanged["is_marked_wild_type"].sum()),
    }


def validate_native_roundtrip(expected, actual, save_format: str) -> None:
    expected_refs = {key: (type(sequence), str(sequence)) for key, sequence in expected.reference_sequences.items()}
    actual_refs = {key: (type(sequence), str(sequence)) for key, sequence in actual.reference_sequences.items()}
    if actual_refs != expected_refs or len(actual) != len(expected):
        raise RuntimeError(f"Native {save_format} output restored different references or a different record count")

    expected_records, actual_records = {}, {}
    for dataset, records in ((expected, expected_records), (actual, actual_records)):
        for index, mutations in enumerate(dataset.mutation_sets):
            key = (dataset.get_mutation_set_reference(index), str(mutations))
            try:
                label = float(dataset.get_mutation_set_label(index))
            except (TypeError, ValueError, OverflowError) as error:
                raise RuntimeError(f"Non-numeric label after {save_format} round-trip for {key}") from error
            if not math.isfinite(label):
                raise RuntimeError(f"Non-finite label after {save_format} round-trip for {key}")
            records.setdefault(key, []).append(label)

    if expected_records.keys() != actual_records.keys():
        raise RuntimeError(f"Native {save_format} output restored different reference/variant associations")
    tolerance = 2e-15 if save_format == "dataframe" else 0.0
    for key, labels in expected_records.items():
        restored = actual_records[key]
        if len(labels) != len(restored) or any(not math.isclose(left, right, rel_tol=tolerance, abs_tol=0.0) for left, right in zip(sorted(labels), sorted(restored))):
            raise RuntimeError(f"Native {save_format} output restored different labels or repeated measurements for {key}")


def write_cleaner(script_path: Path, run: dict, result: dict, command_options: dict, custom_source: bytes | None) -> None:
    dependency_lines = []
    option_lines = []
    for flag, occurrences in command_options.items():
        if flag == "--custom-cleaner":
            continue
        for values in occurrences:
            tokens = [repr(flag)]
            for value in values:
                if isinstance(value, Path):
                    name = f"{flag[2:].replace('-', '_')}_path"
                    dependency_lines.append(f"{name} = {ascii(str(value.absolute()))}")
                    tokens.append(name)
                else:
                    tokens.append(ascii(str(value)))
            option_lines.append(f"    {', '.join(tokens)},")

    parameters = run["effective_parameters"]
    statistics = result["reference_statistics"]
    rows_summary = ", ".join(
        f"{statistics[key]} {label}"
        for key, label in (
            ("input_rows", "input"),
            ("invalid_rows", "invalid"),
            ("filtered_rows", "filtered"),
            ("wild_type_or_unchanged_rows", "WT or unchanged"),
            ("task_invalid_label_rows", "invalid label"),
            ("merged_measurements", "merged"),
            ("unselected_measurements", "unselected"),
            ("output_rows", "exported"),
            ("export_failed_rows", "export failures"),
        )
    )
    notes = [
        f"MutCleaner version: {run['mutcleaner_version']}",
        f"Git commit: {run['git_commit'] or 'unavailable'}",
        f"Original input: {ascii(parameters['input_path'])}",
        f"Original output: {ascii(parameters['output_path'])}",
        f"Result path: {ascii(result['path'])}",
        f"Reference ID: {ascii(result['reference_id'])}",
        f"Label column: {ascii(result['label_column'])}",
        f"Conditions: {ascii(result['conditions'])}",
        f"Steps: {', '.join(run['steps'])}",
        f"Rows: {rows_summary}",
        f"Verification: {', '.join(f'{key}={value}' for key, value in result['verification'].items())}",
    ]
    custom_imports = []
    custom_definition = []
    replay = ['    subprocess.run(["mutcleaner", "clean", *cleaning_options, "--", args.input_path, args.output_path], check=True)']
    if custom_source is not None:
        custom = run["custom"]
        notes.extend(
            [
                f"Custom cleaner SHA-256: {custom['sha256']}",
                f"Custom cleaner: {custom['input_rows']} to {custom['output_rows']} rows, {custom['input_columns']} to {custom['output_columns']} columns",
            ]
        )
        custom_imports = ["import tempfile", "", "from pathlib import Path"]
        custom_definition = [f"custom_source = {custom_source!r}", ""]
        replay = [
            "    with tempfile.TemporaryDirectory() as temporary_directory:",
            '        custom_path = Path(temporary_directory) / "custom.py"',
            "        custom_path.write_bytes(custom_source)",
            '        subprocess.run(["mutcleaner", "clean", *cleaning_options, "--custom-cleaner", str(custom_path), "--", args.input_path, args.output_path], check=True)',
        ]
    lines = [
        CLEANER_HEADER.format(save_format=result["save_format"]),
        "",
        "Reproduce the complete original cleaning task, including all its outputs.",
        "Usage: python cleaner.py INPUT_PATH NEW_OUTPUT_PATH",
        "Use the same MutCleaner build on PATH. Source builds can differ from PyPI",
        "releases with the same version. Edit any dependency paths below as needed.",
        '"""',
        "",
        "import argparse",
        "import subprocess",
        *custom_imports,
        "",
        *(f"# {note}" for note in notes),
        "",
        *dependency_lines,
        *([""] if dependency_lines else []),
        *custom_definition,
        "cleaning_options = [",
        *option_lines,
        "]",
        "",
        'if __name__ == "__main__":',
        '    parser = argparse.ArgumentParser(description="Reproduce the complete MutCleaner cleaning task.")',
        '    parser.add_argument("input_path")',
        '    parser.add_argument("output_path", help="New output directory")',
        "    args = parser.parse_args()",
        *replay,
        "",
    ]
    with script_path.open("x", encoding="utf-8", newline="\n") as script:
        script.write("\n".join(lines))


def write_results(output_path: Path, tasks: list, run: dict, observations: tuple, sequence_options: dict, save_format: str, command_options: dict, custom_source: bytes | None) -> dict:
    if not tasks:
        return {**run, "status": "no_variants", "finished_at": datetime.now(timezone.utc).isoformat(), "outputs": []}

    #######################################################################
    # export and verify all tasks in a private directory
    #######################################################################

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_reports = []
    with tempfile.TemporaryDirectory(prefix=".mutcleaner-", dir=output_path.parent) as temporary_directory:
        working_path = Path(temporary_directory)
        staging_path = working_path / "results"
        verify_path = staging_path if save_format == "mutcleaner" else working_path / "standard_roundtrip"
        for task in tasks:
            task_path = verify_path / task["path"]
            dataset = task["dataset"]
            dataset.save_by_reference(task_path)
            reloaded = mutcleaner.MutationDataset.load_by_reference(task_path, offset=0, **sequence_options)
            validation = reloaded.validate_against_references()
            if validation["invalid_mutation_sets"] or validation["position_mismatches"]:
                raise RuntimeError(f"Saved dataset did not pass reference validation: {task['path']}")

            task_observations = observations
            for condition, value in task["condition_aliases"].items():
                task_observations = tuple(frame.loc[frame[condition].eq(value)] for frame in task_observations)
            observation_indices = tuple(frame.groupby("reference_id", sort=False).indices for frame in task_observations)
            prepared_rows = task["rows"]
            prepared_indices = prepared_rows.groupby("reference_id", sort=False).indices
            expected_columns = ["mutation_name", "variant_sequence", "label"]
            verified_rows = 0
            seen_references = set()
            references = reloaded.reference_sequences
            for reference_path in sorted(task_path.iterdir()):
                metadata = json.loads((reference_path / "metadata.json").read_text(encoding="utf-8"))
                reference_id = str(metadata["reference_id"])
                reference_rows = prepared_rows.iloc[prepared_indices.get(reference_id, [])]
                exported = pd.read_csv(reference_path / "data.csv", keep_default_na=False, dtype={"mutation_name": "string", "mutated_sequence": "string", "label": "float64"}, float_precision="round_trip")
                failed_path = reference_path / "failed_data.csv"
                export_failures = len(pd.read_csv(failed_path)) if failed_path.is_file() else 0
                expected = Counter(reference_rows[expected_columns].itertuples(index=False, name=None))
                actual = Counter(exported[["mutation_name", "mutated_sequence", "label"]].itertuples(index=False, name=None))
                if reference_id not in references or str(references[reference_id]) != str(dataset.reference_sequences[reference_id]):
                    raise RuntimeError(f"Saved reference differs from the input reference: {reference_id}")
                if actual - expected or len(exported) + export_failures != len(reference_rows) or len(exported) != metadata["num_mutation_sets"]:
                    raise RuntimeError(f"Saved sequences, labels, or row counts differ from the prepared dataset: {reference_id}")

                verified_rows += len(exported)
                seen_references.add(reference_id)
                statistics = summarize_rows(*(frame.iloc[indices.get(reference_id, [])] for frame, indices in zip(task_observations, observation_indices, strict=True)))
                before_aggregation = int(task["before_aggregation"].get(reference_id, 0))
                invalid_task_rows = int(task["invalid_labels"].get(reference_id, 0))
                no_baseline_rows = int(task["missing_baseline"].get(reference_id, 0))
                statistics.update(
                    {
                        "task_rows_before_aggregation": before_aggregation,
                        "task_invalid_label_rows": invalid_task_rows,
                        "task_missing_wild_type_baseline_rows": no_baseline_rows,
                        "task_missing_nearest_condition_rows": int(task["missing_targets"].get(reference_id, 0)),
                        "task_finite_wild_type_measurements": int(task["finite_wt_by_reference"].get(reference_id, 0)),
                        "merged_measurements": before_aggregation - len(reference_rows) if task["aggregate"] == "mean" else 0,
                        "unselected_measurements": before_aggregation - len(reference_rows) if task["aggregate"] in {"first", "nearest"} else 0,
                        "output_rows": len(exported),
                        "export_failed_rows": export_failures,
                        "mutations": metadata["total_mutations"],
                        "covered_positions": metadata["covered_positions"],
                        "coverage_percentage": metadata["coverage_percentage"],
                    }
                )
                relative_path = reference_path.relative_to(verify_path)
                if save_format != "mutcleaner":
                    if export_failures:
                        raise RuntimeError(f"Cannot publish a partial {save_format} export for reference {reference_id}")
                    destination_path = staging_path / relative_path
                    destination_path.mkdir(parents=True, exist_ok=False)
                    native_path = destination_path / "dataset"
                    reference_dataset = dataset.filter_by_reference(reference_id)
                    reference_dataset.save(str(native_path), save_type=save_format)
                    native = mutcleaner.MutationDataset.load(str(native_path), load_type=save_format)
                    validate_native_roundtrip(reference_dataset, native, save_format)
                output_reports.append(
                    {
                        "path": str(relative_path),
                        "save_format": save_format,
                        "reference_id": reference_id,
                        "label_column": task["label_column"],
                        "conditions": task["conditions"],
                        "reference_statistics": statistics,
                        "verification": {"reference": "passed", "variant_sequences": "passed", "labels": "passed", "readback": "passed", "export": "partial" if export_failures else "passed"},
                        "status": "partial" if export_failures else "success",
                    }
                )

            if seen_references != set(dataset.reference_sequences) or len(reloaded) != verified_rows:
                raise RuntimeError(f"Reloaded dataset is incomplete: {task['path']}")

        status = "partial" if any(report["status"] == "partial" for report in output_reports) else "success"
        finished_at = datetime.now(timezone.utc).isoformat()
        export_steps = ["save_and_verify_standard"] if save_format == "mutcleaner" else ["verify_standard_intermediate", f"save_and_verify_{save_format}"]
        run = {**run, "steps": [*run["steps"], *export_steps]}
        for result in output_reports:
            write_cleaner(staging_path / result["path"] / "cleaner.py", run, result, command_options, custom_source)

        #######################################################################
        # claim the output path and publish verified results
        #######################################################################

        output_path.mkdir(exist_ok=False)
        for entry in staging_path.iterdir():
            entry.rename(output_path / entry.name)

    return {**run, "status": status, "finished_at": finished_at, "outputs": output_reports}
