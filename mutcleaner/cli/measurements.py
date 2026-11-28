import pandas as pd
import mutcleaner.cli.tables

from pathlib import Path


def measurement_settings(table: pd.DataFrame, arguments, mutation_column: str | None) -> dict:
    tools = mutcleaner.cli.tables
    labels = list(dict.fromkeys(arguments.label_column))
    targets = {column: float(value) for column, value in tools.parse_assignments(arguments.nearest_by, "--nearest-by").items()}
    weights = {column: float(value) for column, value in tools.parse_assignments(arguments.nearest_weight, "--nearest-weight").items()}
    group_columns = list(dict.fromkeys(arguments.group_by + arguments.split_by))
    baseline_columns = list(dict.fromkeys(arguments.wt_group_by + tuple(group_columns) + (tuple(targets) if arguments.subtract_wt else ())))
    if arguments.aggregate == "nearest" and not targets or arguments.aggregate != "nearest" and (targets or weights):
        raise ValueError("--aggregate nearest requires --nearest-by; nearest options apply only to nearest")
    if set(weights) - set(targets) or any(not 0 < value < float("inf") for value in weights.values()) or any(not float("-inf") < value < float("inf") for value in targets.values()):
        raise ValueError("Nearest targets must be finite and weights must be positive, finite target-column weights")
    if set(targets).intersection(group_columns):
        raise ValueError("Nearest target columns must not also be grouping or split columns")
    if arguments.aggregate == "first" and not group_columns and not arguments.assume_comparable:
        raise ValueError("Choose --group-by or --split-by, or use --assume-comparable for a single comparable experiment")
    if not arguments.subtract_wt and (arguments.wt_group_by or arguments.wt_baseline == "mean" or arguments.wt_from_equal):
        raise ValueError("WT baseline options require --subtract-wt")
    if arguments.subtract_wt and mutation_column is None and not arguments.wt_from_equal:
        raise ValueError("--subtract-wt requires a mutation column with measured WT rows or --wt-from-equal")
    if arguments.wt_from_equal and mutation_column is not None:
        raise ValueError("--wt-from-equal is for inputs with a reference and a variant sequence but no annotation column")

    condition_names = list(dict.fromkeys(group_columns + baseline_columns + list(targets)))
    tools.require_columns(table, labels + condition_names)
    conditions = {column: f"_condition_{index}" for index, column in enumerate(condition_names)}
    return {"labels": {column: f"_measurement_{index}" for index, column in enumerate(labels)}, "conditions": conditions, "group_columns": group_columns, "baseline_columns": baseline_columns, "targets": targets, "weights": weights}


def finite_numbers(values: pd.Series) -> pd.Series:
    numbers = pd.to_numeric(values, errors="coerce").astype("float64")
    return numbers.where(numbers.between(float("-inf"), float("inf"), inclusive="neither"))


def measurement_tasks(accepted: pd.DataFrame, unchanged: pd.DataFrame, arguments, settings: dict) -> tuple[list, list, list[str]]:
    aliases = settings["conditions"]
    group_keys = ["reference_id", "mutation_name", *(aliases[column] for column in settings["group_columns"])]
    baseline_keys = ["reference_id", *(aliases[column] for column in settings["baseline_columns"])]
    split_columns = list(dict.fromkeys(arguments.split_by))
    split_keys = [aliases[column] for column in split_columns]
    unsplit_keys = [aliases[column] for column in dict.fromkeys(arguments.group_by) if column not in split_columns]
    variant_keys = ["reference_id", "mutation_name"]
    groups = accepted.groupby(split_keys, sort=False, dropna=False) if split_keys else [((), accepted)]
    tasks = []
    summaries = []
    measurement_steps = dict.fromkeys(("split_conditions", "select_labels", "subtract_wild_type", "aggregate_measurements"), False)

    for group_number, (group_values, group) in enumerate(groups, start=1):
        condition_values = group_values if isinstance(group_values, tuple) else (group_values,)
        conditions = dict(zip(split_columns, condition_values, strict=True))
        if split_columns:
            measurement_steps["split_conditions"] = True
        task_wild_type = unchanged
        for column, value in conditions.items():
            task_wild_type = task_wild_type.loc[task_wild_type[aliases[column]].eq(value)]
        for label_number, (label_column, label_key) in enumerate(settings["labels"].items(), start=1):
            data = group.assign(label=finite_numbers(group[label_key]))
            measurement_steps["select_labels"] = True
            invalid_labels = data.loc[data["label"].isna()].groupby("reference_id").size().to_dict()
            data = data.loc[data["label"].notna()]
            baselines = task_wild_type.assign(label=finite_numbers(task_wild_type[label_key])).dropna(subset=["label"])
            wild_type_counts = baselines.groupby("reference_id").size().to_dict()
            missing_baseline = {}
            if arguments.subtract_wt and not data.empty:
                measured = baselines.groupby(baseline_keys, sort=False, dropna=False)["label"]
                baseline_counts = measured.size()
                baseline = measured.mean() if arguments.wt_baseline == "mean" else measured.first().where(baseline_counts.eq(1))
                baseline_table = pd.DataFrame({"_baseline": baseline, "_baseline_count": baseline_counts}).reset_index()
                data = data.merge(baseline_table, on=baseline_keys, how="left", validate="many_to_one", sort=False)
                if arguments.wt_baseline is None and data["_baseline_count"].gt(1).any():
                    raise ValueError(f"Choose --wt-baseline unique or mean for multiple matched WT measurements in label {label_column!r}")
                data["label"] = finite_numbers(data["label"] - data["_baseline"])
                measurement_steps["subtract_wild_type"] = True
                missing_baseline = data.loc[data["label"].isna()].groupby("reference_id").size().to_dict()
                data = data.loc[data["label"].notna()]

            missing_targets = {}
            nearest = [(aliases[column], target) for column, target in settings["targets"].items()]
            if nearest and not data.empty:
                data = data.copy()
                for column, _ in nearest:
                    data[column] = finite_numbers(data[column])
                missing_target = data[[column for column, _ in nearest]].isna().any(axis=1)
                missing_targets = data.loc[missing_target].groupby("reference_id").size().to_dict()
            before_aggregation = data.groupby("reference_id").size().to_dict()
            repeats = data.duplicated(subset=group_keys).any()
            if repeats and arguments.aggregate is None:
                raise ValueError(f"Choose --aggregate keep, mean, first or nearest for repeated measurements in label {label_column!r}")
            if repeats and arguments.aggregate == "mean" and not settings["group_columns"] and not arguments.assume_comparable:
                raise ValueError("Choose --group-by or --split-by, use --assume-comparable to average comparable repeats, or use --aggregate keep to retain all measurements")
            should_aggregate = arguments.aggregate in {"first", "nearest"} and not data.empty or arguments.aggregate == "mean" and repeats
            if should_aggregate:
                data = mutcleaner.basic_cleaners.aggregate_labels_by_name(data, group_keys, ["label"], strategy=arguments.aggregate, nearest_by=nearest or None, nearest_weights={aliases[column]: weight for column, weight in settings["weights"].items()} or None)
                measurement_steps["aggregate_measurements"] = True
                if finite_numbers(data["label"]).isna().any():
                    raise ValueError(f"Aggregation produced a non-finite value for label {label_column!r}")

            if unsplit_keys:
                combinations = data[variant_keys + unsplit_keys].drop_duplicates()
                conflicts = combinations.loc[combinations.duplicated(variant_keys, keep=False), variant_keys].drop_duplicates()
                if not conflicts.empty:
                    examples = ", ".join(f"{reference}: {mutation}" for reference, mutation in conflicts.head(5).itertuples(index=False, name=None))
                    raise ValueError(f"Different --group-by conditions would be indistinguishable in the output for {examples}. Use --split-by or filter to the desired condition.")

            task_path = Path(".")
            if len(settings["labels"]) > 1:
                task_path /= f"label_{label_number:03d}"
            if split_columns:
                task_path /= f"condition_{group_number:03d}"
            summary = {"label_column": label_column, "conditions": conditions, "path": str(task_path), "candidate_rows": len(group), "finite_wild_type_measurements": len(baselines), "invalid_label_rows": int(sum(invalid_labels.values())), "missing_wild_type_baseline_rows": int(sum(missing_baseline.values())), "missing_nearest_condition_rows": int(sum(missing_targets.values())), "rows_before_aggregation": int(sum(before_aggregation.values())), "rows_for_export": len(data)}
            summaries.append(summary)
            if not data.empty:
                tasks.append({**summary, "path": task_path, "rows": data, "before_aggregation": before_aggregation, "invalid_labels": invalid_labels, "missing_baseline": missing_baseline, "missing_targets": missing_targets, "finite_wt_by_reference": wild_type_counts, "condition_aliases": {aliases[column]: value for column, value in conditions.items()}, "aggregate": arguments.aggregate})
    return tasks, summaries, [step for step, executed in measurement_steps.items() if executed]
