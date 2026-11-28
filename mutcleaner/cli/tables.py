import hashlib
import json
import re
import pandas as pd

from pathlib import Path

TABLE_FORMATS = ("csv", "tsv", "txt", "xlsx", "xls", "parquet")
FASTA_FORMATS = ("fa", "faa", "fna", "fas", "fasta")


def read_table(input_path: Path, file_format, delimiter, encoding, sheet, list_sheets: bool) -> tuple[pd.DataFrame | None, dict]:
    file_format = file_format or input_path.suffix.lower().lstrip(".")
    if file_format not in TABLE_FORMATS:
        raise ValueError(f"Unsupported table format: {file_format or '(none)'}")
    settings = {"format": file_format}
    if file_format in {"csv", "tsv", "txt"}:
        if sheet is not None:
            raise ValueError("--sheet applies only to Excel input")
        delimiter = delimiter if delimiter is not None else ("," if file_format == "csv" else "\t")
        if delimiter == r"\t":
            delimiter = "\t"
        if len(delimiter) != 1:
            raise ValueError("--delimiter must be a single character")
        settings.update(delimiter=delimiter, encoding=encoding or "utf-8-sig")
        table = pd.read_csv(input_path, sep=delimiter, encoding=settings["encoding"], dtype="string", keep_default_na=False)
    elif file_format in {"xls", "xlsx"}:
        if delimiter is not None or encoding is not None:
            raise ValueError("--delimiter and --encoding apply only to text tables")
        with pd.ExcelFile(input_path, engine="xlrd" if file_format == "xls" else "openpyxl") as workbook:
            settings["sheets"] = workbook.sheet_names
            if sheet is None and len(workbook.sheet_names) > 1:
                if list_sheets:
                    return None, settings
                raise ValueError(f"Select --sheet or --sheet-index for a workbook with multiple sheets: {', '.join(workbook.sheet_names)}")
            settings["sheet"] = workbook.sheet_names[0] if sheet is None else sheet
            table = workbook.parse(sheet_name=settings["sheet"], dtype="string", keep_default_na=False)
    else:
        if any(option is not None for option in (delimiter, encoding, sheet)):
            raise ValueError("Parquet input does not accept --delimiter, --encoding or --sheet")
        table = pd.read_parquet(input_path, engine="pyarrow").astype("string")

    table.columns = table.columns.map(str)
    return table.fillna("").reset_index(drop=True), settings


def apply_custom_cleaner(table: pd.DataFrame, script_path: Path) -> tuple[pd.DataFrame, dict, bytes]:
    original_rows = len(table)
    original_columns = table.columns.tolist()
    script_path = script_path.resolve()
    try:
        source = script_path.read_bytes()
        namespace = {"__name__": "__mutcleaner_custom__", "__file__": str(script_path)}
        exec(compile(source, str(script_path), "exec"), namespace)
        cleaner = namespace.get("custom_cleaner")
        if not callable(cleaner):
            raise ValueError("The script must define a callable custom_cleaner(df)")
        result = cleaner(table)
        if not isinstance(result, pd.DataFrame):
            raise ValueError("custom_cleaner(df) must return a pandas DataFrame")
        if not result.columns.is_unique:
            raise ValueError("Custom cleaner output column names must be unique")
        columns = ["" if pd.api.types.is_scalar(column) and pd.isna(column) else str(column) for column in result.columns]
        if any(not column.strip() for column in columns) or len(set(columns)) != len(columns):
            raise ValueError("Custom cleaner output column names must be nonempty and unique after conversion to strings")
        result = result.set_axis(columns, axis="columns").astype("string").fillna("").reset_index(drop=True)
    except (Exception, SystemExit) as error:
        raise RuntimeError(f"Custom cleaner {script_path}: {error}") from error

    statistics = {
        "path": str(script_path),
        "sha256": hashlib.sha256(source).hexdigest(),
        "input_rows": original_rows,
        "input_columns": len(original_columns),
        "output_rows": len(result),
        "output_columns": len(columns),
        "added_columns": [column for column in columns if column not in original_columns],
        "removed_columns": [column for column in original_columns if column not in columns],
        "source_row_basis": "One-based data row in the processed table, after the custom cleaner and field preparation.",
    }
    return result, statistics, source


def require_columns(table: pd.DataFrame, columns) -> None:
    missing = list(dict.fromkeys(column for column in columns if column is not None and column not in table.columns))
    if missing:
        raise ValueError(f"Input is missing columns: {', '.join(missing)}")


def parse_assignments(options, option_name) -> dict:
    parsed = {}
    for option in options:
        key, separator, value = option.partition("=")
        if not separator or not key.strip() or key.strip() in parsed:
            raise ValueError(f"{option_name} requires distinct NAME=VALUE entries: {option!r}")
        parsed[key.strip()] = value
    return parsed


def preprocess(table: pd.DataFrame, options) -> tuple[pd.DataFrame, str | None]:
    table = table.copy()
    mutation_column = options.mutation_column
    rename = parse_assignments(options.rename, "--rename")
    require_columns(table, rename)
    table = table.rename(columns=rename)
    if table.columns.duplicated().any() or any(not name for name in rename.values()):
        raise ValueError("--rename must produce distinct, nonempty column names")

    constants = parse_assignments(options.constant, "--constant")
    for column, value in constants.items():
        if column in table:
            raise ValueError(f"--constant would replace an existing column: {column}")
        table[column] = value

    for target, sources, delimiter in options.join_column:
        columns = sources.split(",")
        require_columns(table, columns)
        if target in table or not target:
            raise ValueError(f"--join-column requires a new column: {target!r}")
        table[target] = table[columns[0]].str.cat([table[column] for column in columns[1:]], sep=delimiter) if len(columns) > 1 else table[columns[0]]

    for source, separator, targets in options.split_column:
        require_columns(table, [source])
        columns = targets.split(",")
        if not separator or any(not column or column in table for column in columns) or len(set(columns)) != len(columns):
            raise ValueError("--split-column requires a separator and distinct new columns")
        values = table[source].str.split(separator, regex=False)
        if values.str.len().gt(len(columns)).any():
            raise ValueError(f"Column {source!r} contains more parts than --split-column names")
        for index, column in enumerate(columns):
            table[column] = values.str.get(index).fillna("")

    for column, pattern, replacement in options.replace:
        require_columns(table, [column])
        table[column] = table[column].str.replace(pattern, replacement, regex=True)

    for target, source, pattern in options.extract:
        require_columns(table, [source])
        if target in table or not target:
            raise ValueError(f"--extract requires a new column: {target!r}")
        table[target] = extract_id(table[source], pattern)

    if options.mutation_parts is not None:
        if mutation_column is not None:
            raise ValueError("Use either --mutation-column or --mutation-parts")
        wt_column, position_column, mutant_column = options.mutation_parts
        require_columns(table, options.mutation_parts)
        mutation_column = "__mutcleaner_input_mutation"
        if mutation_column in table:
            raise ValueError(f"Reserved input column name: {mutation_column}")
        table[mutation_column] = table[wt_column].str.strip() + table[position_column].str.strip() + table[mutant_column].str.strip()
    return table, mutation_column


def extract_id(values: pd.Series, pattern: str) -> pd.Series:
    try:
        compiled = re.compile(pattern)
    except re.PatternError as error:
        raise ValueError(str(error)) from error
    if compiled.groups != 1:
        raise ValueError("ID or field extraction requires exactly one capturing group")
    return values.str.extract(compiled, expand=False).fillna("").str.strip()


def filter_reasons(table: pd.DataFrame, options) -> pd.Series:
    reasons = pd.Series("", index=table.index, dtype="string")
    relations = {"<": pd.Series.lt, "<=": pd.Series.le, ">": pd.Series.gt, ">=": pd.Series.ge}
    for condition in options.filter:
        match = re.fullmatch(r"(.+?)\s*(>=|<=|!=|==|=|>|<)(.*)", condition)
        if match is None:
            raise ValueError(f"Invalid --filter expression: {condition!r}")
        column, operator, target = match.groups()
        column = column.strip()
        require_columns(table, [column])
        cells = table[column].str.strip()
        if operator in relations:
            number = float(target)
            if not float("-inf") < number < float("inf"):
                raise ValueError(f"Non-finite filter target: {condition!r}")
            numeric = pd.to_numeric(cells, errors="coerce")
            keep = relations[operator](numeric, number) & numeric.between(float("-inf"), float("inf"), inclusive="neither")
        else:
            keep = cells.eq(target) if operator in {"=", "=="} else cells.ne(target)
        reasons.loc[reasons.eq("") & ~keep.fillna(False)] = f"Filter excluded: {condition}"

    for column, pattern in parse_assignments(options.exclude, "--exclude").items():
        require_columns(table, [column])
        reasons.loc[reasons.eq("") & table[column].str.contains(pattern, regex=True, na=False)] = f"Pattern excluded: {column}={pattern}"
    for column in options.drop_missing:
        require_columns(table, [column])
        reasons.loc[reasons.eq("") & table[column].str.strip().eq("")] = f"Missing required value: {column}"
    return reasons


def read_reference_table(file_path: Path, id_column, sequence_column) -> pd.DataFrame:
    if file_path.suffix.lower() == ".json":
        if id_column is not None or sequence_column is not None:
            raise ValueError("JSON reference tables already map IDs directly to sequences")
        references = json.loads(file_path.read_text(encoding="utf-8"))
        if not isinstance(references, dict) or any(not isinstance(sequence, str) for sequence in references.values()):
            raise ValueError("Reference JSON must be an object mapping IDs to sequence strings")
        return pd.DataFrame(references.items(), columns=["reference_id", "reference_sequence"], dtype="string")

    if file_path.suffix.lower() not in {".csv", ".tsv"} or id_column is None or sequence_column is None:
        raise ValueError("CSV/TSV reference tables require --reference-table-id-column and --reference-table-sequence-column")
    table, _ = read_table(file_path, None, None, None, None, list_sheets=False)
    require_columns(table, [id_column, sequence_column])
    return pd.DataFrame({"reference_id": table[id_column], "reference_sequence": table[sequence_column]})


def wild_type_mask(annotations: pd.Series, marker: str | None) -> pd.Series:
    if marker is None:
        return pd.Series(False, index=annotations.index, dtype=bool)
    effective_marker = marker.strip()
    if not effective_marker:
        raise ValueError("The wild-type marker must be nonempty")
    return annotations.str.strip().str.casefold().eq(effective_marker.casefold())
