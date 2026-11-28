import argparse
import contextlib
import importlib.resources
import json
import re
import sys
import mutcleaner.cli.clean
import mutcleaner.cli.inspection
import mutcleaner.cli.tables

from itertools import islice
from pathlib import Path
from types import SimpleNamespace

CUSTOM_CLEANER_HELP = "Run custom_cleaner(df) on the full input table before processing. Only use trusted scripts; this executes arbitrary Python code."


class OptionValue(str):
    """Mark a dash-prefixed token already consumed as an option value."""


def unprotect_option_value(value):
    return value[1:] if isinstance(value, OptionValue) else value


def protect_option_values(argv, options):
    arguments = iter(argv)
    protected = []
    for token in arguments:
        if token == "--":
            protected.extend((token, *arguments))
            break
        flag, equals, attached = token.partition("=")
        if flag not in options:
            protected.append(token)
            continue
        _, nargs, _ = options[flag]
        count = 1 if nargs is None else nargs
        if count == 0 or (equals and count == 1):
            protected.append(token)
            continue
        values = ([attached] if equals else []) + list(islice(arguments, count - bool(equals)))
        protected.append(flag)
        protected.extend(OptionValue(f" {value}") if value.startswith("-") else value for value in values)
    return protected


def add_option(group, options, flag, **kwargs):
    action = group.add_argument(flag, **kwargs)
    options[flag] = (action.dest, action.nargs, kwargs.get("action") == "append")


def bounded_int(minimum, maximum):
    def convert(value):
        number = int(unprotect_option_value(value))
        if number < minimum or (maximum is not None and number > maximum):
            bounds = f"{minimum} to {maximum}" if maximum is not None else f"at least {minimum}"
            raise argparse.ArgumentTypeError(f"expected an integer {bounds}")
        return number

    return convert


def guide(topic) -> str:
    """Show the complete guide or a specific topic, including the agent workflow."""
    guide_path = importlib.resources.files("mutcleaner.cli").joinpath("guide.md")
    content = guide_path.read_text(encoding="utf-8")
    if topic is None:
        content = re.sub(r"^<!-- topic: [a-z]+ -->\n", "", content, flags=re.MULTILINE)
    else:
        _, _, section = content.partition(f"<!-- topic: {topic} -->\n")
        content = section.partition("\n<!-- topic: ")[0].rstrip() + "\n"
    return content


def inspect_table(input_path: Path, options) -> dict:
    """Report evidence from a table, FASTA, or directory without modifying it."""
    profiles = mutcleaner.cli.inspection
    is_fasta = options.file_format is None and input_path.suffix.lower().lstrip(".") in mutcleaner.cli.tables.FASTA_FORMATS
    if options.custom_cleaner and (input_path.is_dir() or options.verify or is_fasta):
        raise ValueError("--custom-cleaner requires a table and cannot be used with FASTA, directories, or --verify")
    if input_path.is_dir():
        report = profiles.describe_directory(input_path, options.verify)
    elif options.verify:
        raise ValueError("--verify requires a standard dataset directory")
    elif is_fasta:
        report = profiles.describe_fasta(input_path, options.sample_rows)
    else:
        if options.sheet is not None and options.sheet_index is not None:
            raise ValueError("Use either --sheet or --sheet-index")
        if options.reference_id is not None and options.reference_id_column is not None:
            raise ValueError("Use either --reference-id or --reference-id-column")
        if sum(source is not None for source in (options.reference_fasta, options.reference_table, options.reference_sequence_column)) + options.reference_from_wt > 1:
            raise ValueError("Inspect one reference source at a time")
        if options.reference_id_pattern is not None and options.reference_id_column is None:
            raise ValueError("--reference-id-pattern requires --reference-id-column")
        if options.reference_file_id_pattern is not None and options.reference_table is None and options.reference_fasta is None:
            raise ValueError("--reference-file-id-pattern requires --reference-fasta or --reference-table")
        if options.reference_table is None and (options.reference_table_id_column is not None or options.reference_table_sequence_column is not None):
            raise ValueError("Reference table columns require --reference-table")
        if (options.offset or options.offset_for or options.position_map is not None or options.wild_type_marker is not None or options.no_wild_type_marker) and options.mutation_column is None:
            raise ValueError("Numbering and WT marker options require --mutation-column for inspection")
        if options.wild_type_marker is not None and options.no_wild_type_marker:
            raise ValueError("Use either --wild-type-marker or --no-wild-type-marker")
        if options.reference_from_wt and (options.mutation_column is None or options.variant_sequence_column is None):
            raise ValueError("--reference-from-wt requires mutation and variant sequence columns")
        if options.reference_from_wt and options.no_wild_type_marker:
            raise ValueError("--reference-from-wt requires a WT marker")
        sheet = options.sheet_index if options.sheet_index is not None else options.sheet
        table, settings = mutcleaner.cli.tables.read_table(input_path, options.file_format, options.delimiter, options.encoding, sheet, list_sheets=True)
        if table is None:
            if options.custom_cleaner:
                raise ValueError("--custom-cleaner requires a selected worksheet; specify --sheet or --sheet-index")
            report = {"kind": "workbook", "sheets": settings["sheets"]}
        else:
            custom_statistics = None
            if options.custom_cleaner:
                table, custom_statistics, _ = mutcleaner.cli.tables.apply_custom_cleaner(table, options.custom_cleaner[0])
            if options.profile_rows is None:
                options.profile_rows = len(table)
            report = profiles.describe_table(table, options.sample_rows, options.profile_rows)
            if custom_statistics is not None:
                report["custom"] = custom_statistics
            if any((options.reference_id_column, options.reference_id, options.reference_sequence_column, options.reference_fasta, options.reference_table, options.variant_sequence_column, options.mutation_column, options.label_column, options.condition_column, options.reference_from_wt, options.reference_slice)):
                report["candidate_fields"] = profiles.inspect_candidates(table, options)
        report.update(input_path=str(input_path.resolve()), read_settings=settings)
    return report


def clean(input_path: Path, output_path: Path, options, option_specs) -> dict:
    """Clean measured variants and write verified datasets to a new OUTPUT_PATH."""
    sheet = options["sheet_index"] if options["sheet_index"] is not None else options["sheet"]
    file_format = options["file_format"] or input_path.suffix.lower().lstrip(".")
    if options["custom_cleaner"] and (input_path.is_dir() or file_format not in mutcleaner.cli.tables.TABLE_FORMATS):
        raise ValueError("--custom-cleaner requires a supported table input")
    table, settings = mutcleaner.cli.tables.read_table(input_path, options["file_format"], options["delimiter"], options["encoding"], sheet, list_sheets=False)
    replay_options = {
        **options,
        "file_format": settings["format"],
        "delimiter": settings.get("delimiter"),
        "encoding": settings.get("encoding"),
        "sheet": settings.get("sheet") if isinstance(settings.get("sheet"), str) else None,
        "sheet_index": settings.get("sheet") if isinstance(settings.get("sheet"), int) else None,
    }
    command_options = {}
    for flag, (destination, nargs, repeated) in option_specs.items():
        value = replay_options[destination]
        if value is None:
            continue
        if nargs == 0:
            if value:
                command_options[flag] = [()]
            continue
        occurrences = value if repeated else (value,)
        command_options[flag] = [(item,) if nargs is None else item for item in occurrences]
    arguments = SimpleNamespace(input_path=input_path, output_path=output_path, command_options=command_options, read_settings=settings, **options)
    return mutcleaner.cli.clean.clean_table(table, arguments)


def add_input_arguments(group, options):
    add_option(group, options, "--format", dest="file_format", choices=mutcleaner.cli.tables.TABLE_FORMATS)
    add_option(group, options, "--delimiter", help="Text table field delimiter, or \\t for tab.")
    add_option(group, options, "--encoding", help="Text table encoding.")
    add_option(group, options, "--sheet", help="Excel sheet name.")
    add_option(group, options, "--sheet-index", type=bounded_int(0, None), help="Zero-based Excel sheet index.")
    add_option(group, options, "--custom-cleaner", type=Path, action="append", help=CUSTOM_CLEANER_HELP)


def add_inspect_arguments(inspect_parser, options):
    inspect_parser.add_argument("input_path", type=Path, metavar="INPUT_PATH")
    add_input_arguments(inspect_parser, options)
    add_option(inspect_parser, options, "--sample-rows", type=bounded_int(0, 10), default=3, help="Number of sample rows to display (default: %(default)s).")
    add_option(inspect_parser, options, "--profile-rows", type=bounded_int(1, None), help="Maximum number of leading rows to profile. Defaults to all rows.")
    add_option(inspect_parser, options, "--reference-id-column")
    add_option(inspect_parser, options, "--reference-id")
    add_option(inspect_parser, options, "--reference-id-pattern")
    add_option(inspect_parser, options, "--reference-from-wt", action="store_true", help="Use measured WT rows as the reference sequence source.")
    add_option(inspect_parser, options, "--reference-sequence-column")
    add_option(inspect_parser, options, "--reference-fasta", type=Path)
    add_option(inspect_parser, options, "--reference-table", type=Path)
    add_option(inspect_parser, options, "--reference-table-id-column")
    add_option(inspect_parser, options, "--reference-table-sequence-column")
    add_option(inspect_parser, options, "--reference-file-id-pattern")
    add_option(inspect_parser, options, "--variant-sequence-column")
    add_option(inspect_parser, options, "--mutation-column")
    add_option(inspect_parser, options, "--mutation-separator", help="Separator between substitutions. Detected from the profiled annotations when omitted.")
    add_option(inspect_parser, options, "--wild-type-marker", help="Measured WT annotation. Matching ignores case and edge whitespace. If omitted, WT is checked as a candidate.")
    add_option(inspect_parser, options, "--no-wild-type-marker", action="store_true", help="Do not recognize any annotation as measured WT.")
    add_option(inspect_parser, options, "--sequence-type", choices=["protein", "dna", "rna"], help="Required when checking sequences or substitution annotations.")
    add_option(inspect_parser, options, "--workers", type=bounded_int(1, None), default=4, help="Workers for supported sequence operations (default: %(default)s).")
    add_option(inspect_parser, options, "--offset", type=int, action="append", help="Candidate annotation offset. Defaults to checking 0 and 1 without advanced rules.")
    add_option(inspect_parser, options, "--offset-for", action="append", metavar="REFERENCE=OFFSET", help="Use an offset for a specific reference.")
    add_option(inspect_parser, options, "--position-map", type=Path, help="Table with reference_id, input_position, reference_index.")
    add_option(inspect_parser, options, "--reference-slice", type=int, nargs=2, metavar=("START", "STOP"), help="Zero-based half-open reference slice in residues or codons.")
    add_option(inspect_parser, options, "--include-ambiguous", action="store_true", help="Allow B, Z, J and X in protein sequences.")
    add_option(inspect_parser, options, "--include-sec-pyl", action="store_true", help="Allow U and O in protein sequences.")
    add_option(inspect_parser, options, "--label-column", action="append", help="Candidate measurement column.")
    add_option(inspect_parser, options, "--condition-column", action="append", help="Candidate experiment condition column.")
    add_option(inspect_parser, options, "--verify", action="store_true", help="Read back and verify a standard output directory.")


def add_clean_arguments(clean_parser, options):
    clean_parser.add_argument("input_path", type=Path, metavar="INPUT_PATH")
    clean_parser.add_argument("output_path", type=Path, metavar="OUTPUT_PATH")
    add_input_arguments(clean_parser.add_argument_group("Input format"), options)

    group = clean_parser.add_argument_group("Prepare and filter fields")
    add_option(group, options, "--rename", action="append", metavar="OLD=NEW", help="Rename an input column.")
    add_option(group, options, "--constant", action="append", metavar="COLUMN=VALUE", help="Add a constant column.")
    add_option(group, options, "--join-column", nargs=3, action="append", metavar=("TARGET", "SOURCES", "SEPARATOR"), help="Join comma-separated source columns.")
    add_option(group, options, "--split-column", nargs=3, action="append", metavar=("SOURCE", "SEPARATOR", "TARGETS"), help="Split a field into comma-separated target names.")
    add_option(group, options, "--replace", nargs=3, action="append", metavar=("COLUMN", "PATTERN", "REPLACEMENT"), help="Apply a regex replacement in a column.")
    add_option(group, options, "--extract", nargs=3, action="append", metavar=("TARGET", "SOURCE", "PATTERN"), help="Extract a single regex capture group into a column.")
    add_option(group, options, "--filter", action="append", metavar="COLUMN=VALUE", help="Select rows with =, !=, <, <=, > or >= comparisons.")
    add_option(group, options, "--exclude", action="append", metavar="COLUMN=PATTERN", help="Exclude rows matching a regex.")
    add_option(group, options, "--drop-missing", action="append", metavar="COLUMN", help="Exclude rows with an empty value in a column.")

    group = clean_parser.add_argument_group("Reference sources and IDs")
    add_option(group, options, "--reference-id-column", help="Column identifying the reference for each row.")
    add_option(group, options, "--reference-id", help="Reference identifier for an input with a single reference.")
    add_option(group, options, "--reference-id-pattern", help="Extract the reference ID using one regex capture group.")
    add_option(group, options, "--reference-from-wt", action="store_true", help="Use measured WT rows as the reference sequence source.")
    add_option(group, options, "--infer-reference", action="store_true", help="Reconstruct references from substitution annotations and observed variant sequences.")
    add_option(group, options, "--reference-sequence-column", help="Column containing reference sequences.")
    add_option(group, options, "--reference-fasta", type=Path, help="FASTA with reference IDs matching the table.")
    add_option(group, options, "--reference-table", type=Path, help="CSV, TSV or JSON with reference sequences.")
    add_option(group, options, "--reference-table-id-column", help="Reference table ID column for CSV/TSV.")
    add_option(group, options, "--reference-table-sequence-column", help="Reference table sequence column for CSV/TSV.")
    add_option(group, options, "--reference-file-id-pattern", help="Extract IDs from the FASTA or reference table IDs.")
    add_option(group, options, "--reference-sequence", help="Known reference sequence for an input with one reference.")
    add_option(group, options, "--reference-slice", type=int, nargs=2, metavar=("START", "STOP"), help="Zero-based half-open reference slice in residues or codons.")

    group = clean_parser.add_argument_group("Sequences, annotations and coordinates")
    add_option(group, options, "--variant-sequence-column", help="Column containing observed variant sequences.")
    add_option(group, options, "--mutation-column", help="Column containing substitution annotations.")
    add_option(group, options, "--mutation-parts", nargs=3, metavar=("WT", "POSITION", "MUTANT"), help="Join these three columns into an annotation.")
    add_option(group, options, "--offset", type=int, help="Position of the first reference residue or codon in the annotations.")
    add_option(group, options, "--offset-for", action="append", metavar="REFERENCE=OFFSET", help="Use an offset for a specific reference.")
    add_option(group, options, "--position-map", type=Path, help="Table with reference_id, input_position, reference_index.")
    add_option(group, options, "--mutation-separator", help="Required with substitution annotations. For single substitutions only, pass a comma.")
    add_option(group, options, "--sequence-type", choices=["protein", "dna", "rna"], required=True, help="Interpret the input as amino acids or DNA/RNA codons.")
    add_option(group, options, "--workers", type=bounded_int(1, None), default=4, help="Workers for supported sequence operations (default: %(default)s).")
    add_option(group, options, "--include-ambiguous", action="store_true", help="Allow B, Z, J and X in protein sequences.")
    add_option(group, options, "--include-sec-pyl", action="store_true", help="Allow U and O in protein sequences.")
    add_option(group, options, "--filter-stop", action="store_true", help="Exclude substitutions introducing a DNA/RNA stop codon.")
    add_option(group, options, "--to-protein", action="store_true", help="Translate measured codon variants before grouping and export.")
    add_option(group, options, "--min-substitutions", type=bounded_int(1, None), help="Minimum number of substitutions per record.")
    add_option(group, options, "--max-substitutions", type=bounded_int(1, None), help="Maximum number of substitutions per record.")
    add_option(group, options, "--mutation-effect", choices=["synonymous", "missense", "nonsense"], action="append", help="Retain records whose substitutions all have an allowed effect.")

    group = clean_parser.add_argument_group("Measurements and output grouping")
    add_option(group, options, "--label-column", required=True, action="append", help="Numerical measurement to keep. Repeat for separate outputs.")
    add_option(group, options, "--group-by", action="append", metavar="COLUMN", help="Keep experiment groups separate during aggregation.")
    add_option(group, options, "--split-by", action="append", metavar="COLUMN", help="Write separate outputs for distinct condition combinations.")
    add_option(group, options, "--aggregate", choices=["keep", "mean", "first", "nearest"], help="Required when repeated measurements remain in a group.")
    add_option(group, options, "--assume-comparable", action="store_true", help="Confirm that repeats are comparable without condition columns.")
    add_option(group, options, "--nearest-by", action="append", metavar="COLUMN=TARGET", help="Target condition for nearest, in priority order.")
    add_option(group, options, "--nearest-weight", action="append", metavar="COLUMN=WEIGHT", help="Weight one distance in the native priority-ordered tuple.")
    add_option(group, options, "--wild-type-marker", help="Measured WT annotation. Matching ignores case and edge whitespace.")
    add_option(group, options, "--no-wild-type-marker", action="store_true", help="Declare that no annotation marks measured WT. Choose one WT marker option when using annotations.")
    add_option(group, options, "--subtract-wt", action="store_true", help="Subtract a measured WT baseline before variant aggregation.")
    add_option(group, options, "--wt-from-equal", action="store_true", help="Use equal reference and variant sequences as measured WT.")
    add_option(group, options, "--wt-group-by", action="append", metavar="COLUMN", help="Extra grouping column for matching a WT baseline.")
    add_option(group, options, "--wt-baseline", choices=["unique", "mean"], help="Required if a matched baseline group has multiple WT measurements.")
    add_option(group, options, "--save-format", choices=["mutcleaner", "dataframe", "pickle"], default="mutcleaner", help="Choose a native MutCleaner export format (default: %(default)s).")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="mutcleaner", description="Inspect and clean measured substitution datasets from local files.", allow_abbrev=False)
    parser.add_argument("--version", action="version", version=f"%(prog)s {mutcleaner.__version__}")
    subcommands = parser.add_subparsers(dest="command")
    guide_parser = subcommands.add_parser("guide", help="Show the complete guide or a specific topic, including the agent workflow.", allow_abbrev=False)
    guide_parser.add_argument("topic", nargs="?", choices=("inspect", "clean", "measurements", "output", "agent"))

    inspect_help = "Report evidence from a table, FASTA, or directory without modifying it."
    inspect_parser = subcommands.add_parser("inspect", help=inspect_help, description=inspect_help, allow_abbrev=False)
    inspect_options = {}
    add_inspect_arguments(inspect_parser, inspect_options)

    clean_help = "Clean measured variants and write verified datasets to a new OUTPUT_PATH."
    clean_parser = subcommands.add_parser("clean", help=clean_help, description=clean_help, allow_abbrev=False)
    clean_options = {}
    add_clean_arguments(clean_parser, clean_options)

    for command_parser in (inspect_parser, clean_parser):
        command_parser.register("type", None, unprotect_option_value)
        command_parser.register("type", int, lambda value: int(unprotect_option_value(value)))
        command_parser.register("type", Path, lambda value: Path(unprotect_option_value(value)))

    argv = list(sys.argv[1:] if argv is None else argv)
    option_specs = {"inspect": inspect_options, "clean": clean_options}
    if argv and argv[0] in option_specs:
        argv = [argv[0], *protect_option_values(argv[1:], option_specs[argv[0]])]
    options = vars(parser.parse_args(argv))
    command = options.pop("command")
    if command is None:
        parser.print_help(file=sys.stderr)
        return 2
    if command == "guide":
        sys.stdout.write(guide(options["topic"]))
        return 0

    selected_parser = inspect_parser if command == "inspect" else clean_parser
    custom_paths = options["custom_cleaner"] or ()
    if len(custom_paths) > 1:
        selected_parser.error("argument --custom-cleaner: specify this option at most once.")
    if custom_paths and custom_paths[0].is_dir():
        selected_parser.error(f"argument --custom-cleaner: {custom_paths[0]} is a directory.")

    input_path = options.pop("input_path")
    command_options = option_specs[command]
    for destination, nargs, repeated in command_options.values():
        value = options[destination]
        if repeated:
            options[destination] = tuple(tuple(item) if nargs else item for item in value or ())
        elif nargs and value is not None:
            options[destination] = tuple(value)

    if command == "clean":
        if (options["reference_id_column"] is None) == (options["reference_id"] is None):
            selected_parser.error("Specify exactly one of --reference-id-column or --reference-id.")
        if options["sheet"] is not None and options["sheet_index"] is not None:
            selected_parser.error("Use either --sheet or --sheet-index.")

    try:
        with contextlib.redirect_stdout(sys.stderr):
            if command == "inspect":
                report = inspect_table(input_path, SimpleNamespace(**options))
            else:
                report = clean(input_path, options.pop("output_path"), options, command_options)
    except (ImportError, KeyError, OSError, re.PatternError, RuntimeError, ValueError) as error:
        if command == "clean":
            sys.stdout.write(json.dumps({"status": "error", "error": str(error)}, ensure_ascii=False) + "\n")
        sys.stderr.write(f"Error: {error}\n")
        return 1

    sys.stdout.write(json.dumps(report, ensure_ascii=False, allow_nan=False, indent=2) + "\n")
    if command == "inspect":
        return int(report.get("verification") == "failed")
    return int(report["status"] != "success")
