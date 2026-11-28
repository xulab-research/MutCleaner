import ast
import hashlib
import json
import os
import runpy
import subprocess
import sys
import pandas as pd
import pytest

from pathlib import Path
from types import SimpleNamespace

from mutcleaner.cli.main import main


@pytest.fixture
def invoke_cli(capsys):
    def invoke(arguments):
        try:
            exit_code = main(arguments)
        except SystemExit as error:
            exit_code = error.code
        capture = capsys.readouterr()
        return SimpleNamespace(exit_code=exit_code, stdout=capture.out, stderr=capture.err, output=capture.out + capture.err)

    return invoke


@pytest.fixture
def raw_path(tmp_path):
    path = tmp_path / "raw.csv"
    path.write_text("protein_id,reference,encoded,discard\nR1,ACD,sub(A;1;V),x\nR_SKIP,ACD,sub(D;3;N),y\nR2,ACD,sub(C;2;V),z\n", encoding="utf-8")
    return path


def clean_options(label_column):
    return [
        "--reference-id-column",
        "protein_id",
        "--reference-sequence-column",
        "reference",
        "--sequence-type",
        "protein",
        "--mutation-column",
        "mutation",
        "--mutation-separator",
        ",",
        "--no-wild-type-marker",
        "--offset",
        "1",
        "--label-column",
        label_column,
        "--workers",
        "1",
    ]


def test_custom_inspection_cleaning_and_replay(raw_path, tmp_path, monkeypatch, invoke_cli):
    script_path = tmp_path / "dataset_rules.py"
    script_source = '''# coding: latin-1
def custom_cleaner(df):
    """Retain café data."""
    print("custom ran")
    selected = df.iloc[[2, 0, 2, 0]].drop(columns="discard")
    mutation = selected["encoded"].str.replace(r"^sub\\(([A-Z]);(\\d+);([A-Z])\\)$", r"\\1\\2\\3", regex=True)
    return selected.assign(mutation=mutation, measurement=[1.5, None, 3.5, 6.0])
'''.encode("latin-1")
    script_path.write_bytes(script_source)
    digest = hashlib.sha256(script_source).hexdigest()
    inspect_flags = ["--sequence-type", "protein", "--mutation-column", "mutation", "--reference-id-column", "protein_id", "--label-column", "measurement"]
    result = invoke_cli(["inspect", str(raw_path), "--custom-cleaner", str(script_path), "--profile-rows", "1", *inspect_flags])
    assert result.exit_code == 0, result.output
    assert result.stderr.count("custom ran") == 1
    report = json.loads(result.stdout)
    assert report["row_count"] == 4
    assert report["profiled_rows"] == 1
    assert report["sample_rows"][0]["mutation"] == "C2V"
    assert report["candidate_fields"]["profiled_rows"] == 1
    assert report["custom"]["input_rows"] == 3
    assert report["custom"]["output_rows"] == 4
    assert report["custom"]["added_columns"] == ["mutation", "measurement"]
    assert report["custom"]["removed_columns"] == ["discard"]
    assert report["custom"]["sha256"] == digest

    result = invoke_cli(["inspect", str(raw_path), "--custom-cleaner", str(script_path)])
    report = json.loads(result.stdout)
    assert report["profiled_rows"] == 4
    assert report["full_table_profiled"]
    assert report["sample_rows"][1]["measurement"] == ""

    output_path = tmp_path / "output"
    replacements = ["--replace=protein_id", "^R", "--R", "--replace", "protein_id", "^--", ""]
    result = invoke_cli(["clean", str(raw_path), str(output_path), "--custom-cleaner", str(script_path), "--rename", "measurement=score", *replacements, "--aggregate", "mean", "--assume-comparable", *clean_options("score")])
    assert result.exit_code == 0, result.output
    assert result.stderr.count("custom ran") == 1
    report = json.loads(result.stdout)
    assert report["custom"]["sha256"] == digest
    assert report["custom"]["source_row_basis"].startswith("One-based")
    assert report["steps"][:3] == ["read", "custom", "prepare_fields"]
    assert report["global_statistics"]["input_rows"] == 4
    assert "custom_cleaner" not in report["effective_parameters"]
    assert pd.read_csv(output_path / "R2" / "data.csv")["label"].tolist() == [2.5]
    assert pd.read_csv(output_path / "R1" / "data.csv")["label"].tolist() == [6.0]

    script_path.rename(tmp_path / "source_no_longer_at_original_path.py")
    original_run = subprocess.run
    source_root = Path(__file__).resolve().parents[1]
    environment = dict(os.environ, PYTHONPATH=os.pathsep.join(filter(None, [str(source_root), os.environ.get("PYTHONPATH")])))
    restored_paths = []

    def run_mutcleaner(command, *, check):
        assert command[:2] == ["mutcleaner", "clean"]
        restored_path = Path(command[command.index("--custom-cleaner") + 1])
        assert restored_path.read_bytes() == script_source
        restored_paths.append(restored_path)
        result = original_run([sys.executable, "-m", "mutcleaner.cli", *command[1:]], check=check, capture_output=True, text=True, env=environment)
        assert json.loads(result.stdout)["custom"]["sha256"] == digest
        return result

    monkeypatch.setattr(subprocess, "run", run_mutcleaner)
    for old_output, new_output in [(output_path, tmp_path / "replayed"), (tmp_path / "replayed", tmp_path / "replayed_again")]:
        cleaner = old_output / "R1" / "cleaner.py"
        generated = cleaner.read_text(encoding="utf-8")
        assert digest in generated
        assert str(script_path) not in generated
        source_assignment = next(node for node in ast.parse(generated).body if isinstance(node, ast.Assign) and node.targets[0].id == "custom_source")
        assert ast.literal_eval(source_assignment.value) == script_source
        monkeypatch.setattr(sys, "argv", [str(cleaner), str(raw_path), str(new_output)])
        runpy.run_path(str(cleaner), run_name="__main__")
        assert not restored_paths[-1].exists()
        assert {path.name for path in new_output.iterdir()} == {"R1", "R2"}
        for identifier in ("R1", "R2"):
            for filename in ("data.csv", "wt.fasta", "metadata.json"):
                assert (new_output / identifier / filename).read_bytes() == (output_path / identifier / filename).read_bytes()
            assert {path.name for path in (new_output / identifier).iterdir()} == {"data.csv", "wt.fasta", "metadata.json", "cleaner.py"}


@pytest.mark.parametrize(
    ("source", "message"),
    [
        (None, "No such file"),
        ("pass", "define a callable custom_cleaner"),
        ("def custom_cleaner(df):\n    return [df]", "must return a pandas DataFrame"),
        ("def custom_cleaner(df):\n    raise KeyError('missing field')", "missing field"),
        ("def custom_cleaner(df):\n    return df.set_axis(['id', 'id', 'value', 'other'], axis=1)", "must be unique"),
        ("def custom_cleaner(df):\n    return df.set_axis([1, '1', 'value', 'other'], axis=1)", "unique after conversion"),
        ("def custom_cleaner(df):\n    return df.set_axis([' ', 'id', 'value', 'other'], axis=1)", "must be nonempty"),
        ("def custom_cleaner(df):\n    return df.set_axis([None, 'id', 'value', 'other'], axis=1)", "must be nonempty"),
        ("def custom_cleaner(df):\n    return ???", "invalid syntax"),
        ("raise SystemExit(0)", "Custom cleaner"),
    ],
)
def test_custom_failures_are_reported(raw_path, tmp_path, source, message, invoke_cli):
    script_path = tmp_path / "rules.py"
    if source is not None:
        script_path.write_text(source, encoding="utf-8")
    result = invoke_cli(["inspect", str(raw_path), "--custom-cleaner", str(script_path)])
    assert result.exit_code == 1 and message in result.stderr
    output_path = tmp_path / "output"
    result = invoke_cli(["clean", str(raw_path), str(output_path), "--custom-cleaner", str(script_path), *clean_options("score")])
    assert result.exit_code == 1 and message in json.loads(result.stdout)["error"]
    assert not output_path.exists()


def test_custom_requires_a_single_table(raw_path, tmp_path, invoke_cli):
    script_path = tmp_path / "rules.py"
    script_path.write_text("def custom_cleaner(df):\n    return df\n", encoding="utf-8")
    options = ["--custom-cleaner", str(script_path)]
    fasta_path = tmp_path / "reference.fasta"
    fasta_path.write_text(">R\nACD\n", encoding="utf-8")
    workbook_path = tmp_path / "input.xlsx"
    with pd.ExcelWriter(workbook_path, engine="openpyxl") as writer:
        for sheet in ["a", "b"]:
            pd.DataFrame({"reference": ["ACD"]}).to_excel(writer, sheet_name=sheet, index=False)
    for command, message in [
        (["inspect", str(fasta_path), *options], "requires a table"),
        (["inspect", str(tmp_path), "--verify", *options], "requires a table"),
        (["inspect", str(workbook_path), *options], "requires a selected worksheet"),
        (["inspect", str(raw_path), *options, *options], "at most once"),
        (["clean", str(raw_path), str(tmp_path / "duplicate"), *options, *options, *clean_options("score")], "at most once"),
        (["clean", str(fasta_path), str(tmp_path / "fasta_output"), *options, *clean_options("score")], "requires a supported table"),
        (["clean", str(workbook_path), str(tmp_path / "xlsx_output"), *options, *clean_options("score")], "Select --sheet"),
    ]:
        result = invoke_cli(command)
        assert result.exit_code != 0 and message in result.output, result.output
    result = invoke_cli(["inspect", str(workbook_path), "--sheet", "b", *options])
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout)["custom"]["output_rows"] == 1


def test_clean_without_custom(tmp_path, invoke_cli):
    source_path = tmp_path / "plain.csv"
    source_path.write_text("protein_id,reference,mutation,score\nR1,ACD,C2V,1.5\n", encoding="utf-8")
    output_path = tmp_path / "plain_result"
    inspection = invoke_cli(["inspect", str(source_path)])
    assert "custom" not in json.loads(inspection.stdout)
    result = invoke_cli(["clean", str(source_path), str(output_path), *clean_options("score")])
    assert result.exit_code == 0, result.output
    report = json.loads(result.stdout)
    assert "custom" not in report and "custom" not in report["steps"]
    assert "custom_source" not in (output_path / "R1" / "cleaner.py").read_text(encoding="utf-8")
    result = invoke_cli(["inspect", str(output_path), "--verify"])
    assert json.loads(result.stdout)["verification"] == "passed"


def test_inspect_counts_each_repeated_annotation_and_compares_its_own_variant(tmp_path, invoke_cli):
    source_path = tmp_path / "variants.csv"
    source_path.write_text("protein_id,reference,mutation,variant,score\nR1,ACD,A1V,VCD,1.5\nR1,ACD,A1V,CCD,2.5\nR2,ACE,A1V,VCE,3.5\n", encoding="utf-8")
    result = invoke_cli(["inspect", str(source_path), "--sequence-type", "protein", "--reference-id-column", "protein_id", "--reference-sequence-column", "reference", "--mutation-column", "mutation", "--variant-sequence-column", "variant", "--label-column", "score", "--offset", "1", "--workers", "1"])
    assert result.exit_code == 0, result.output
    (candidate,) = json.loads(result.stdout)["candidate_fields"]["offset_candidates"]
    assert candidate["parsed_rows"] == candidate["reference_match_rows"] == 3
    assert candidate["variant_rows_checked"] == 3
    assert candidate["variant_match_rows"] == 2
    assert candidate["repeats"]["groups_with_repeats"] == 1
    assert candidate["repeats"]["groups_with_different_measurements"] == {"score": 1}


def test_clean_rejects_stop_and_mismatched_protein_sequences(tmp_path, invoke_cli):
    source_path = tmp_path / "variants.csv"
    source_path.write_text("protein_id,reference,mutation,variant,score\nR1,ACD,C2*,A*D,1.5\nR1,ACD,D3V,ACA,2.5\nR1,ACD,C2V,AVD,3.5\n", encoding="utf-8")
    output_path = tmp_path / "result"
    result = invoke_cli(["clean", str(source_path), str(output_path), *clean_options("score"), "--variant-sequence-column", "variant"])
    assert result.exit_code == 0, result.output
    report = json.loads(result.stdout)
    assert report["global_statistics"]["invalid_rows"] == 2
    assert any("Provided mutated sequence does not match" in reason for reason in report["global_statistics"]["invalid_reasons"])
    exported = pd.read_csv(output_path / "R1" / "data.csv")
    assert exported[["mutation_name", "mutated_sequence", "label"]].values.tolist() == [["C1V", "AVD", 3.5]]
    result = invoke_cli(["inspect", str(output_path), "--verify"])
    assert json.loads(result.stdout)["verification"] == "passed"


@pytest.mark.parametrize(
    ("sequence_type", "reference", "stop_variant", "missense_variant", "stop_annotation", "missense_annotation"),
    [
        ("dna", "ATGGCT", "ATGTAA", "ATGGTT", "GCT1TAA", "GCT1GTT"),
        ("rna", "AUGGCU", "AUGUAA", "AUGGUU", "GCU1UAA", "GCU1GUU"),
    ],
)
def test_clean_codon_stop_output_and_protein_validation(tmp_path, invoke_cli, sequence_type, reference, stop_variant, missense_variant, stop_annotation, missense_annotation):
    source_path = tmp_path / "codons.csv"
    source_path.write_text(f"reference,variant,score\n{reference},{stop_variant},1.5\n{reference},{missense_variant},2.5\n", encoding="utf-8")
    options = ["--sequence-type", sequence_type, "--reference-id", "R1", "--reference-sequence-column", "reference", "--variant-sequence-column", "variant", "--label-column", "score", "--workers", "1"]

    for directory, extra_options, expected, invalid_rows, filtered_rows in (
        ("codon", [], [[stop_annotation, stop_variant, 1.5], [missense_annotation, missense_variant, 2.5]], 0, 0),
        ("codon_filtered", ["--filter-stop"], [[missense_annotation, missense_variant, 2.5]], 0, 1),
        ("protein", ["--to-protein"], [["A1V", "MV", 2.5]], 1, 0),
    ):
        output_path = tmp_path / directory
        result = invoke_cli(["clean", str(source_path), str(output_path), *options, *extra_options])
        assert result.exit_code == 0, result.output
        statistics = json.loads(result.stdout)["global_statistics"]
        assert statistics["invalid_rows"] == invalid_rows
        assert statistics["filtered_rows"] == filtered_rows
        exported = pd.read_csv(output_path / "R1" / "data.csv")
        assert exported[["mutation_name", "mutated_sequence", "label"]].values.tolist() == expected
        result = invoke_cli(["inspect", str(output_path), "--verify"])
        assert json.loads(result.stdout)["verification"] == "passed"
