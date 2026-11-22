import json
import pandas as pd
import pytest

from mutcleaner.core.dataset import MutationDataset
from mutcleaner.core.mutation import AminoAcidMutation, AminoAcidMutationSet
from mutcleaner.core.sequence import ProteinSequence


class TestSaveByReferenceCoverageMetadata:
    """Test that metadata coverage separates submitted and successful records."""

    def test_metadata_reports_successful_rows_not_submitted_rows(self, tmp_path):
        """Coverage must describe data.csv, while submitted_ fields keep the input scope."""
        dataset = MutationDataset(name="coverage")
        dataset.add_reference_sequence("ref1", ProteinSequence("ACDEFG", name="coverage_protein"))
        dataset.add_mutation_set(
            AminoAcidMutationSet([AminoAcidMutation("A", 0, "V"), AminoAcidMutation("C", 1, "G"), AminoAcidMutation("D", 2, "A")]),
            "ref1",
            label=1.0,
        )
        dataset.add_mutation_set(
            AminoAcidMutationSet([AminoAcidMutation("A", 3, "V"), AminoAcidMutation("A", 4, "V"), AminoAcidMutation("A", 5, "V")]),
            "ref1",
            label=2.0,
        )
        dataset.save(str(tmp_path / "output"), save_type="mutcleaner")
        reference_dir = tmp_path / "output" / "ref1"
        metadata = json.loads((reference_dir / "metadata.json").read_text())
        data_rows = pd.read_csv(reference_dir / "data.csv")
        failed_rows = pd.read_csv(reference_dir / "failed_data.csv")

        assert metadata["covered_positions"] == 3
        assert metadata["uncovered_positions"] == 3
        assert metadata["coverage_percentage"] == pytest.approx(50.0)
        assert metadata["submitted_covered_positions"] == 6
        assert metadata["submitted_coverage_percentage"] == pytest.approx(100.0)
        assert metadata["num_mutation_sets"] == 2
        assert metadata["num_mutation_sets"] == metadata["num_successful_mutation_sets"] + metadata["num_failed_mutation_sets"]
        assert len(data_rows) == metadata["num_successful_mutation_sets"]
        assert len(failed_rows) == metadata["num_failed_mutation_sets"]

    def test_metadata_coverage_equals_submitted_coverage_without_failures(self, tmp_path):
        """Successful and submitted coverage agree when every row is exported."""
        dataset = MutationDataset(name="coverage")
        dataset.add_reference_sequence("ref1", ProteinSequence("ACDEFG", name="coverage_protein"))
        dataset.add_mutation_set(
            AminoAcidMutationSet([AminoAcidMutation("A", 0, "V"), AminoAcidMutation("C", 1, "G"), AminoAcidMutation("D", 2, "A")]),
            "ref1",
            label=1.0,
        )
        dataset.add_mutation_set(
            AminoAcidMutationSet([AminoAcidMutation("E", 3, "A"), AminoAcidMutation("F", 4, "A"), AminoAcidMutation("G", 5, "A")]),
            "ref1",
            label=2.0,
        )
        dataset.save(str(tmp_path / "output"), save_type="mutcleaner")
        reference_dir = tmp_path / "output" / "ref1"
        metadata = json.loads((reference_dir / "metadata.json").read_text())
        data_rows = pd.read_csv(reference_dir / "data.csv")

        assert metadata["num_failed_mutation_sets"] == 0
        assert metadata["num_successful_mutation_sets"] == len(data_rows)
        assert metadata["covered_positions"] == metadata["submitted_covered_positions"] == 6
        assert metadata["coverage_percentage"] == metadata["submitted_coverage_percentage"] == pytest.approx(100.0)
        assert not (reference_dir / "failed_data.csv").exists()
