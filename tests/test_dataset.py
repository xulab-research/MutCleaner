import json
import tempfile
from pathlib import Path

import pandas as pd
import pytest

from mutcleaner.cleaners.basic_cleaners import apply_mutations_to_sequences, convert_to_mutation_dataset_format, validate_mutations
from mutcleaner.core.dataset import MutationDataset
from mutcleaner.core.mutation import AminoAcidMutation, AminoAcidMutationSet, CodonMutation, CodonMutationSet
from mutcleaner.core.sequence import DNASequence, ProteinSequence, RNASequence


class TestMutationDatasetInit:
    """Test MutationDataset initialization."""

    def test_dataset_creation(self):
        """Test dataset creation."""
        dataset = MutationDataset()

        assert dataset.name is None
        assert len(dataset.reference_sequences) == 0
        assert len(dataset.mutation_sets) == 0
        assert len(dataset) == 0

    def test_dataset_creation_with_name(self):
        """Test dataset creation with a name."""
        dataset = MutationDataset(name="test_dataset")

        assert dataset.name == "test_dataset"
        assert len(dataset) == 0


class TestReferenceSequenceManagement:
    """Test reference sequence management."""

    def test_add_reference_sequence(self):
        """Test adding a reference sequence."""
        dataset = MutationDataset()
        protein_seq = ProteinSequence("ACDEFGHIKLMNPQRSTVWY", name="test_protein")

        dataset.add_reference_sequence("ref1", protein_seq)

        assert "ref1" in dataset.reference_sequences
        assert dataset.reference_sequences["ref1"] == protein_seq
        assert dataset.list_reference_sequences() == ["ref1"]

    def test_add_duplicate_reference_sequence(self):
        """Test adding a duplicate reference sequence ID."""
        dataset = MutationDataset()
        protein_seq = ProteinSequence("ACDEFG", name="test_protein")

        dataset.add_reference_sequence("ref1", protein_seq)

        with pytest.raises(
            ValueError, match="Reference sequence with ID 'ref1' already exists"
        ):
            dataset.add_reference_sequence("ref1", protein_seq)

    def test_get_reference_sequence(self):
        """Test retrieving a reference sequence."""
        dataset = MutationDataset()
        protein_seq = ProteinSequence("ACDEFG", name="test_protein")
        dataset.add_reference_sequence("ref1", protein_seq)

        retrieved_seq = dataset.get_reference_sequence("ref1")
        assert retrieved_seq == protein_seq

    def test_get_nonexistent_reference_sequence(self):
        """Test retrieving a nonexistent reference sequence."""
        dataset = MutationDataset()

        with pytest.raises(
            ValueError, match="Reference sequence with ID 'nonexistent' not found"
        ):
            dataset.get_reference_sequence("nonexistent")

    def test_remove_reference_sequence(self):
        """Test removing a reference sequence."""
        dataset = MutationDataset()
        protein_seq = ProteinSequence("ACDEFG", name="test_protein")
        dataset.add_reference_sequence("ref1", protein_seq)

        dataset.remove_reference_sequence("ref1")

        assert "ref1" not in dataset.reference_sequences
        assert len(dataset.list_reference_sequences()) == 0

    def test_remove_referenced_sequence(self):
        """Test removing a reference sequence referenced by a mutation set."""
        dataset = MutationDataset()
        protein_seq = ProteinSequence("ACDEFG", name="test_protein")
        dataset.add_reference_sequence("ref1", protein_seq)

        # Add a mutation set that references this sequence
        mutation = AminoAcidMutation("A", 0, "V")
        mutation_set = AminoAcidMutationSet([mutation])
        dataset.add_mutation_set(mutation_set, "ref1")

        with pytest.raises(
            ValueError, match="Cannot remove sequence 'ref1' as it is referenced"
        ):
            dataset.remove_reference_sequence("ref1")


class TestMutationSetManagement:
    """Test mutation set management."""

    def test_add_mutation_set(self):
        """Test adding a mutation set."""
        dataset = MutationDataset()
        protein_seq = ProteinSequence("ACDEFG", name="test_protein")
        dataset.add_reference_sequence("ref1", protein_seq)

        mutation = AminoAcidMutation("A", 0, "V")
        mutation_set = AminoAcidMutationSet([mutation], name="test_set")

        dataset.add_mutation_set(mutation_set, "ref1", label=1.0)

        assert len(dataset) == 1
        assert dataset.mutation_sets[0] == mutation_set
        assert dataset.get_mutation_set_reference(0) == "ref1"
        assert dataset.get_mutation_set_label(0) == 1.0

    def test_add_mutation_set_without_reference(self):
        """Test adding a mutation set without a reference sequence."""
        dataset = MutationDataset()

        mutation = AminoAcidMutation("A", 0, "V")
        mutation_set = AminoAcidMutationSet([mutation])

        with pytest.raises(
            ValueError, match="Reference sequence with ID 'ref1' not found"
        ):
            dataset.add_mutation_set(mutation_set, "ref1")

    def test_add_multiple_mutation_sets(self):
        """Test adding multiple mutation sets."""
        dataset = MutationDataset()
        protein_seq = ProteinSequence("ACDEFG", name="test_protein")
        dataset.add_reference_sequence("ref1", protein_seq)

        mutation_sets = [
            AminoAcidMutationSet([AminoAcidMutation("A", 0, "V")]),
            AminoAcidMutationSet([AminoAcidMutation("C", 1, "G")]),
            AminoAcidMutationSet([AminoAcidMutation("D", 2, "E")]),
        ]
        reference_ids = ["ref1", "ref1", "ref1"]
        labels = [1.0, 2.0, 3.0]

        dataset.add_mutation_sets(mutation_sets, reference_ids, labels)

        assert len(dataset) == 3
        for i in range(3):
            assert dataset.get_mutation_set_label(i) == labels[i]

    def test_add_mutation_sets_mismatched_lengths(self):
        """Test mismatched lengths when adding mutation sets."""
        dataset = MutationDataset()
        protein_seq = ProteinSequence("ACDEFG", name="test_protein")
        dataset.add_reference_sequence("ref1", protein_seq)

        mutation_sets = [AminoAcidMutationSet([AminoAcidMutation("A", 0, "V")])]
        reference_ids = ["ref1", "ref1"]  # Mismatched length

        with pytest.raises(ValueError, match="Number of reference_ids must match"):
            dataset.add_mutation_sets(mutation_sets, reference_ids)

    def test_set_mutation_set_reference(self):
        """Test setting the reference sequence for a mutation set."""
        dataset = MutationDataset()
        protein_seq1 = ProteinSequence("ACDEFG", name="test_protein1")
        protein_seq2 = ProteinSequence("GHIKLM", name="test_protein2")
        dataset.add_reference_sequence("ref1", protein_seq1)
        dataset.add_reference_sequence("ref2", protein_seq2)

        mutation_set = AminoAcidMutationSet([AminoAcidMutation("A", 0, "V")])
        dataset.add_mutation_set(mutation_set, "ref1")

        # Change the reference sequence
        dataset.set_mutation_set_reference(0, "ref2")
        assert dataset.get_mutation_set_reference(0) == "ref2"

    def test_remove_mutation_set(self):
        """Test removing a mutation set."""
        dataset = MutationDataset()
        protein_seq = ProteinSequence("ACDEFG", name="test_protein")
        dataset.add_reference_sequence("ref1", protein_seq)

        # Add three mutation sets
        for i in range(3):
            mutation_set = AminoAcidMutationSet(
                [AminoAcidMutation("A", 0, f"{chr(82+i)}")]
            )
            dataset.add_mutation_set(mutation_set, "ref1", label=i)

        # Remove the middle mutation set
        dataset.remove_mutation_set(1)

        assert len(dataset) == 2
        assert dataset.get_mutation_set_label(0) == 0
        assert dataset.get_mutation_set_label(1) == 2  # The index has been updated


class TestValidation:
    """Test validation functionality."""

    def test_validate_against_references_valid(self):
        """Test validation of valid mutations."""
        dataset = MutationDataset()
        protein_seq = ProteinSequence("ACDEFG", name="test_protein")
        dataset.add_reference_sequence("ref1", protein_seq)

        # Add a valid mutation
        mutation = AminoAcidMutation("A", 0, "V")
        mutation_set = AminoAcidMutationSet([mutation], name="valid_set")
        dataset.add_mutation_set(mutation_set, "ref1")

        results = dataset.validate_against_references()

        assert len(results["valid_mutation_sets"]) == 1
        assert len(results["invalid_mutation_sets"]) == 0
        assert len(results["position_mismatches"]) == 0

    def test_validate_against_references_position_mismatch(self):
        """Test validation of position mismatches."""
        dataset = MutationDataset()
        protein_seq = ProteinSequence("ACDEFG", name="test_protein")
        dataset.add_reference_sequence("ref1", protein_seq)

        # Add a mutation with a wild-type mismatch
        mutation = AminoAcidMutation("V", 0, "L")  # Position 0 should be A instead of V
        mutation_set = AminoAcidMutationSet([mutation], name="mismatch_set")
        dataset.add_mutation_set(mutation_set, "ref1")

        results = dataset.validate_against_references()

        assert len(results["position_mismatches"]) == 1
        assert results["position_mismatches"][0]["expected"] == "A"
        assert results["position_mismatches"][0]["found"] == "V"

    def test_validate_against_references_out_of_bounds(self):
        """Test validation of out-of-bounds positions."""
        dataset = MutationDataset()
        protein_seq = ProteinSequence("ACDEFG", name="test_protein")
        dataset.add_reference_sequence("ref1", protein_seq)

        # Add an out-of-bounds mutation
        mutation = AminoAcidMutation("A", 10, "V")  # The sequence length is only 6
        mutation_set = AminoAcidMutationSet([mutation], name="out_of_bounds_set")
        dataset.add_mutation_set(mutation_set, "ref1")

        results = dataset.validate_against_references()

        assert len(results["invalid_mutation_sets"]) == 1
        assert "exceeds sequence length" in results["invalid_mutation_sets"][0]["error"]


class TestDataConversion:
    """Test data conversion functionality."""

    def test_to_dataframe(self):
        """Test conversion to a DataFrame."""
        dataset = MutationDataset()
        protein_seq = ProteinSequence("ACDEFG", name="test_protein")
        dataset.add_reference_sequence("ref1", protein_seq)

        # Add a mutation set
        mutations = [AminoAcidMutation("A", 0, "V"), AminoAcidMutation("C", 1, "G")]
        mutation_set = AminoAcidMutationSet(mutations, name="test_set")
        dataset.add_mutation_set(mutation_set, "ref1", label=1.5)

        df = dataset.to_dataframe()

        assert isinstance(df, pd.DataFrame)
        assert len(df) == 2  # Two mutations
        assert df.iloc[0]["mutation_set_name"] == "test_set"
        assert df.iloc[0]["reference_id"] == "ref1"
        assert df.iloc[0]["label"] == 1.5
        assert df.iloc[0]["mutation_string"] == "A0V"
        assert df.iloc[1]["mutation_string"] == "C1G"

    def test_convert_codon_to_amino_acid_sets(self):
        """Test converting codon mutation sets to amino acid mutation sets."""
        dataset = MutationDataset()
        dna_seq = DNASequence("ATGGCCTAT", name="test_dna")
        dataset.add_reference_sequence("ref1", dna_seq)

        # Add a codon mutation set
        codon_mutation = CodonMutation("ATG", 0, "CTG")  # M -> L
        codon_set = CodonMutationSet([codon_mutation])
        dataset.add_mutation_set(codon_set, "ref1", label=2.0)

        # Convert
        converted_dataset = dataset.convert_codon_to_amino_acid_sets(
            convert_labels=True
        )

        assert len(converted_dataset) == 1
        converted_set = converted_dataset.mutation_sets[0]
        assert isinstance(converted_set, AminoAcidMutationSet)
        assert len(converted_set) == 1
        assert str(converted_set.mutations[0]) == "M0L"
        assert converted_dataset.get_mutation_set_label(0) == 2.0


class TestFiltering:
    """Test filtering functionality."""

    def test_filter_by_reference(self):
        """Test filtering by reference sequence."""
        dataset = MutationDataset()

        # Add two reference sequences
        protein_seq1 = ProteinSequence("ACDEFG", name="protein1")
        protein_seq2 = ProteinSequence("GHIKLM", name="protein2")
        dataset.add_reference_sequence("ref1", protein_seq1)
        dataset.add_reference_sequence("ref2", protein_seq2)

        # Add mutation sets for each reference sequence
        for i in range(3):
            mutation_set = AminoAcidMutationSet([AminoAcidMutation("A", 0, "V")])
            dataset.add_mutation_set(mutation_set, "ref1")

        for i in range(2):
            mutation_set = AminoAcidMutationSet([AminoAcidMutation("G", 0, "A")])
            dataset.add_mutation_set(mutation_set, "ref2")

        # Filter
        filtered = dataset.filter_by_reference("ref1")

        assert len(filtered) == 3
        assert "ref1" in filtered.reference_sequences
        assert "ref2" not in filtered.reference_sequences

    def test_filter_by_mutation_type(self):
        """Test filtering by mutation type."""
        dataset = MutationDataset()
        protein_seq = ProteinSequence("ACDEFG", name="test_protein")
        dataset.add_reference_sequence("ref1", protein_seq)

        # Add an amino acid mutation
        aa_mutation = AminoAcidMutation("A", 0, "V")
        aa_set = AminoAcidMutationSet([aa_mutation])
        dataset.add_mutation_set(aa_set, "ref1")

        # Add a DNA sequence and a codon mutation
        dna_seq = DNASequence("ATGGCCTAT", name="test_dna")
        dataset.add_reference_sequence("ref2", dna_seq)
        codon_mutation = CodonMutation("ATG", 0, "CTG")
        codon_set = CodonMutationSet([codon_mutation])
        dataset.add_mutation_set(codon_set, "ref2")

        # Filter amino acid mutations
        filtered = dataset.filter_by_mutation_type(AminoAcidMutation)

        assert len(filtered) == 1
        assert isinstance(filtered.mutation_sets[0], AminoAcidMutationSet)

    def test_filter_by_effect_type(self):
        """Test filtering by effect type."""
        dataset = MutationDataset()
        protein_seq = ProteinSequence("ACDEFG", name="test_protein")
        dataset.add_reference_sequence("ref1", protein_seq)

        # Add mutations with different effect types
        mutations = [
            AminoAcidMutation("A", 0, "V"),  # missense
            AminoAcidMutation("C", 1, "C"),  # synonymous
            AminoAcidMutation("D", 2, "*"),  # nonsense
        ]

        for mut in mutations:
            mutation_set = AminoAcidMutationSet([mut])
            dataset.add_mutation_set(mutation_set, "ref1")

        # Filter missense mutations
        filtered = dataset.filter_by_effect_type("missense")

        assert len(filtered) == 1
        assert filtered.mutation_sets[0].mutations[0].is_missense()


class TestStatistics:
    """Test statistics functionality."""

    def test_get_statistics(self):
        """Test retrieving statistics."""
        dataset = MutationDataset()
        protein_seq = ProteinSequence("ACDEFG", name="test_protein")
        dataset.add_reference_sequence("ref1", protein_seq)

        # Add single-mutation and multi-mutation sets
        single_mutation = AminoAcidMutationSet([AminoAcidMutation("A", 0, "V")])
        dataset.add_mutation_set(single_mutation, "ref1")

        multiple_mutations = AminoAcidMutationSet(
            [AminoAcidMutation("A", 0, "V"), AminoAcidMutation("C", 1, "G")]
        )
        dataset.add_mutation_set(multiple_mutations, "ref1")

        stats = dataset.get_statistics()

        assert stats["total_mutation_sets"] == 2
        assert stats["total_mutations"] == 3
        assert stats["single_mutation_sets"] == 1
        assert stats["multiple_mutation_sets"] == 1
        assert stats["average_mutations_per_set"] == 1.5

    def test_get_position_coverage(self):
        """Test retrieving position coverage statistics."""
        dataset = MutationDataset()
        protein_seq = ProteinSequence("ACDEFG", name="test_protein")
        dataset.add_reference_sequence("ref1", protein_seq)

        # Add mutations covering different positions
        positions = [0, 1, 3, 5]  # Skip positions 2 and 4
        for pos in positions:
            mutation = AminoAcidMutation(protein_seq.get_residue(pos), pos, "A")
            mutation_set = AminoAcidMutationSet([mutation])
            dataset.add_mutation_set(mutation_set, "ref1")

        coverage = dataset.get_position_coverage("ref1")

        assert coverage["sequence_length"] == 6
        assert coverage["covered_positions"] == 4
        assert coverage["uncovered_positions"] == 2
        assert coverage["coverage_percentage"] == pytest.approx(66.67, rel=0.01)
        assert coverage["position_list"] == [0, 1, 3, 5]


class TestSaveLoad:
    """Test save and load functionality."""

    @pytest.mark.parametrize("format_mutations, expected_annotation, next_offset", [(True, "A0V,G2S", 0), (False, "A2V,G4S", 2)])
    def test_offset_through_cleaning_and_save_load(self, format_mutations, expected_annotation, next_offset, tmp_path):
        """Preserve both sites when cleaning and reloading mutations numbered from 2."""
        raw = pd.DataFrame({"name": ["example"], "sequence": ["ACG"], "mut_info": ["A2V,G4S"], "score": [1.0]})
        validated = validate_mutations(raw, format_mutations=format_mutations, offset=2, num_workers=1)
        assert validated.side["failed"].empty
        assert validated.main["mut_info"].tolist() == [expected_annotation]

        applied = apply_mutations_to_sequences(validated.main, offset=next_offset, num_workers=1)
        assert applied.side["failed"].empty
        assert applied.main["mut_seq"].tolist() == ["VCS"]

        variants, references = convert_to_mutation_dataset_format(applied.main, sequence_column="sequence", offset=next_offset)
        assert variants["position"].tolist() == [0, 2]
        dataset = MutationDataset.from_dataframe(variants, references)
        dataset.save_by_reference(tmp_path)

        loaded = MutationDataset.load(str(tmp_path))
        assert len(loaded) == 1
        assert str(loaded.mutation_sets[0]) == "A0V,G2S"

    def test_load_by_reference_with_external_offset(self, tmp_path):
        """Read annotations numbered from 2 from a reference directory."""
        reference_dir = tmp_path / "example"
        reference_dir.mkdir()
        (reference_dir / "wt.fasta").write_text(">example\nACG\n")
        variants = pd.DataFrame({"mutation_name": ["A2V,G4S"], "mutated_sequence": ["VCS"], "label": [1.0]})
        variants.to_csv(reference_dir / "data.csv", index=False)

        dataset = MutationDataset.load_by_reference(tmp_path, offset=2)
        assert len(dataset) == 1
        assert dataset.mutation_sets[0].get_positions() == [0, 2]

    def test_save_load_pickle(self):
        """Test saving and loading in pickle format."""
        # Create a dataset
        dataset = MutationDataset(name="test_dataset")
        protein_seq = ProteinSequence("ACDEFG", name="test_protein")
        dataset.add_reference_sequence("ref1", protein_seq)

        mutation = AminoAcidMutation("A", 0, "V")
        mutation_set = AminoAcidMutationSet([mutation], name="test_set")
        dataset.add_mutation_set(mutation_set, "ref1", label=1.0)

        # Save and load
        with tempfile.NamedTemporaryFile(suffix=".pkl", delete=False) as tmp:
            dataset.save(tmp.name, save_type="pickle")
            loaded_dataset = MutationDataset.load(tmp.name)

        # Validate
        assert loaded_dataset.name == "test_dataset"
        assert len(loaded_dataset) == 1
        assert loaded_dataset.get_mutation_set_label(0) == 1.0

        # Clean up
        Path(tmp.name).unlink()

    def test_save_load_dataframe(self):
        """Test saving and loading in dataframe format."""
        # Create a dataset
        dataset = MutationDataset(name="test_dataset")
        protein_seq = ProteinSequence("ACDEFG", name="test_protein")
        dataset.add_reference_sequence("ref1", protein_seq)

        mutation = AminoAcidMutation("A", 0, "V")
        mutation_set = AminoAcidMutationSet([mutation], name="test_set")
        dataset.add_mutation_set(mutation_set, "ref1", label=2.5)
        dataset.metadata["test_key"] = "test_value"

        # Save and load
        with tempfile.TemporaryDirectory() as tmpdir:
            base_path = Path(tmpdir) / "test_data"
            dataset.save(str(base_path), save_type="dataframe")
            loaded_dataset = MutationDataset.load(str(base_path), load_type="dataframe")

        # Validate
        assert loaded_dataset.name == "test_dataset"
        assert len(loaded_dataset) == 1
        assert loaded_dataset.get_mutation_set_label(0) == 2.5
        assert loaded_dataset.metadata["test_key"] == "test_value"

    def test_save_load_by_reference(self):
        """Test saving and loading by reference sequence."""
        # Create a dataset
        dataset = MutationDataset(name="test_dataset")

        # Add two reference sequences
        protein_seq1 = ProteinSequence("ACDEFG", name="protein1")
        protein_seq2 = ProteinSequence("GHIKLM", name="protein2")
        dataset.add_reference_sequence("ref1", protein_seq1)
        dataset.add_reference_sequence("ref2", protein_seq2)

        # Add multiple mutation sets for the first reference sequence
        mutation1 = AminoAcidMutation("A", 0, "V")
        mutation_set1 = AminoAcidMutationSet([mutation1])
        dataset.add_mutation_set(mutation_set1, "ref1", label=1.0)

        # Add a multi-mutation set
        mutation1_multi = AminoAcidMutation("A", 0, "V")
        mutation2_multi = AminoAcidMutation("C", 1, "G")
        mutation_set1_multi = AminoAcidMutationSet([mutation1_multi, mutation2_multi])
        dataset.add_mutation_set(mutation_set1_multi, "ref1", label=3.0)

        # Add a mutation for the second reference sequence
        mutation2 = AminoAcidMutation("G", 0, "A")
        mutation_set2 = AminoAcidMutationSet([mutation2])
        dataset.add_mutation_set(mutation_set2, "ref2", label=2.0)

        # Add an unlabeled mutation set
        mutation3 = AminoAcidMutation("H", 1, "R")
        mutation_set3 = AminoAcidMutationSet([mutation3])
        dataset.add_mutation_set(mutation_set3, "ref2")  # Unlabeled

        # Save and load
        with tempfile.TemporaryDirectory() as tmpdir:
            dataset.save_by_reference(tmpdir)

            # Validate the saved file structure
            base_path = Path(tmpdir)
            ref_dirs = [d for d in base_path.iterdir() if d.is_dir()]
            assert (
                len(ref_dirs) == 2
            ), f"Expected 2 reference directories, got {len(ref_dirs)}"

            # Validate that each directory contains required files
            for ref_dir in ref_dirs:
                assert (
                    ref_dir / "data.csv"
                ).exists(), f"Missing data.csv in {ref_dir.name}"
                assert (
                    ref_dir / "wt.fasta"
                ).exists(), f"Missing wt.fasta in {ref_dir.name}"
                assert (
                    ref_dir / "metadata.json"
                ).exists(), f"Missing metadata.json in {ref_dir.name}"

                # Validate the metadata.json format
                with open(ref_dir / "metadata.json") as f:
                    metadata = json.load(f)
                    required_fields = [
                        "reference_id",
                        "sequence_name",
                        "sequence_type",
                        "sequence_length",
                        "num_mutation_sets",
                        "total_mutations",
                        "covered_positions",
                        "coverage_percentage",
                        "num_unique_labels",
                        "has_unlabeled",
                        "dataset_name",
                    ]
                    for field in required_fields:
                        assert field in metadata, f"Missing field {field} in metadata"

            # Load the dataset
            loaded_dataset = MutationDataset.load_by_reference(tmpdir, "loaded_dataset", offset=0)

        # Validate basic attributes
        assert loaded_dataset.name == "loaded_dataset"
        assert (
            len(loaded_dataset) == 4
        ), f"Expected 4 mutation sets, got {len(loaded_dataset)}"
        assert len(loaded_dataset.reference_sequences) == 2
        assert "ref1" in loaded_dataset.reference_sequences
        assert "ref2" in loaded_dataset.reference_sequences

        # Validate reference sequence contents
        loaded_seq1 = loaded_dataset.reference_sequences["ref1"]
        loaded_seq2 = loaded_dataset.reference_sequences["ref2"]
        assert str(loaded_seq1) == "ACDEFG"
        assert str(loaded_seq2) == "GHIKLM"
        assert loaded_seq1.name == "protein1"
        assert loaded_seq2.name == "protein2"

        # Validate mutation sets and labels
        ref1_sets = []
        ref2_sets = []

        for i, (mutation_set, ref_id) in enumerate(loaded_dataset):
            label = loaded_dataset.mutation_set_labels.get(i, "")

            if ref_id == "ref1":
                ref1_sets.append((mutation_set, label))
            elif ref_id == "ref2":
                ref2_sets.append((mutation_set, label))

        # Validate mutation sets for ref1
        assert (
            len(ref1_sets) == 2
        ), f"Expected 2 mutation sets for ref1, got {len(ref1_sets)}"

        # Find single and multiple mutations
        single_mut = None
        multi_mut = None
        for mut_set, label in ref1_sets:
            if len(mut_set) == 1:
                single_mut = (mut_set, label)
            elif len(mut_set) == 2:
                multi_mut = (mut_set, label)

        assert single_mut is not None, "Single mutation not found for ref1"
        assert multi_mut is not None, "Multi mutation not found for ref1"
        assert single_mut[1] == 1.0, f"Expected label 1.0, got {single_mut[1]}"
        assert multi_mut[1] == 3.0, f"Expected label '3.0', got {multi_mut[1]}"

        # Validate mutation sets for ref2
        assert (
            len(ref2_sets) == 2
        ), f"Expected 2 mutation sets for ref2, got {len(ref2_sets)}"

        # Check labels
        labels = [label for _, label in ref2_sets]
        assert 2.0 in labels, "Label 2.0 not found in ref2"
        assert pd.isna(labels).any(), "Empty label (unlabeled) not found in ref2"

        # Validate statistics
        stats = loaded_dataset.get_statistics()
        assert stats["total_mutation_sets"] == 4
        assert stats["total_mutations"] == 5  # 1 + 2 + 1 + 1
        assert stats["num_reference_sequences"] == 2

        print("All tests passed!")

    def test_save_load_error_handling(self):
        """Test error handling."""
        # Test an empty directory
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                MutationDataset.load_by_reference(tmpdir, offset=0)
                assert False, "Should raise ValueError for empty directory"
            except ValueError as e:
                assert "No reference directories found" in str(e)

        # Test a nonexistent directory
        try:
            MutationDataset.load_by_reference("/nonexistent/path", offset=0)
            assert False, "Should raise FileNotFoundError"
        except FileNotFoundError:
            pass

        print("Error handling tests passed!")

    def test_sanitize_filename(self):
        """Test filename sanitization."""
        # Test filenames containing special characters
        assert MutationDataset._sanitize_filename("test:name") == "test_name"
        assert MutationDataset._sanitize_filename("test/name") == "test_name"
        assert MutationDataset._sanitize_filename("test<>name") == "test__name"

        # Test an empty string
        assert MutationDataset._sanitize_filename("") == "unnamed"

        # Test an overly long filename
        long_name = "a" * 300
        sanitized = MutationDataset._sanitize_filename(long_name)
        assert len(sanitized) == 200


class TestSpecialMethods:
    """Test special methods."""

    def test_len(self):
        """Test the __len__ method."""
        dataset = MutationDataset()
        assert len(dataset) == 0

        protein_seq = ProteinSequence("ACDEFG", name="test_protein")
        dataset.add_reference_sequence("ref1", protein_seq)

        for i in range(5):
            mutation = AminoAcidMutation("A", 0, "V")
            mutation_set = AminoAcidMutationSet([mutation])
            dataset.add_mutation_set(mutation_set, "ref1")

        assert len(dataset) == 5

    def test_iter(self):
        """Test the __iter__ method."""
        dataset = MutationDataset()
        protein_seq = ProteinSequence("ACDEFG", name="test_protein")
        dataset.add_reference_sequence("ref1", protein_seq)

        # Add a mutation set
        mutation_sets_added = []
        for i in range(3):
            mutation = AminoAcidMutation("A", i, "V")
            mutation_set = AminoAcidMutationSet([mutation], name=f"set_{i}")
            dataset.add_mutation_set(mutation_set, "ref1")
            mutation_sets_added.append(mutation_set)

        # Test iteration
        for i, (mutation_set, ref_id) in enumerate(dataset):
            assert mutation_set.name == f"set_{i}"
            assert ref_id == "ref1"

    def test_str(self):
        """Test the __str__ method."""
        dataset = MutationDataset(name="test_dataset")
        protein_seq = ProteinSequence("ACDEFG", name="test_protein")
        dataset.add_reference_sequence("ref1", protein_seq)

        # Add mutations
        mutations = [AminoAcidMutation("A", 0, "V"), AminoAcidMutation("C", 1, "G")]
        mutation_set = AminoAcidMutationSet(mutations)
        dataset.add_mutation_set(mutation_set, "ref1")

        str_repr = str(dataset)
        assert "MutationDataset(test_dataset)" in str_repr
        assert "1 reference sequences" in str_repr
        assert "1 mutation sets" in str_repr
        assert "2 mutations" in str_repr


class TestEdgeCases:
    """Test edge cases."""

    def test_empty_dataset_operations(self):
        """Test operations on an empty dataset."""
        dataset = MutationDataset()

        # Test statistics for an empty dataset
        stats = dataset.get_statistics()
        assert stats["total_mutation_sets"] == 0
        assert stats["total_mutations"] == 0

        # Test DataFrame conversion for an empty dataset
        df = dataset.to_dataframe()
        assert df.empty

        # Test validation for an empty dataset
        results = dataset.validate_against_references()
        assert len(results["valid_mutation_sets"]) == 0

    def test_mutation_set_index_out_of_range(self):
        """Test out-of-range mutation set indices."""
        dataset = MutationDataset()

        with pytest.raises(ValueError, match="Mutation set index 0 out of range"):
            dataset.get_mutation_set_reference(0)

        with pytest.raises(ValueError, match="Mutation set index 0 out of range"):
            dataset.get_mutation_set_label(0)

        with pytest.raises(ValueError, match="Mutation set index 0 out of range"):
            dataset.set_mutation_set_label(0, 1.0)

    def test_from_dataframe_validation(self):
        """Test validation when creating from a DataFrame."""
        # Empty DataFrame
        with pytest.raises(ValueError, match="DataFrame cannot be empty"):
            MutationDataset.from_dataframe(pd.DataFrame(), {})

        # Missing required columns
        df = pd.DataFrame({"col1": [1, 2, 3]})
        with pytest.raises(ValueError, match="DataFrame missing required columns"):
            MutationDataset.from_dataframe(df, {})

    def test_mixed_sequence_types(self):
        """Test a dataset with mixed sequence types."""
        dataset = MutationDataset()

        # Add different types of sequences
        protein_seq = ProteinSequence("ACDEFG", name="protein")
        dna_seq = DNASequence("ATGGCCTAT", name="dna")
        rna_seq = RNASequence("AUGGCCUAU", name="rna")

        dataset.add_reference_sequence("protein_ref", protein_seq)
        dataset.add_reference_sequence("dna_ref", dna_seq)
        dataset.add_reference_sequence("rna_ref", rna_seq)

        # Add corresponding mutations
        aa_mutation = AminoAcidMutation("A", 0, "V")
        aa_set = AminoAcidMutationSet([aa_mutation])
        dataset.add_mutation_set(aa_set, "protein_ref")

        codon_mutation = CodonMutation("ATG", 0, "CTG")
        codon_set = CodonMutationSet([codon_mutation])
        dataset.add_mutation_set(codon_set, "dna_ref")

        # Validate
        assert len(dataset) == 2
        assert len(dataset.reference_sequences) == 3

        # Test statistics
        stats = dataset.get_statistics()
        assert stats["num_reference_sequences"] == 3
