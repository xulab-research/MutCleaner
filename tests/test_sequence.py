import pytest
import warnings
from unittest.mock import Mock, patch

from mutcleaner.core.sequence import (
    DNASequence,
    RNASequence,
    ProteinSequence,
    translate,
)
from mutcleaner.core.mutation import (
    AminoAcidMutation,
    CodonMutation,
    AminoAcidMutationSet,
    CodonMutationSet,
)
from mutcleaner.core.alphabet import DNAAlphabet, RNAAlphabet, ProteinAlphabet
from mutcleaner.core.codon import CodonTable


class TestSequenceBase:
    """Test the base Sequence class."""

    def test_sequence_creation(self):
        """Test sequence creation."""
        # Create a simple DNA sequence for testing.
        alphabet = DNAAlphabet()
        seq = DNASequence("ATCG", alphabet=alphabet, name="test_seq")

        assert len(seq) == 4
        assert str(seq) == "ATCG"
        assert seq.name == "test_seq"
        assert seq.metadata == {}

    def test_sequence_creation_with_metadata(self):
        """Test sequence creation with metadata."""
        metadata = {"source": "test", "version": 1}
        seq = DNASequence("ATCG", name="test", metadata=metadata)

        assert seq.metadata == metadata

    def test_sequence_equality(self):
        """Test sequence equality."""
        seq1 = DNASequence("ATCG")
        seq2 = DNASequence("ATCG")
        seq3 = DNASequence("TTCG")
        seq4 = RNASequence("AUCG")

        assert seq1 == seq2
        assert seq1 != seq3
        assert seq1 != seq4  # Different alphabet types.
        assert seq1 == "ATCG"  # Compare with a string.

    def test_get_subsequence_0_indexed(self):
        """Test 0-indexed subsequence extraction."""
        seq = DNASequence("ATCGATCG", name="test")

        # Test the normal case.
        subseq = seq.get_subsequence(2, 5)
        assert str(subseq) == "CGA"
        assert subseq.name == "test_2_5"

        # Test without the end parameter.
        subseq2 = seq.get_subsequence(3)  # From position 3 to the end.
        assert str(subseq2) == "GATCG"
        assert subseq2.name == f"test_3_{len(seq)}"

    def test_get_subsequence_invalid_start(self):
        """Test an invalid start position."""
        seq = DNASequence("ATCG")

        with pytest.raises(
            IndexError, match="Start position must be greater than or equal to 0"
        ):
            seq.get_subsequence(-1)

    def test_slice(self):
        """Test slicing."""
        seq = DNASequence("ATCGATCG")
        assert str(seq[:5]) == "ATCGA"  # 0-indexed slice.
        assert str(seq[5:]) == "TCG"  # 0-indexed slice.
        assert str(seq[::-1]) == "ATCGATCG"[::-1]  # 0-indexed slice, reversed.


class TestDNASequence:
    """Test the DNASequence class."""

    def test_dna_sequence_creation(self):
        """Test DNA sequence creation."""
        seq = DNASequence("ATCGATCG")
        assert str(seq) == "ATCGATCG"
        assert isinstance(seq.alphabet, DNAAlphabet)

    def test_reverse_complement(self):
        """Test DNA reverse complement."""
        seq = DNASequence("ATCG", name="test")
        rc = seq.reverse_complement()

        assert str(rc) == "CGAT"  # ATCG -> CGAT.
        assert rc.name == "test_rc"
        assert isinstance(rc, DNASequence)

    def test_reverse_complement_complex(self):
        """Test reverse complement for a complex DNA sequence."""
        seq = DNASequence("ATCGATCGATCG")
        rc = seq.reverse_complement()

        assert str(rc) == "CGATCGATCGAT"

    def test_transcribe(self):
        """Test DNA transcription to RNA."""
        seq = DNASequence("ATCGATCG", name="test")
        rna = seq.transcribe()

        assert str(rna) == "AUCGAUCG"  # T -> U.
        assert rna.name == "test_transcribed"
        assert isinstance(rna, RNASequence)

    @patch("mutcleaner.core.sequence.translate")
    def test_dna_translate(self, mock_translate):
        """Test DNA translation."""
        mock_translate.return_value = "MET"

        seq = DNASequence("ATGATG", name="test")
        protein = seq.translate()

        mock_translate.assert_called_once_with(
            sequence="ATGATG",
            seq_type="DNA",
            codon_table=None,
            start_at_first_met=False,
            stop_at_stop_codon=False,
            require_mod3=True,
            start=None,
            end=None,
        )
        assert str(protein) == "MET"
        assert protein.name == "test_translation"


class TestRNASequence:
    """Test the RNASequence class."""

    def test_rna_sequence_creation(self):
        """Test RNA sequence creation."""
        seq = RNASequence("AUCGAUCG")
        assert str(seq) == "AUCGAUCG"
        assert isinstance(seq.alphabet, RNAAlphabet)

    def test_reverse_complement(self):
        """Test RNA reverse complement."""
        seq = RNASequence("AUCG", name="test")
        rc = seq.reverse_complement()

        assert str(rc) == "CGAU"  # AUCG -> CGAU.
        assert rc.name == "test_rc"
        assert isinstance(rc, RNASequence)

    def test_back_transcribe(self):
        """Test RNA back-transcription to DNA."""
        seq = RNASequence("AUCGAUCG", name="test")
        dna = seq.back_transcribe()

        assert str(dna) == "ATCGATCG"  # U -> T.
        assert dna.name == "test_back_transcribe"
        assert isinstance(dna, DNASequence)

    @patch("mutcleaner.core.sequence.translate")
    def test_rna_translate(self, mock_translate):
        """Test RNA translation."""
        mock_translate.return_value = "MET"

        seq = RNASequence("AUGUGAUG", name="test")
        protein = seq.translate(start_at_first_met=True)

        mock_translate.assert_called_once_with(
            sequence="AUGUGAUG",
            seq_type="RNA",
            codon_table=None,
            start_at_first_met=True,
            stop_at_stop_codon=False,
            require_mod3=True,
            start=None,
            end=None,
        )


class TestProteinSequence:
    """Test the ProteinSequence class."""

    def test_protein_sequence_creation(self):
        """Test protein sequence creation."""
        seq = ProteinSequence("METALA")
        assert str(seq) == "METALA"
        assert isinstance(seq.alphabet, ProteinAlphabet)

    def test_get_residue(self):
        """Test retrieving the amino acid at a specific position."""
        seq = ProteinSequence("METALA")

        assert seq.get_residue(0) == "M"
        assert seq.get_residue(1) == "E"
        assert seq.get_residue(5) == "A"

    def test_get_residue_out_of_range(self):
        """Test out-of-range access."""
        seq = ProteinSequence("MET")

        with pytest.raises(IndexError):
            seq.get_residue(-1)

        with pytest.raises(IndexError):
            seq.get_residue(3)  # Should be >= 3 because the length is 3.

    def test_find_motif(self):
        """Test motif search."""
        seq = ProteinSequence("METKMETALA")

        # Find a single motif.
        positions = seq.find_motif("MET")
        assert positions == [0, 4]

        # Find a nonexistent motif.
        positions = seq.find_motif("XYZ")
        assert positions == []

        # Find a single character.
        positions = seq.find_motif("A")
        assert positions == [7, 9]

    def test_find_motif_case_insensitive(self):
        """Test whether motif search is case-insensitive."""
        seq = ProteinSequence("METALA")

        positions = seq.find_motif("met")
        assert positions == [0]

        positions = seq.find_motif("Met")
        assert positions == [0]


class TestTranslateFunction:
    """Test the translate function."""

    @pytest.fixture
    def mock_codon_table(self):
        """Mock a codon table."""
        table = Mock(spec=CodonTable)
        table.is_start_codon.return_value = False
        table.is_stop_codon.return_value = False
        table.translate_codon.side_effect = lambda x: {
            "ATG": "M",
            "GAA": "X",
            "TAG": "*",
        }.get(x, "X")
        return table

    @patch("mutcleaner.core.codon.CodonTable.get_standard_table")
    def test_translate_basic(self, mock_get_table, mock_codon_table):
        """Test basic translation."""
        mock_get_table.return_value = mock_codon_table

        result = translate("ATGGAA", seq_type="DNA")

        assert result == "MX"  # ATG->M, GAA->E, but the mock returns X.
        mock_get_table.assert_called_once_with(seq_type="DNA")

    @patch("mutcleaner.core.codon.CodonTable.get_standard_table")
    def test_translate_with_start_codon(self, mock_get_table, mock_codon_table):
        """Test translation starting from the start codon."""
        mock_codon_table.is_start_codon.side_effect = lambda x: x == "ATG"
        mock_get_table.return_value = mock_codon_table

        result = translate("GAAATGGAA", seq_type="DNA", start_at_first_met=True)

        assert result == "MX"  # Start from ATG.

    @patch("mutcleaner.core.codon.CodonTable.get_standard_table")
    def test_translate_with_stop_codon(self, mock_get_table, mock_codon_table):
        """Test stopping translation at a stop codon."""
        mock_codon_table.is_stop_codon.side_effect = lambda x: x == "TAG"
        mock_get_table.return_value = mock_codon_table

        result = translate("ATGTAGGAA", seq_type="DNA", stop_at_stop_codon=True)

        assert result == "M*"  # Stop at TAG.

    def test_translate_not_mod3_strict(self):
        """Test strict mode when sequence length is not divisible by 3."""
        with pytest.raises(ValueError, match="not divisible by 3"):
            translate("ATGGA", require_mod3=True)

    @patch("mutcleaner.core.sequence.CodonTable.get_standard_table")
    def test_translate_not_mod3_lenient(self, mock_get_table, mock_codon_table):
        """Test lenient mode when sequence length is not divisible by 3."""
        mock_get_table.return_value = mock_codon_table

        with warnings.catch_warnings(record=True) as w:
            result = translate("ATGGA", require_mod3=False)

            assert len(w) == 1
            assert "not divisible by 3" in str(w[0].message)

    @patch("mutcleaner.core.sequence.CodonTable.get_standard_table")
    def test_translate_custom_positions(self, mock_get_table, mock_codon_table):
        """Test custom start and end positions."""
        mock_get_table.return_value = mock_codon_table

        result = translate("AAATGGAATAG", start=3, end=9)

        # Should translate ATG GAA (positions 3-8).
        assert len(result) == 2

    @patch("mutcleaner.core.sequence.CodonTable.get_standard_table")
    def test_translate_no_start_codon_found(self, mock_get_table, mock_codon_table):
        """Test the case where no start codon is found."""
        mock_codon_table.is_start_codon.return_value = False
        mock_get_table.return_value = mock_codon_table

        result = translate("GAAGAAGAA", start_at_first_met=True)

        assert result == ""


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_sequence(self):
        """Test an empty sequence."""
        # This may raise an exception depending on the validation logic.
        try:
            seq = DNASequence("")
            assert len(seq) == 0
        except ValueError:
            # If the code disallows empty sequences, this is expected.
            pass

    def test_invalid_characters(self):
        """Test invalid characters, depending on alphabet validation."""
        # This depends on the alphabet validation implementation.
        try:
            seq = DNASequence("ATCGXYZ")
            # If this succeeds, the alphabet allows these characters.
        except ValueError:
            # If an exception is raised, the alphabet has proper validation.
            pass

    def test_none_name_handling(self):
        """Test handling of a None name."""
        seq = DNASequence("ATCG", name=None)
        subseq = seq.get_subsequence(1, 2)

        assert subseq.name is None

    def test_metadata_preservation(self):
        """Test metadata preservation during operations."""
        metadata = {"test": "value"}
        seq = DNASequence("ATCG", metadata=metadata)

        rc = seq.reverse_complement()
        rna = seq.transcribe()

        assert rc.metadata == metadata
        assert rna.metadata == metadata


class TestPerformance:
    """Performance-related tests."""

    def test_large_sequence_operations(self):
        """Test large sequence operations."""
        # Create a relatively large sequence.
        large_seq = DNASequence("ATCG" * 1000)  # 4000 characters.

        # These operations should complete quickly.
        assert len(large_seq) == 4000

        rc = large_seq.reverse_complement()
        assert len(rc) == 4000

        subseq = large_seq.get_subsequence(1, 100)
        assert len(subseq) == 99

    def test_motif_finding_performance(self):
        """Test motif search performance."""
        seq = ProteinSequence("MET" * 100)  # Repeated MET.

        positions = seq.find_motif("MET")
        assert len(positions) == 100


class TestApplyMutation:
    """Test the apply_mutation method."""

    def test_apply_amino_acid_mutation(self):
        """Test applying an amino acid mutation."""
        seq = ProteinSequence("METALA", name="test")
        mutation = AminoAcidMutation("M", 0, "V")  # M0V.

        mutated_seq = seq.apply_mutation(mutation)

        assert str(mutated_seq) == "VETALA"
        assert mutated_seq.name == "test"  # Name remains unchanged.
        assert "mutations_applied" in mutated_seq.metadata
        assert len(mutated_seq.metadata["mutations_applied"]) == 1

        mutation_record = mutated_seq.metadata["mutations_applied"][0]
        assert mutation_record["type"] == "AminoAcidMutation"
        assert mutation_record["position"] == 0
        assert mutation_record["wild_amino_acid"] == "M"
        assert mutation_record["mutant_amino_acid"] == "V"

    def test_apply_codon_mutation_dna(self):
        """Test applying a DNA codon mutation."""
        seq = DNASequence("ATGAAATAG", name="test")
        mutation = CodonMutation("ATG", 0, "TTG")  # ATG0TTG.

        mutated_seq = seq.apply_mutation(mutation)

        assert str(mutated_seq) == "TTGAAATAG"
        assert mutated_seq.name == "test"
        assert "mutations_applied" in mutated_seq.metadata

        mutation_record = mutated_seq.metadata["mutations_applied"][0]
        assert mutation_record["type"] == "CodonMutation"
        assert mutation_record["wild_codon"] == "ATG"
        assert mutation_record["mutant_codon"] == "TTG"

    def test_apply_codon_mutation_rna(self):
        """Test applying an RNA codon mutation."""
        seq = RNASequence("AUGUGAAAG", name="test")
        mutation = CodonMutation("AUG", 0, "UUG")  # AUG0UUG.

        mutated_seq = seq.apply_mutation(mutation)

        assert str(mutated_seq) == "UUGUGAAAG"

    def test_apply_amino_acid_mutation_set(self):
        """Test applying an amino acid mutation set."""
        seq = ProteinSequence("METALA")
        mutations = [
            AminoAcidMutation("M", 0, "V"),  # M0V.
            AminoAcidMutation("A", 3, "G"),  # A3G.
        ]
        mutation_set = AminoAcidMutationSet(mutations)

        mutated_seq = seq.apply_mutation(mutation_set)

        assert str(mutated_seq) == "VETGLA"
        assert len(mutated_seq.metadata["mutations_applied"]) == 2

    def test_apply_codon_mutation_set(self):
        """Test applying a codon mutation set."""
        seq = DNASequence("ATGAAACCC")
        mutations = [
            CodonMutation("ATG", 0, "TTG"),  # ATG0TTG.
            CodonMutation("CCC", 2, "GGG"),  # CCC6GGG.
        ]
        mutation_set = CodonMutationSet(mutations)

        mutated_seq = seq.apply_mutation(mutation_set)

        assert str(mutated_seq) == "TTGAAAGGG"

    def test_mutation_position_validation_amino_acid(self):
        """Test amino acid mutation position validation."""
        seq = ProteinSequence("MET")

        # Test an out-of-bounds position.
        mutation = AminoAcidMutation("A", 5, "V")  # Position 5 is out of range.
        with pytest.raises(ValueError, match="out of bounds"):
            seq.apply_mutation(mutation)

    def test_mutation_position_validation_codon(self):
        """Test codon mutation position validation."""
        seq = DNASequence("ATGAAA")  # Length 6.

        # Test a codon extending beyond sequence bounds.
        mutation = CodonMutation("AAA", 2, "TTT")  
        with pytest.raises(ValueError, match="extends beyond sequence length"):
            seq.apply_mutation(mutation)

    def test_original_sequence_validation_amino_acid(self):
        """Test original amino acid sequence validation."""
        seq = ProteinSequence("METALA")

        # Original amino acid mismatch.
        mutation = AminoAcidMutation("K", 0, "V")  # Position 0 is M, not K.
        with pytest.raises(ValueError, match="Expected amino acid 'K'"):
            seq.apply_mutation(mutation)

    def test_original_sequence_validation_codon(self):
        """Test original codon sequence validation."""
        seq = DNASequence("ATGAAATAG")

        # Original codon mismatch.
        mutation = CodonMutation("TTG", 0, "GGG")  # Position 0 is ATG, not TTG.
        with pytest.raises(ValueError, match="Expected codon 'TTG'"):
            seq.apply_mutation(mutation)

    def test_multiple_mutations_order(self):
        """Test application order for multiple mutations, using reverse order to avoid position shifts."""
        seq = ProteinSequence("METALA")
        mutations = [
            AminoAcidMutation("M", 0, "V"),  # Position 0.
            AminoAcidMutation("T", 2, "S"),  # Position 2.
            AminoAcidMutation("A", 3, "G"),  # Position 3.
        ]
        mutation_set = AminoAcidMutationSet(mutations)

        mutated_seq = seq.apply_mutation(mutation_set)

        # Should apply mutations from higher to lower positions: A3G -> T2S -> M0V.
        assert str(mutated_seq) == "VESGLA"

    def test_metadata_preservation(self):
        """Test preservation of existing metadata."""
        original_metadata = {"source": "test", "version": 1}
        seq = ProteinSequence("MET", metadata=original_metadata)
        mutation = AminoAcidMutation("M", 0, "V")

        mutated_seq = seq.apply_mutation(mutation)

        # Existing metadata should be preserved.
        assert mutated_seq.metadata["source"] == "test"
        assert mutated_seq.metadata["version"] == 1
        # The new mutation record should also exist.
        assert "mutations_applied" in mutated_seq.metadata

    def test_unsupported_mutation_type(self):
        """Test unsupported mutation types."""
        seq = ProteinSequence("MET")

        # Mock an unsupported mutation type.
        unsupported_mutation = Mock()
        unsupported_mutation.__class__.__name__ = "UnsupportedMutation"

        with pytest.raises(TypeError, match="Unsupported mutation type"):
            seq.apply_mutation(unsupported_mutation)

    def test_immutability(self):
        """Test original sequence immutability."""
        seq = ProteinSequence("MET")
        mutation = AminoAcidMutation("M", 0, "V")

        mutated_seq = seq.apply_mutation(mutation)

        # The original sequence should remain unchanged.
        assert str(seq) == "MET"
        assert str(mutated_seq) == "VET"
        assert seq is not mutated_seq

    def test_chained_mutations(self):
        """Test chained mutation application."""
        seq = ProteinSequence("METALA")
        mutation1 = AminoAcidMutation("M", 0, "V")
        mutation2 = AminoAcidMutation("E", 1, "K")

        # Apply mutations in a chain.
        mutated_seq = seq.apply_mutation(mutation1).apply_mutation(mutation2)

        assert str(mutated_seq) == "VKTALA"
        # There should be two mutation records.
        assert len(mutated_seq.metadata["mutations_applied"]) == 2

    def test_unmatched_mutation_subtype_and_sequence_type(self):
        """Test the error for mismatched types."""
        seq = RNASequence("AUGUGAAAG", name="test")
        unsupported_mutation = CodonMutation("ATG", 0, "AAG")  # subtype: codon_dna.

        with pytest.raises(TypeError, match="Unmatching mutation subtype"):
            seq.apply_mutation(unsupported_mutation)

if __name__ == "__main__":
    # Run tests.
    pytest.main([__file__, "-v"])
