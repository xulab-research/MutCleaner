import pytest

from mutcleaner.core.mutation import (
    BaseMutation,
    AminoAcidMutation,
    CodonMutation,
    MutationSet,
    AminoAcidMutationSet,
    CodonMutationSet,
)
from mutcleaner.core.alphabet import ProteinAlphabet
from mutcleaner.core.codon import CodonTable


class TestBaseMutation:
    """Test the base Mutation class."""

    def test_base_mutation_is_abstract(self):
        """Test that BaseMutation is an abstract class."""
        with pytest.raises(TypeError):
            BaseMutation(position=1)  # type: ignore

    def test_invalid_position(self):
        """Test invalid mutation positions."""
        with pytest.raises(ValueError, match="Position must be non-negative"):
            AminoAcidMutation("A", -1, "V")


class TestAminoAcidMutation:
    """Test the amino acid mutation class."""

    def test_amino_acid_mutation_creation(self):
        """Test amino acid mutation creation."""
        mutation = AminoAcidMutation("A", 123, "V")

        assert mutation.wild_amino_acid == "A"
        assert mutation.position == 123
        assert mutation.mutant_amino_acid == "V"
        assert mutation.type == "amino_acid"
        assert str(mutation) == "A123V"

    def test_amino_acid_mutation_with_metadata(self):
        """Test amino acid mutation creation with metadata."""
        metadata = {"source": "test", "confidence": 0.95}
        mutation = AminoAcidMutation("A", 123, "V", metadata=metadata)

        assert mutation.metadata == metadata

    def test_amino_acid_mutation_with_custom_alphabet(self):
        """Test amino acid mutation creation with a custom alphabet."""
        alphabet = ProteinAlphabet(include_stop=False)
        mutation = AminoAcidMutation("A", 123, "V", alphabet=alphabet)

        assert mutation.alphabet == alphabet

    def test_amino_acid_mutation_case_insensitive(self):
        """Test that amino acid mutations are case-insensitive."""
        mutation = AminoAcidMutation("a", 123, "v")

        assert mutation.wild_amino_acid == "A"
        assert mutation.mutant_amino_acid == "V"

    def test_amino_acid_mutation_validation(self):
        """Test amino acid mutation validation."""
        # Test invalid amino acids.
        with pytest.raises(ValueError, match="Invalid amino acid mutation"):
            AminoAcidMutation("X", 123, "V")

        # Test invalid positions.
        with pytest.raises(ValueError, match="Position must be non-negative"):
            AminoAcidMutation("A", -1, "V")

    def test_amino_acid_mutation_types(self):
        """Test amino acid mutation type classification."""
        # Synonymous mutation.
        synonymous = AminoAcidMutation("A", 123, "A")
        assert synonymous.is_synonymous()
        assert not synonymous.is_missense()
        assert not synonymous.is_nonsense()
        assert synonymous.get_mutation_category() == "synonymous"
        assert synonymous.effect_type == "synonymous"

        # Missense mutation.
        missense = AminoAcidMutation("A", 123, "V")
        assert not missense.is_synonymous()
        assert missense.is_missense()
        assert not missense.is_nonsense()
        assert missense.get_mutation_category() == "missense"

        # Nonsense mutation.
        nonsense = AminoAcidMutation("A", 123, "*")
        assert not nonsense.is_synonymous()
        assert not nonsense.is_missense()
        assert nonsense.is_nonsense()
        assert nonsense.get_mutation_category() == "nonsense"

    def test_amino_acid_mutation_equality(self):
        """Test amino acid mutation equality."""
        mutation1 = AminoAcidMutation("A", 123, "V")
        mutation2 = AminoAcidMutation("A", 123, "V")
        mutation3 = AminoAcidMutation("A", 124, "V")
        mutation4 = AminoAcidMutation("A", 123, "L")

        assert mutation1 == mutation2
        assert mutation1 != mutation3
        assert mutation1 != mutation4

    def test_amino_acid_mutation_hash(self):
        """Test that amino acid mutations are hashable."""
        mutation1 = AminoAcidMutation("A", 123, "V")
        mutation2 = AminoAcidMutation("A", 123, "V")
        mutation3 = AminoAcidMutation("A", 124, "V")

        assert hash(mutation1) == hash(mutation2)
        assert hash(mutation1) != hash(mutation3)

        # Test that mutations can be used as set elements.
        mutation_set = {mutation1, mutation2, mutation3}
        assert len(mutation_set) == 2

    def test_amino_acid_mutation_from_string_one_letter(self):
        """Test parsing an amino acid mutation from a one-letter string."""
        mutation = AminoAcidMutation.from_string("A123V")

        assert mutation.wild_amino_acid == "A"
        assert mutation.position == 122
        assert mutation.mutant_amino_acid == "V"

    def test_amino_acid_mutation_from_string_three_letter(self):
        """Test parsing an amino acid mutation from a three-letter string."""
        mutation = AminoAcidMutation.from_string("Ala123Val")

        assert mutation.wild_amino_acid == "A"
        assert mutation.position == 122
        assert mutation.mutant_amino_acid == "V"

    def test_amino_acid_mutation_from_string_with_stop(self):
        """Test parsing an amino acid mutation containing a stop codon."""
        mutation = AminoAcidMutation.from_string("A123*")

        assert mutation.wild_amino_acid == "A"
        assert mutation.position == 122
        assert mutation.mutant_amino_acid == "*"
        assert mutation.is_nonsense()

    def test_amino_acid_mutation_from_string_invalid_format(self):
        """Test parsing invalid amino acid mutation formats."""
        invalid_formats = [
            "A123",  # Missing the mutant amino acid.
            "123V",  # Missing the wild-type amino acid.
            "ABC123V",  # Invalid format.
            "A-123V",  # Invalid character.
            "",  # Empty string.
            "   ",  # Whitespace-only string.
        ]

        for invalid_format in invalid_formats:
            with pytest.raises(ValueError, match="Invalid mutation format"):
                AminoAcidMutation.from_string(invalid_format)

    def test_amino_acid_mutation_from_string_invalid_amino_acid(self):
        """Test parsing invalid amino acid codes."""
        with pytest.raises(ValueError, match="Unknown three-letter amino acid code"):
            AminoAcidMutation.from_string("Xxx123Val")


class TestCodonMutation:
    """Test the codon mutation class."""

    def test_codon_mutation_creation_dna(self):
        """Test DNA codon mutation creation."""
        mutation = CodonMutation("ATG", 1, "TAA")

        assert mutation.wild_codon == "ATG"
        assert mutation.position == 1
        assert mutation.mutant_codon == "TAA"
        assert mutation.seq_type == "DNA"
        assert mutation.type == "codon_dna"
        assert str(mutation) == "ATG1TAA"

    def test_codon_mutation_creation_rna(self):
        """Test RNA codon mutation creation."""
        mutation = CodonMutation("AUG", 1, "UAA")

        assert mutation.wild_codon == "AUG"
        assert mutation.position == 1
        assert mutation.mutant_codon == "UAA"
        assert mutation.seq_type == "RNA"
        assert mutation.type == "codon_rna"

    def test_codon_mutation_creation_both(self):
        """Test codon mutation creation without T or U bases."""
        mutation = CodonMutation("ACG", 1, "GCA")

        assert mutation.seq_type == "Both"
        assert mutation.type == "codon_both"

    def test_codon_mutation_mixed_tu_error(self):
        """Test the error raised for mixed T and U bases."""
        with pytest.raises(ValueError, match="Codons cannot contain both T and U"):
            CodonMutation("ATG", 1, "UAA")

    def test_codon_mutation_case_insensitive(self):
        """Test that codon mutations are case-insensitive."""
        mutation = CodonMutation("atg", 1, "taa")

        assert mutation.wild_codon == "ATG"
        assert mutation.mutant_codon == "TAA"

    def test_codon_mutation_validation(self):
        """Test codon mutation validation."""
        # Test invalid codon lengths.
        with pytest.raises(ValueError, match="Invalid codon mutation"):
            CodonMutation("AT", 1, "TAA")

        with pytest.raises(ValueError, match="Invalid codon mutation"):
            CodonMutation("ATG", 1, "TA")

        # Test invalid characters.
        with pytest.raises(ValueError, match="Invalid codon mutation"):
            CodonMutation("ATX", 1, "TAA")

        # Test invalid positions.
        with pytest.raises(ValueError, match="Position must be non-negative"):
            CodonMutation("ATG", -1, "TAA")

        with pytest.raises(ValueError, match="Invalid codon mutation"):
            CodonMutation("ATG", 1.4, "TAA")  # type: ignore

    def test_codon_mutation_to_amino_acid(self):
        """Test converting a codon mutation to an amino acid mutation."""
        codon_mutation = CodonMutation("ATG", 1, "TAA")  # Met -> Stop
        aa_mutation = codon_mutation.to_amino_acid_mutation()

        assert aa_mutation.wild_amino_acid == "M"
        assert aa_mutation.position == 1
        assert aa_mutation.mutant_amino_acid == "*"
        assert aa_mutation.is_nonsense()

    def test_codon_mutation_to_amino_acid_with_custom_table(self):
        """Test conversion using a custom codon table."""
        codon_table = CodonTable.get_standard_table("DNA")
        codon_mutation = CodonMutation("ATG", 1, "TTG")  # Met -> Leu
        aa_mutation = codon_mutation.to_amino_acid_mutation(codon_table)

        assert aa_mutation.wild_amino_acid == "M"
        assert aa_mutation.mutant_amino_acid == "L"

    def test_codon_mutation_from_string(self):
        """Test parsing a codon mutation from a string."""
        mutation = CodonMutation.from_string("ATG123TAA")

        assert mutation.wild_codon == "ATG"
        assert mutation.position == 122
        assert mutation.mutant_codon == "TAA"

    def test_codon_mutation_from_string_rna(self):
        """Test parsing a codon mutation from an RNA string."""
        mutation = CodonMutation.from_string("AUG123UAA")

        assert mutation.wild_codon == "AUG"
        assert mutation.position == 122
        assert mutation.mutant_codon == "UAA"
        assert mutation.seq_type == "RNA"

    def test_codon_mutation_from_string_invalid_format(self):
        """Test parsing invalid codon mutation formats."""
        invalid_formats = [
            "ATG123",  # Missing the mutant codon.
            "123TAA",  # Missing the wild-type codon.
            "AT123TAA",  # Invalid codon length.
            "ATG123TAAA",  # Invalid codon length.
            "",  # Empty string.
            "ATX123TAA",  # Invalid character.
        ]

        for invalid_format in invalid_formats:
            with pytest.raises(ValueError, match="Invalid codon mutation format"):
                CodonMutation.from_string(invalid_format)


class TestMutationSet:
    """Test the mutation set class."""

    def test_mutation_set_creation(self):
        """Test mutation set creation."""
        mutations = [
            AminoAcidMutation("A", 1, "V"),
            AminoAcidMutation("L", 2, "P"),
        ]
        mutation_set = MutationSet(mutations, AminoAcidMutation)

        assert len(mutation_set) == 2
        assert mutation_set.mutation_type == AminoAcidMutation
        assert mutation_set.mutation_subtype == "amino_acid"

    def test_mutation_set_auto_type_inference(self):
        """Test automatic mutation type inference."""
        mutations = [AminoAcidMutation("A", 1, "V")]
        mutation_set = MutationSet(mutations, None)

        assert mutation_set.mutation_type == AminoAcidMutation

    def test_mutation_set_empty_error(self):
        """Test the error raised for an empty mutation set."""
        with pytest.raises(
            ValueError, match="MutationSet must contain at least one mutation"
        ):
            MutationSet([], AminoAcidMutation)

    def test_mutation_set_type_consistency(self):
        """Test mutation type consistency."""
        mutations = [
            AminoAcidMutation("A", 1, "V"),
            CodonMutation("ATG", 2, "TAA"),  # Different type.
        ]

        with pytest.raises(ValueError, match="All mutations must be of type"):
            MutationSet(mutations, AminoAcidMutation)

    def test_mutation_set_subtype_consistency(self):
        """Test mutation subtype consistency."""
        mutations = [
            CodonMutation("ATG", 1, "TAA"),  # DNA codon.
            CodonMutation("AUG", 2, "UAA"),  # RNA codon.
        ]

        with pytest.raises(
            ValueError, match="Codon mutations contain incompatible DNA and RNA types"
        ):
            MutationSet(mutations, CodonMutation)

    def test_mutation_set_duplicate_positions(self):
        """Test the error raised for duplicate positions."""
        mutations = [
            AminoAcidMutation("A", 1, "V"),
            AminoAcidMutation("L", 1, "P"),  # Same position.
        ]

        with pytest.raises(ValueError, match="Duplicate mutations at positions"):
            MutationSet(mutations, AminoAcidMutation)

    def test_mutation_set_sorting(self):
        """Test sorting mutations by position."""
        mutations = [
            AminoAcidMutation("L", 3, "P"),
            AminoAcidMutation("A", 1, "V"),
            AminoAcidMutation("G", 2, "S"),
        ]
        mutation_set = MutationSet(mutations, AminoAcidMutation)

        positions = [m.position for m in mutation_set.mutations]
        assert positions == [1, 2, 3]

    def test_mutation_set_iteration(self):
        """Test mutation set iteration."""
        mutations = [
            AminoAcidMutation("A", 1, "V"),
            AminoAcidMutation("L", 2, "P"),
        ]
        mutation_set = MutationSet(mutations, AminoAcidMutation)

        mutation_list = list(mutation_set)
        assert len(mutation_list) == 2
        assert mutation_list[0].position == 1
        assert mutation_list[1].position == 2

    def test_mutation_set_add_mutation(self):
        """Test adding a mutation."""
        mutations = [AminoAcidMutation("A", 1, "V")]
        mutation_set = MutationSet(mutations, AminoAcidMutation)

        new_mutation = AminoAcidMutation("L", 2, "P")
        mutation_set.add_mutation(new_mutation)

        assert len(mutation_set) == 2
        assert mutation_set.get_mutation_at(2) == new_mutation

    def test_mutation_set_add_wrong_type(self):
        """Test adding a mutation with the wrong type."""
        mutations = [AminoAcidMutation("A", 1, "V")]
        mutation_set = MutationSet(mutations, AminoAcidMutation)

        codon_mutation = CodonMutation("ATG", 2, "TAA")
        with pytest.raises(ValueError, match="Mutation must be of type"):
            mutation_set.add_mutation(codon_mutation)  # type: ignore

    def test_mutation_set_add_duplicate_position(self):
        """Test adding a mutation with a duplicate position."""
        mutations = [AminoAcidMutation("A", 1, "V")]
        mutation_set = MutationSet(mutations, AminoAcidMutation)

        duplicate_mutation = AminoAcidMutation("A", 1, "L")
        with pytest.raises(ValueError, match="Mutation already exists at position"):
            mutation_set.add_mutation(duplicate_mutation)

    def test_mutation_set_remove_mutation(self):
        """Test removing a mutation."""
        mutations = [
            AminoAcidMutation("A", 1, "V"),
            AminoAcidMutation("L", 2, "P"),
        ]
        mutation_set = MutationSet(mutations, AminoAcidMutation)

        removed = mutation_set.remove_mutation(1)
        assert removed is True
        assert len(mutation_set) == 1
        assert not mutation_set.has_mutation_at(1)

        # Remove a nonexistent position.
        removed = mutation_set.remove_mutation(999)
        assert removed is False

    def test_mutation_set_get_mutation_at(self):
        """Test getting the mutation at a specified position."""
        mutation = AminoAcidMutation("A", 1, "V")
        mutations = [mutation]
        mutation_set = MutationSet(mutations, AminoAcidMutation)

        assert mutation_set.get_mutation_at(1) == mutation
        assert mutation_set.get_mutation_at(999) is None

    def test_mutation_set_positions(self):
        """Test retrieving mutation positions."""
        mutations = [
            AminoAcidMutation("A", 1, "V"),
            AminoAcidMutation("L", 3, "P"),
            AminoAcidMutation("G", 2, "S"),
        ]
        mutation_set = MutationSet(mutations, AminoAcidMutation)

        positions = mutation_set.get_positions()
        assert positions == [1, 2, 3]  # Should be sorted.

        positions_set = mutation_set.get_positions_set()
        assert positions_set == {1, 2, 3}

    def test_mutation_set_validation(self):
        """Test mutation set validation."""
        mutations = [AminoAcidMutation("A", 1, "V")]
        mutation_set = MutationSet(mutations, AminoAcidMutation)

        assert mutation_set.validate_all() is True

    def test_mutation_set_categories(self):
        """Test mutation category counts."""
        mutations = [
            AminoAcidMutation("A", 1, "A"),  # Synonymous.
            AminoAcidMutation("L", 2, "P"),  # Missense.
            AminoAcidMutation("G", 3, "*"),  # Nonsense.
        ]
        mutation_set = MutationSet(mutations, AminoAcidMutation)

        categories = mutation_set.get_mutation_categories()
        assert categories["synonymous"] == 1
        assert categories["missense"] == 1
        assert categories["nonsense"] == 1

    def test_mutation_set_filter_by_category(self):
        """Test filtering mutations by category."""
        mutations = [
            AminoAcidMutation("A", 1, "A"),  # Synonymous.
            AminoAcidMutation("L", 2, "P"),  # Missense.
            AminoAcidMutation("G", 3, "*"),  # Nonsense.
        ]
        mutation_set = MutationSet(mutations, AminoAcidMutation)

        missense_mutations = mutation_set.filter_by_category("missense")
        assert len(missense_mutations) == 1
        assert missense_mutations[0].position == 2

    def test_mutation_set_from_string_single(self):
        """Test creating a single-mutation set from a string."""
        mutation_set = MutationSet.from_string("A123V")

        assert len(mutation_set) == 1
        assert mutation_set.mutation_type == AminoAcidMutation
        assert mutation_set.mutations[0].position == 122

    def test_mutation_set_from_string_multiple_comma(self):
        """Test creating a mutation set from a comma-separated string."""
        mutation_set = MutationSet.from_string("A123V,L456P,G789S")

        assert len(mutation_set) == 3
        assert mutation_set.mutation_type == AminoAcidMutation
        positions = [m.position for m in mutation_set.mutations]
        assert positions == [122, 455, 788]

    def test_mutation_set_from_string_multiple_semicolon(self):
        """Test creating a mutation set from a semicolon-separated string."""
        mutation_set = MutationSet.from_string("A123V;L456P;G789S")

        assert len(mutation_set) == 3
        positions = [m.position for m in mutation_set.mutations]
        assert positions == [122, 455, 788]

    def test_mutation_set_from_string_with_spaces(self):
        """Test strings containing spaces."""
        mutation_set = MutationSet.from_string(" A123V , L456P , G789S ")

        assert len(mutation_set) == 3

    def test_mutation_set_from_string_custom_separator(self):
        """Test custom separators."""
        mutation_set = MutationSet.from_string("A123V|L456P|G789S", sep="|")

        assert len(mutation_set) == 3

    def test_mutation_set_from_string_specified_type(self):
        """Test specifying a mutation type."""
        mutation_set = MutationSet.from_string("ATG123TAA", mutation_type=CodonMutation)

        assert len(mutation_set) == 1
        assert mutation_set.mutation_type == CodonMutation

    def test_mutation_set_from_string_with_name_metadata(self):
        """Test string-based creation with a name and metadata."""
        metadata = {"source": "test"}
        mutation_set = MutationSet.from_string(
            "A123V", name="test_set", metadata=metadata
        )

        assert mutation_set.name == "test_set"
        assert mutation_set.metadata == metadata

    def test_mutation_set_from_string_empty_error(self):
        """Test the error raised for an empty input string."""
        with pytest.raises(ValueError, match="Input string cannot be empty"):
            MutationSet.from_string("")

        with pytest.raises(ValueError, match="Input string cannot be empty"):
            MutationSet.from_string("   ")

    def test_mutation_set_from_string_parse_error(self):
        """Test parsing errors."""
        with pytest.raises(ValueError, match="Some mutations could not be parsed"):
            MutationSet.from_string("A123V,INVALID,L456P")

    def test_mutation_set_from_string_no_valid_mutations(self):
        """Test input strings with no valid mutations."""
        with pytest.raises(ValueError, match="No valid mutations could be parsed"):
            MutationSet.from_string("INVALID1,INVALID2")

    def test_mutation_set_str_repr(self):
        """Test string representations."""
        mutations = [AminoAcidMutation("A", 1, "V")]

        # Without a name.
        mutation_set = MutationSet(mutations, AminoAcidMutation)
        assert str(mutation_set) in "MutationSet: A1V"

        # With a name.
        named_set = MutationSet(mutations, AminoAcidMutation, name="test")
        assert str(named_set) in "MutationSet(test): A1V"

        # repr.
        repr_str = repr(mutation_set)
        assert "MutationSet(mutations=" in repr_str


class TestAminoAcidMutationSet:
    """Test the amino acid mutation set class."""

    def test_amino_acid_mutation_set_creation(self):
        """Test amino acid mutation set creation."""
        mutations = [
            AminoAcidMutation("A", 1, "V"),
            AminoAcidMutation("L", 2, "P"),
        ]
        aa_set = AminoAcidMutationSet(mutations)

        assert len(aa_set) == 2
        assert aa_set.mutation_type == AminoAcidMutation

    def test_amino_acid_mutation_set_effect_filters(self):
        """Test filtering by effect type."""
        mutations = [
            AminoAcidMutation("A", 1, "A"),  # Synonymous.
            AminoAcidMutation("L", 2, "P"),  # Missense.
            AminoAcidMutation("G", 3, "*"),  # Nonsense.
            AminoAcidMutation("T", 4, "S"),  # Missense.
        ]
        aa_set = AminoAcidMutationSet(mutations)

        synonymous = aa_set.get_synonymous_mutations()
        assert len(synonymous) == 1
        assert synonymous[0].position == 1

        missense = aa_set.get_missense_mutations()
        assert len(missense) == 2
        assert {m.position for m in missense} == {2, 4}

        nonsense = aa_set.get_nonsense_mutations()
        assert len(nonsense) == 1
        assert nonsense[0].position == 3

    def test_amino_acid_mutation_set_stop_codon_check(self):
        """Test stop-codon mutation checks."""
        mutations_with_stop = [
            AminoAcidMutation("A", 1, "V"),
            AminoAcidMutation("L", 2, "*"),
        ]
        aa_set_with_stop = AminoAcidMutationSet(mutations_with_stop)
        assert aa_set_with_stop.has_stop_codon_mutations() is True

        mutations_without_stop = [
            AminoAcidMutation("A", 1, "V"),
            AminoAcidMutation("L", 2, "P"),
        ]
        aa_set_without_stop = AminoAcidMutationSet(mutations_without_stop)
        assert aa_set_without_stop.has_stop_codon_mutations() is False

    def test_amino_acid_mutation_set_count_by_effect(self):
        """Test counting mutations by effect type."""
        mutations = [
            AminoAcidMutation("A", 1, "A"),  # Synonymous.
            AminoAcidMutation("L", 2, "P"),  # Missense.
            AminoAcidMutation("G", 3, "*"),  # Nonsense.
            AminoAcidMutation("T", 4, "S"),  # Missense.
        ]
        aa_set = AminoAcidMutationSet(mutations)

        counts = aa_set.count_by_effect_type()
        assert counts["synonymous"] == 1
        assert counts["missense"] == 2
        assert counts["nonsense"] == 1


class TestCodonMutationSet:
    """Test the codon mutation set class."""

    def test_codon_mutation_set_creation(self):
        """Test codon mutation set creation."""
        mutations = [
            CodonMutation("ATG", 1, "TAA"),
            CodonMutation("TTG", 2, "TAG"),
        ]
        codon_set = CodonMutationSet(mutations)

        assert len(codon_set) == 2
        assert codon_set.mutation_type == CodonMutation
        assert codon_set.seq_type == "DNA"

    def test_codon_mutation_set_seq_type(self):
        """Test the sequence type property."""
        # DNA mutations.
        dna_mutations = [CodonMutation("ATG", 1, "TAA")]
        dna_set = CodonMutationSet(dna_mutations)
        assert dna_set.seq_type == "DNA"

        # RNA mutations.
        rna_mutations = [CodonMutation("AUG", 1, "UAA")]
        rna_set = CodonMutationSet(rna_mutations)
        assert rna_set.seq_type == "RNA"

        # Both mutations.
        both_mutations = [CodonMutation("ACG", 1, "GCA")]
        both_set = CodonMutationSet(both_mutations)
        assert both_set.seq_type == "Both"

    def test_codon_mutation_set_to_amino_acid(self):
        """Test converting to an amino acid mutation set."""
        mutations = [
            CodonMutation("ATG", 1, "TAA"),  # Met -> Stop
            CodonMutation("TTG", 2, "CTG"),  # Leu -> Leu
        ]
        codon_set = CodonMutationSet(mutations, name="test_codons")

        aa_set = codon_set.to_amino_acid_mutation_set()

        assert isinstance(aa_set, AminoAcidMutationSet)
        assert len(aa_set) == 2
        assert aa_set.name == "test_codons_aa"

        # Check the conversion results.
        aa_mutations = list(aa_set)
        assert aa_mutations[0].wild_amino_acid == "M"
        assert aa_mutations[0].mutant_amino_acid == "*"
        assert aa_mutations[1].wild_amino_acid == "L"
        assert aa_mutations[1].mutant_amino_acid == "L"

    def test_codon_mutation_set_to_amino_acid_custom_table(self):
        """Test conversion using a custom codon table."""
        mutations = [CodonMutation("ATG", 1, "TTG")]
        codon_set = CodonMutationSet(mutations)

        codon_table = CodonTable.get_standard_table("DNA")
        aa_set = codon_set.to_amino_acid_mutation_set(codon_table)

        aa_mutation = list(aa_set)[0]
        assert aa_mutation.wild_amino_acid == "M"
        assert aa_mutation.mutant_amino_acid == "L"


class TestMutationInference:
    """Test mutation type inference."""

    def test_infer_amino_acid_mutation(self):
        """Test inferring an amino acid mutation."""
        mutation = MutationSet._infer_and_create_mutation("A123V")
        assert isinstance(mutation, AminoAcidMutation)

    def test_infer_codon_mutation(self):
        """Test inferring a codon mutation."""
        mutation = MutationSet._infer_and_create_mutation("ATG123TAA")
        assert isinstance(mutation, CodonMutation)

    def test_infer_unknown_mutation(self):
        """Test inferring an unknown mutation format."""
        with pytest.raises(ValueError, match="Could not parse mutation string"):
            MutationSet._infer_and_create_mutation("INVALID_FORMAT")


class TestSeparatorGuessing:
    """Test separator guessing."""

    def test_guess_separator_comma(self):
        """Test guessing a comma separator."""
        sep = MutationSet._guess_sep("A123V,L456P,G789S")
        assert sep == ","

    def test_guess_separator_semicolon(self):
        """Test guessing a semicolon separator."""
        sep = MutationSet._guess_sep("A123V;L456P;G789S")
        assert sep == ";"

    def test_guess_separator_pipe(self):
        """Test guessing a pipe separator."""
        sep = MutationSet._guess_sep("A123V|L456P|G789S")
        assert sep == "|"

    def test_guess_separator_none(self):
        """Test the no-separator case."""
        sep = MutationSet._guess_sep("A123V")
        assert sep is None

    def test_guess_separator_empty(self):
        """Test an empty string."""
        sep = MutationSet._guess_sep("")
        assert sep is None

    def test_guess_separator_priority(self):
        """Test separator priority."""
        # When multiple separators have the same count, choose the higher-priority one.
        sep = MutationSet._guess_sep("A123V;L456P,G789S")  # One semicolon and one comma.
        assert sep == ";"  # Semicolon has higher priority.
