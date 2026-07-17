# mutcleaner/core/dataset.py
from __future__ import annotations

import pickle
import json
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from typing import Any, Dict, TYPE_CHECKING

from .mutation import (
    MutationSet,
    AminoAcidMutationSet,
    CodonMutationSet,
    AminoAcidMutation,
    CodonMutation,
)
from .sequence import (
    DNASequence,
    ProteinSequence,
    RNASequence,
    load_sequences_from_fasta,
)

if TYPE_CHECKING:
    from typing import List, Literal, Optional, Sequence, Type, Union

    from .mutation import BaseMutation
    from .sequence import BaseSequence

__all__ = ["MutationDataset"]


def __dir__() -> List[str]:
    return __all__


SEQUENCE_TYPE_MAP = {
    "ProteinSequence": ProteinSequence,
    "DNASequence": DNASequence,
    "RNASequence": RNASequence,
    "DNA": DNASequence,
    "RNA": RNASequence,
    "Protein": ProteinSequence,
}


class MutationDataset:
    """
    Dataset container for cleaned mutation data with multiple reference sequences.

    All mutation sets must be linked to a reference sequence when added to the dataset.
    This ensures data integrity and enables proper validation and analysis.
    """

    def __init__(self, name: Optional[str] = None):
        """
        Initialize a new MutationDataset.

        Parameters
        ----------
        name : Optional[str], default=None
            Optional name for the dataset

        Note:
            All mutation sets added to this dataset must be linked to a reference sequence.
            Use add_reference_sequence() first, then add_mutation_set() with reference_id.
        """
        self.name = name
        self.reference_sequences: Dict[str, BaseSequence] = (
            {}
        )  # sequence_id -> sequence
        self.mutation_sets: List[MutationSet] = []
        self.mutation_set_references: Dict[int, str] = (
            {}
        )  # mutation_set_index -> sequence_id
        self.mutation_set_labels: Dict[int, Any] = {}  # mutation_set_index -> label
        self.metadata: Dict[str, Any] = {}
        self._df: Optional[pd.DataFrame] = None

    def __len__(self) -> int:
        return len(self.mutation_sets)

    def __iter__(self):
        """
        Iterate over mutation sets and their reference sequence IDs.

        Yields
        ------
        Tuple[MutationSet, str]: (mutation_set, reference_id) pairs

        Examples
        --------
        >>> for mutation_set, ref_id in dataset:
        ...    print(f"Processing {len(mutation_set)} mutations for {ref_id}")
        ...    ref_seq = dataset.get_reference_sequence(ref_id)
        ...    # ... analysis code
        """
        for i, mutation_set in enumerate(self.mutation_sets):
            reference_id = self.mutation_set_references[i]
            yield mutation_set, reference_id

    def __str__(self) -> str:
        stats = self.get_statistics()
        ref_count = stats["num_reference_sequences"]
        ref_info = (
            f" ({ref_count} reference sequences)"
            if ref_count > 0
            else " (no references)"
        )

        return (
            f"MutationDataset({self.name}){ref_info}: "
            f"{stats['total_mutation_sets']} mutation sets, "
            f"{stats['total_mutations']} mutations"
        )

    def add_reference_sequence(self, sequence_id: str, sequence: BaseSequence):
        """Add a reference sequence with a unique identifier"""
        if sequence_id in self.reference_sequences:
            raise ValueError(
                f"Reference sequence with ID '{sequence_id}' already exists"
            )

        self.reference_sequences[sequence_id] = sequence
        self._df = None  # Reset cached DataFrame

    def remove_reference_sequence(self, sequence_id: str):
        """Remove a reference sequence"""
        if sequence_id not in self.reference_sequences:
            raise ValueError(f"Reference sequence with ID '{sequence_id}' not found")

        # Check if any mutation sets reference this sequence
        referencing_sets = [
            idx
            for idx, ref_id in self.mutation_set_references.items()
            if ref_id == sequence_id
        ]
        if referencing_sets:
            raise ValueError(
                f"Cannot remove sequence '{sequence_id}' as it is referenced by "
                f"{len(referencing_sets)} mutation sets. Remove the mutation sets first."
            )

        del self.reference_sequences[sequence_id]
        self._df = None

    def get_reference_sequence(self, sequence_id: str) -> BaseSequence:
        """Get a reference sequence by ID"""
        if sequence_id not in self.reference_sequences:
            raise ValueError(f"Reference sequence with ID '{sequence_id}' not found")
        return self.reference_sequences[sequence_id]

    def list_reference_sequences(self) -> List[str]:
        """Get list of all reference sequence IDs"""
        return list(self.reference_sequences.keys())

    def add_mutation_set(
        self,
        mutation_set: MutationSet,
        reference_id: str,
        label: Optional[float] = None,
    ):
        """Add a mutation set to the dataset, linking to a reference sequence"""
        if reference_id not in self.reference_sequences:
            raise ValueError(f"Reference sequence with ID '{reference_id}' not found")

        mutation_set_index = len(self.mutation_sets)
        self.mutation_sets.append(mutation_set)
        self.mutation_set_references[mutation_set_index] = reference_id
        self.mutation_set_labels[mutation_set_index] = label

        self._df = None  # Reset cached DataFrame

    def add_mutation_sets(
        self,
        mutation_sets: Sequence[MutationSet],
        reference_ids: Sequence[str],
        labels: Optional[Sequence[float]] = None,
    ):
        """Add multiple mutation sets to the dataset"""
        if len(reference_ids) != len(mutation_sets):
            raise ValueError(
                "Number of reference_ids must match number of mutation_sets"
            )

        if labels is not None and len(labels) != len(mutation_sets):
            raise ValueError("Number of labels must match number of mutation_sets")

        for i, (mutation_set, ref_id) in enumerate(zip(mutation_sets, reference_ids)):
            label = labels[i] if labels is not None else None
            self.add_mutation_set(mutation_set, ref_id, label)

    def set_mutation_set_reference(self, mutation_set_index: int, reference_id: str):
        """Set the reference sequence for a specific mutation set"""
        if mutation_set_index >= len(self.mutation_sets):
            raise ValueError(f"Mutation set index {mutation_set_index} out of range")
        if reference_id not in self.reference_sequences:
            raise ValueError(f"Reference sequence with ID '{reference_id}' not found")

        self.mutation_set_references[mutation_set_index] = reference_id
        self._df = None

    def get_mutation_set_reference(self, mutation_set_index: int) -> str:
        """Get the reference sequence ID for a specific mutation set"""
        if mutation_set_index >= len(self.mutation_sets):
            raise ValueError(f"Mutation set index {mutation_set_index} out of range")
        return self.mutation_set_references[mutation_set_index]

    def set_mutation_set_label(self, mutation_set_index: int, label: float):
        """Set the label for a specific mutation set"""
        if mutation_set_index >= len(self.mutation_sets):
            raise ValueError(f"Mutation set index {mutation_set_index} out of range")
        self.mutation_set_labels[mutation_set_index] = label
        self._df = None

    def get_mutation_set_label(self, mutation_set_index: int) -> Any:
        """Get the label for a specific mutation set"""
        if mutation_set_index >= len(self.mutation_sets):
            raise ValueError(f"Mutation set index {mutation_set_index} out of range")
        return self.mutation_set_labels.get(mutation_set_index)

    def remove_mutation_set(self, mutation_set_index: int):
        """Remove a mutation set from the dataset"""
        if mutation_set_index >= len(self.mutation_sets):
            raise ValueError(f"Mutation set index {mutation_set_index} out of range")

        # Remove the mutation set
        del self.mutation_sets[mutation_set_index]

        # Update the reference mapping (shift indices)
        new_references = {}
        new_labels = {}
        for idx, ref_id in self.mutation_set_references.items():
            if idx < mutation_set_index:
                new_references[idx] = ref_id
                new_labels[idx] = self.mutation_set_labels.get(idx)
            elif idx > mutation_set_index:
                new_references[idx - 1] = ref_id
                new_labels[idx - 1] = self.mutation_set_labels.get(idx)
            # Skip the removed index

        self.mutation_set_references = new_references
        self.mutation_set_labels = new_labels
        self._df = None

    def validate_against_references(self) -> Dict[str, Any]:
        """Validate mutations against their reference sequences"""
        validation_results = {
            "valid_mutation_sets": [],
            "invalid_mutation_sets": [],
            "position_mismatches": [],
        }

        for i, mutation_set in enumerate(self.mutation_sets):
            set_name = mutation_set.name or f"MutationSet_{i}"

            # All mutation sets must have a reference sequence
            reference_id = self.mutation_set_references[i]  # This should always exist
            reference_sequence = self.reference_sequences[reference_id]
            set_valid = True

            for mutation in mutation_set.mutations:
                # Check if mutation position is within sequence bounds
                if isinstance(mutation, AminoAcidMutation):
                    if mutation.position >= len(reference_sequence):
                        validation_results["invalid_mutation_sets"].append(
                            {
                                "mutation_set": set_name,
                                "reference_id": reference_id,
                                "mutation": str(mutation),
                                "error": f"Position {mutation.position} exceeds sequence length (0-indexed)",
                            }
                        )
                        set_valid = False
                        continue
                
                elif isinstance(mutation, CodonMutation):
                    if mutation.nucleotide_end > len(reference_sequence):
                        validation_results["invalid_mutation_sets"].append(
                            {
                                "mutation_set": set_name,
                                "reference_id": reference_id,
                                "mutation": str(mutation),
                                "error": f"Position {mutation.position} exceeds sequence length (0-indexed)",
                            }
                        )
                        set_valid = False
                        continue  
                # Check if wild type matches reference for amino acid mutations
                if isinstance(mutation, AminoAcidMutation) and isinstance(
                    reference_sequence, ProteinSequence
                ):
                    try:
                        ref_residue = reference_sequence.get_residue(mutation.position)
                        if ref_residue != mutation.wild_amino_acid:
                            validation_results["position_mismatches"].append(
                                {
                                    "mutation_set": set_name,
                                    "reference_id": reference_id,
                                    "mutation": str(mutation),
                                    "expected": ref_residue,
                                    "found": mutation.wild_amino_acid,
                                    "position": mutation.position,
                                }
                            )
                    except IndexError:
                        validation_results["invalid_mutation_sets"].append(
                            {
                                "mutation_set": set_name,
                                "reference_id": reference_id,
                                "mutation": str(mutation),
                                "error": f"Position {mutation.position} out of range",
                            }
                        )
                        set_valid = False

                # Check codon mutations for nucleotide sequences
                elif isinstance(mutation, CodonMutation) and isinstance(
                    reference_sequence, (DNASequence, RNASequence)
                ):
                    try:
                        # Assuming position is codon position, get the codon at this position
                        start_pos = mutation.nucleotide_start
                        if start_pos + 3 <= len(reference_sequence):
                            ref_codon = str(
                                reference_sequence[start_pos : start_pos + 3]
                            ).upper()
                            if ref_codon != mutation.wild_codon:
                                validation_results["position_mismatches"].append(
                                    {
                                        "mutation_set": set_name,
                                        "reference_id": reference_id,
                                        "mutation": str(mutation),
                                        "expected": ref_codon,
                                        "found": mutation.wild_codon,
                                        "position": mutation.position,
                                    }
                                )
                        else:
                            validation_results["invalid_mutation_sets"].append(
                                {
                                    "mutation_set": set_name,
                                    "reference_id": reference_id,
                                    "mutation": str(mutation),
                                    "error": f"Codon position {mutation.position} exceeds sequence bounds",
                                }
                            )
                            set_valid = False
                    except Exception as e:
                        validation_results["invalid_mutation_sets"].append(
                            {
                                "mutation_set": set_name,
                                "reference_id": reference_id,
                                "mutation": str(mutation),
                                "error": f"Error validating codon: {str(e)}",
                            }
                        )
                        set_valid = False

            if set_valid:
                validation_results["valid_mutation_sets"].append(
                    {
                        "mutation_set": set_name,
                        "reference_id": reference_id,
                    }
                )

        return validation_results

    def to_dataframe(self) -> pd.DataFrame:
        """Convert dataset to pandas DataFrame"""
        if self._df is None:
            data = []

            for i, mutation_set in tqdm(
                enumerate(self.mutation_sets), desc="Converting dataset to DataFrame: "
            ):
                reference_id = self.mutation_set_references[
                    i
                ]  # This should always exist
                reference_sequence = self.reference_sequences[reference_id]
                label = self.mutation_set_labels.get(i)

                base_data = {
                    "mutation_set_id": i,
                    "mutation_set_name": mutation_set.name,
                    "reference_id": reference_id,
                    "reference_sequence_name": reference_sequence.name,
                    "reference_sequence_length": len(reference_sequence),
                    "reference_sequence_type": type(reference_sequence).__name__,
                    "num_mutations": len(mutation_set),
                    "is_single_mutation": mutation_set.is_single_mutation(),
                    "is_valid": mutation_set.validate_all(),
                    "positions": ",".join(map(str, mutation_set.get_positions())),
                    "mutation_subtype": mutation_set.mutation_subtype,
                    "label": label,
                }

                # Add mutation set metadata
                for key, value in mutation_set.metadata.items():
                    base_data[f"set_{key}"] = value

                # Add mutation information
                for j, mutation in enumerate(mutation_set.mutations):
                    mutation_data = base_data.copy()
                    mutation_data.update(
                        {
                            "mutation_id": j,
                            "mutation_type": mutation.type,
                            "mutation_string": str(mutation),
                            "position": mutation.position,
                            "mutation_category": mutation.get_mutation_category(),
                        }
                    )

                    # Add amino acid mutation-specific data
                    if isinstance(mutation, AminoAcidMutation):
                        mutation_data.update(
                            {
                                "wild_amino_acid": mutation.wild_amino_acid,
                                "mutant_amino_acid": mutation.mutant_amino_acid,
                                "effect_type": mutation.effect_type,
                                "is_synonymous": mutation.is_synonymous(),
                                "is_nonsense": mutation.is_nonsense(),
                                "is_missense": mutation.is_missense(),
                            }
                        )

                        # Add reference residue if available
                        if isinstance(reference_sequence, ProteinSequence):
                            try:
                                ref_residue = reference_sequence.get_residue(
                                    mutation.position
                                )
                                mutation_data["reference_residue"] = ref_residue
                                mutation_data["wild_type_matches_reference"] = (
                                    ref_residue == mutation.wild_amino_acid
                                )
                            except IndexError:
                                mutation_data["reference_residue"] = None
                                mutation_data["wild_type_matches_reference"] = False

                    # Add codon mutation-specific data
                    elif isinstance(mutation, CodonMutation):
                        mutation_data.update(
                            {
                                "wild_codon": mutation.wild_codon,
                                "mutant_codon": mutation.mutant_codon,
                                "seq_type": mutation.seq_type,
                            }
                        )

                        # Add reference codon if available
                        if isinstance(reference_sequence, (DNASequence, RNASequence)):
                            try:
                                # Get the codon at this position (assuming position is codon position)
                                start_pos = mutation.nucleotide_start
                                if start_pos + 3 <= len(reference_sequence):
                                    ref_codon = str(
                                        reference_sequence[start_pos : start_pos + 3]
                                    )
                                    mutation_data["reference_codon"] = ref_codon
                                    mutation_data["wild_codon_matches_reference"] = (
                                        ref_codon.upper() == mutation.wild_codon
                                    )
                                else:
                                    mutation_data["reference_codon"] = None
                                    mutation_data["wild_codon_matches_reference"] = (
                                        False
                                    )
                            except Exception:
                                mutation_data["reference_codon"] = None
                                mutation_data["wild_codon_matches_reference"] = False

                    # Add mutation metadata
                    for key, value in mutation.metadata.items():
                        mutation_data[f"mutation_{key}"] = value

                    data.append(mutation_data)

            self._df = pd.DataFrame(data)

        return self._df

    def filter_by_reference(self, reference_id: str) -> "MutationDataset":
        """Filter dataset to only include mutation sets from a specific reference sequence"""
        if reference_id not in self.reference_sequences:
            raise ValueError(f"Reference sequence with ID '{reference_id}' not found")

        filtered_sets = []
        filtered_references = []
        filtered_labels = []

        for i, mutation_set in tqdm(
            enumerate(self.mutation_sets), desc="Filtering by reference: "
        ):
            if self.mutation_set_references[i] == reference_id:
                filtered_sets.append(mutation_set)
                filtered_references.append(reference_id)
                filtered_labels.append(self.mutation_set_labels.get(i))

        filtered_dataset = MutationDataset(
            name=f"{self.name}_{reference_id}" if self.name else reference_id
        )
        filtered_dataset.add_reference_sequence(
            reference_id, self.reference_sequences[reference_id]
        )
        filtered_dataset.add_mutation_sets(
            filtered_sets, filtered_references, filtered_labels
        )

        return filtered_dataset

    def filter_by_mutation_type(
        self, mutation_type: Type[BaseMutation]
    ) -> "MutationDataset":
        """Filter dataset by mutation type"""
        filtered_sets = []
        filtered_references = []
        filtered_labels = []

        for i, mutation_set in tqdm(
            enumerate(self.mutation_sets), desc="Filtering by mutation type: "
        ):
            # Filter mutations by type
            filtered_mutations = [
                m for m in mutation_set.mutations if isinstance(m, mutation_type)
            ]

            if filtered_mutations:
                # Create new mutation set with filtered mutations
                if mutation_type == AminoAcidMutation:
                    new_set = AminoAcidMutationSet(
                        mutations=filtered_mutations,  # type: ignore
                        name=(
                            f"{mutation_set.name}_filtered"
                            if mutation_set.name
                            else "filtered"
                        ),
                        metadata=mutation_set.metadata.copy(),
                    )
                elif mutation_type == CodonMutation:
                    new_set = CodonMutationSet(
                        mutations=filtered_mutations,  # type: ignore
                        name=(
                            f"{mutation_set.name}_filtered"
                            if mutation_set.name
                            else "filtered"
                        ),
                        metadata=mutation_set.metadata.copy(),
                    )
                else:
                    new_set = MutationSet(
                        mutations=filtered_mutations,
                        mutation_type=mutation_type,
                        name=(
                            f"{mutation_set.name}_filtered"
                            if mutation_set.name
                            else "filtered"
                        ),
                        metadata=mutation_set.metadata.copy(),
                    )

                filtered_sets.append(new_set)
                ref_id = self.mutation_set_references[i]
                filtered_references.append(ref_id)
                filtered_labels.append(self.mutation_set_labels.get(i))

        filtered_dataset = MutationDataset(
            name=f"{self.name}_filtered" if self.name else "filtered"
        )
        # Copy all reference sequences that are still needed
        needed_refs = set(filtered_references)
        for ref_id in needed_refs:
            if ref_id is not None:
                filtered_dataset.add_reference_sequence(
                    ref_id, self.reference_sequences[ref_id]
                )

        filtered_dataset.add_mutation_sets(
            filtered_sets, filtered_references, filtered_labels
        )
        return filtered_dataset

    def filter_by_effect_type(self, effect_type: str) -> "MutationDataset":
        """Filter dataset by amino acid mutation effect type (synonymous, missense, nonsense)"""
        filtered_sets = []
        filtered_references = []
        filtered_labels = []

        for i, mutation_set in tqdm(
            enumerate(self.mutation_sets), desc="Filtering by effect type: "
        ):
            # Filter amino acid mutations by effect type
            filtered_mutations = []
            for mutation in mutation_set.mutations:
                if isinstance(mutation, AminoAcidMutation):
                    if mutation.effect_type == effect_type:
                        filtered_mutations.append(mutation)

            if filtered_mutations:
                new_set = AminoAcidMutationSet(
                    mutations=filtered_mutations,
                    name=(
                        f"{mutation_set.name}_{effect_type}"
                        if mutation_set.name
                        else effect_type
                    ),
                    metadata=mutation_set.metadata.copy(),
                )
                filtered_sets.append(new_set)
                ref_id = self.mutation_set_references[i]  # Always exists now
                filtered_references.append(ref_id)
                filtered_labels.append(self.mutation_set_labels.get(i))

        filtered_dataset = MutationDataset(
            name=f"{self.name}_{effect_type}" if self.name else effect_type
        )
        # Copy all reference sequences that are still needed
        needed_refs = set(filtered_references)
        for ref_id in needed_refs:
            if ref_id is not None:
                filtered_dataset.add_reference_sequence(
                    ref_id, self.reference_sequences[ref_id]
                )

        filtered_dataset.add_mutation_sets(
            filtered_sets, filtered_references, filtered_labels
        )
        return filtered_dataset

    def get_statistics(self) -> Dict[str, Any]:
        """Get basic statistics about the dataset"""
        total_sets = len(self.mutation_sets)
        total_mutations = sum(len(ms) for ms in self.mutation_sets)
        single_mutation_sets = sum(
            1 for ms in self.mutation_sets if ms.is_single_mutation()
        )
        multiple_mutation_sets = total_sets - single_mutation_sets

        mutation_types = {}
        mutation_categories = {}
        effect_types = {}
        reference_stats = {}

        # Statistics by reference sequence
        for ref_id, sequence in tqdm(
            self.reference_sequences.items(), desc="Statistics - ref seq: "
        ):
            reference_stats[ref_id] = {
                "sequence_name": sequence.name,
                "sequence_length": len(sequence),
                "sequence_type": type(sequence).__name__,
                "mutation_sets": 0,
                "mutations": 0,
            }

        for i, mutation_set in tqdm(
            enumerate(self.mutation_sets), desc="Statistics - mutation sets: "
        ):
            ref_id = self.mutation_set_references[i]  # Always exists now
            if ref_id in reference_stats:
                reference_stats[ref_id]["mutation_sets"] += 1
                reference_stats[ref_id]["mutations"] += len(mutation_set)

            for mutation in mutation_set.mutations:
                # Count mutation types
                mut_type = mutation.type
                mutation_types[mut_type] = mutation_types.get(mut_type, 0) + 1

                # Count mutation categories
                category = mutation.get_mutation_category()
                mutation_categories[category] = mutation_categories.get(category, 0) + 1

                # Count effect types for amino acid mutations
                if isinstance(mutation, AminoAcidMutation):
                    effect = mutation.effect_type
                    effect_types[effect] = effect_types.get(effect, 0) + 1

        stats = {
            "total_mutation_sets": total_sets,
            "total_mutations": total_mutations,
            "single_mutation_sets": single_mutation_sets,
            "multiple_mutation_sets": multiple_mutation_sets,
            "mutation_types": mutation_types,
            "mutation_categories": mutation_categories,
            "effect_types": effect_types,
            "average_mutations_per_set": (
                total_mutations / total_sets if total_sets > 0 else 0
            ),
            "reference_sequences": reference_stats,
            "num_reference_sequences": len(self.reference_sequences),
        }

        return stats

    def get_position_coverage(
        self, reference_id: Optional[str] = None
    ) -> Dict[str, Any]:
        """Get statistics about position coverage across reference sequences"""
        if reference_id is not None:
            if reference_id not in self.reference_sequences:
                raise ValueError(
                    f"Reference sequence with ID '{reference_id}' not found"
                )
            return self._get_single_sequence_coverage(reference_id)
        else:
            # Get coverage for all sequences
            coverage_stats = {}
            for ref_id in self.reference_sequences:
                coverage_stats[ref_id] = self._get_single_sequence_coverage(ref_id)
            return coverage_stats

    def _get_single_sequence_coverage(self, reference_id: str) -> Dict[str, Any]:
        """Get position coverage for a single reference sequence"""
        sequence = self.reference_sequences[reference_id]

        mutations = [
            mutation
            for i, mutation_set in enumerate(self.mutation_sets)
            if self.mutation_set_references[i] == reference_id
            for mutation in mutation_set.mutations
        ]

        all_positions = {mutation.position for mutation in mutations}

        if mutations and all(
            isinstance(mutation, AminoAcidMutation)
            for mutation in mutations
        ):
            if not isinstance(sequence, ProteinSequence):
                raise ValueError(
                    f"Amino-acid mutations require a ProteinSequence, "
                    f"but reference {reference_id!r} is "
                    f"{type(sequence).__name__}"
                )

            position_unit = "residue"
            position_space_length = len(sequence)

        elif mutations and all(
            isinstance(mutation, CodonMutation)
            for mutation in mutations
        ):
            if not isinstance(sequence, (DNASequence, RNASequence)):
                raise ValueError(
                    f"Codon mutations require a DNASequence or RNASequence, "
                    f"but reference {reference_id!r} is "
                    f"{type(sequence).__name__}"
                )

            if len(sequence) % 3 != 0:
                raise ValueError(
                    f"Reference sequence {reference_id!r} has length "
                    f"{len(sequence)}, which is not divisible by 3"
                )

            position_unit = "codon"
            position_space_length = len(sequence) // 3

        elif mutations:
            mutation_types = {
                type(mutation).__name__
                for mutation in mutations
            }
            raise ValueError(
                f"Cannot calculate position coverage for reference "
                f"{reference_id!r} with mixed mutation types: "
                f"{sorted(mutation_types)}"
            )

        else:
            # No mutation sets are associated with this reference.
            if isinstance(sequence, ProteinSequence):
                position_unit = "residue"
                position_space_length = len(sequence)
            elif isinstance(sequence, (DNASequence, RNASequence)):
                if len(sequence) % 3 != 0:
                    raise ValueError(
                        f"Reference sequence {reference_id!r} has length "
                        f"{len(sequence)}, which is not divisible by 3"
                    )

                position_unit = "codon"
                position_space_length = len(sequence) // 3
            else:
                raise TypeError(
                    f"Unsupported reference sequence type: "
                    f"{type(sequence).__name__}"
                )

        invalid_positions = sorted(
            position
            for position in all_positions
            if position < 0 or position >= position_space_length
        )

        if invalid_positions:
            raise ValueError(
                f"Reference {reference_id!r} contains out-of-range "
                f"{position_unit} positions: {invalid_positions[:10]}"
            )

        covered_positions = len(all_positions)
        uncovered_positions = position_space_length - covered_positions

        coverage_percentage = (
            covered_positions / position_space_length * 100
            if position_space_length > 0
            else 0.0
        )

        return {
            "reference_id": reference_id,
            "sequence_name": sequence.name,
            "sequence_length": len(sequence),
            "sequence_type": type(sequence).__name__,
            "position_unit": position_unit,
            "position_space_length": position_space_length,
            "covered_positions": covered_positions,
            "uncovered_positions": uncovered_positions,
            "coverage_percentage": coverage_percentage,
            "position_list": sorted(all_positions),
        }

    def convert_codon_to_amino_acid_sets(
        self, convert_labels: bool = False
    ) -> "MutationDataset":
        """
        Convert all codon mutation sets to amino acid mutation sets

        Parameters
        ----------
        convert_labels : bool, default: False
            Whether to save the labels with the mutation sets (default: False)
        """
        converted_sets = []
        converted_references = []
        converted_labels = []

        for i, mutation_set in enumerate(self.mutation_sets):
            if isinstance(mutation_set, CodonMutationSet):
                aa_set = mutation_set.to_amino_acid_mutation_set()
                converted_sets.append(aa_set)
            else:
                converted_sets.append(mutation_set)

            ref_id = self.mutation_set_references[i]
            converted_references.append(ref_id)
            converted_labels.append(self.get_mutation_set_label(i))

        converted_dataset = MutationDataset(
            name=f"{self.name}_aa_converted" if self.name else "aa_converted"
        )

        # Copy all reference sequences
        for ref_id, sequence in self.reference_sequences.items():
            converted_dataset.add_reference_sequence(ref_id, sequence)

        if not convert_labels:
            converted_dataset.add_mutation_sets(converted_sets, converted_references)
        else:
            converted_dataset.add_mutation_sets(
                converted_sets, converted_references, converted_labels
            )
        return converted_dataset

    @staticmethod
    def _sanitize_filename(name: str) -> str:
        """Sanitize a string to be used as a filename or directory name"""
        import re

        # Replace invalid characters with underscores
        sanitized = re.sub(r'[<>:"/\\|?*]', "_", name)
        # Remove leading/trailing whitespace and dots
        sanitized = sanitized.strip(". ")
        # Ensure it's not empty
        if not sanitized:
            sanitized = "unnamed"
        # Limit length to avoid filesystem issues
        if len(sanitized) > 200:
            sanitized = sanitized[:200]
        return sanitized

    def save_by_reference(self, base_dir: Union[str, Path]) -> None:
        """
        Save dataset by reference_id, creating separate folders for each reference.

        For each reference_id, creates:

        - {base_dir}/{reference_id}/data.csv: mutation data with columns [mutation_name, mutated_sequence, label]
        - {base_dir}/{reference_id}/wt.fasta: wild-type reference sequence
        - {base_dir}/{reference_id}/metadata.json: statistics and metadata for this reference

        Parameters
        ----------
        base_dir : Union[str, Path]
            Base directory to create reference folders in
        """
        base_path = Path(base_dir)
        base_path.mkdir(parents=True, exist_ok=True)

        tqdm.write(f"Saving dataset by reference to: {base_path}")

        # Group mutation sets and statistics by reference ID.
        ref_data = {}

        for i, mutation_set in tqdm(
            enumerate(self.mutation_sets),
            total=len(self.mutation_sets),
            desc="Grouping mutation sets",
        ):
            ref_id = self.mutation_set_references[i]

            if ref_id not in ref_data:
                ref_data[ref_id] = {
                    "mutation_sets": [],
                    "total_mutations": 0,
                    "covered_positions": set(),
                    "mutation_types": set(),
                    "unique_labels": set(),
                }

            data = ref_data[ref_id]

            data["mutation_sets"].append((i, mutation_set))
            data["total_mutations"] += len(mutation_set)
            data["covered_positions"].update(
                mutation_set.get_positions()
            )
            data["mutation_types"].update(
                type(mutation)
                for mutation in mutation_set.mutations
            )

            label = self.mutation_set_labels.get(i)

            if label is None or pd.isna(label):
                data["unique_labels"].add("unlabeled")
            else:
                data["unique_labels"].add(str(label))

        # Process each reference independently.
        for ref_id, data in tqdm(
            ref_data.items(),
            total=len(ref_data),
            desc="Saving data by reference",
        ):
            sanitized_ref_id = self._sanitize_filename(str(ref_id))
            ref_dir = base_path / sanitized_ref_id
            ref_dir.mkdir(parents=True, exist_ok=True)

            ref_sequence = self.reference_sequences[ref_id]
            mutation_types = data["mutation_types"]

            # Determine the coordinate system and number of valid positions.
            if mutation_types == {AminoAcidMutation}:
                if not isinstance(ref_sequence, ProteinSequence):
                    raise ValueError(
                        f"Amino-acid mutations for reference {ref_id!r} "
                        f"require a ProteinSequence, but found "
                        f"{type(ref_sequence).__name__}"
                    )

                position_unit = "residue"
                position_space_length = len(ref_sequence)

            elif mutation_types == {CodonMutation}:
                if not isinstance(
                    ref_sequence,
                    (DNASequence, RNASequence),
                ):
                    raise ValueError(
                        f"Codon mutations for reference {ref_id!r} "
                        f"require a DNASequence or RNASequence, but found "
                        f"{type(ref_sequence).__name__}"
                    )

                if len(ref_sequence) % 3 != 0:
                    raise ValueError(
                        f"Reference sequence {ref_id!r} has length "
                        f"{len(ref_sequence)}, which is not divisible by 3"
                    )

                position_unit = "codon"
                position_space_length = len(ref_sequence) // 3

            else:
                mutation_type_names = sorted(
                    mutation_type.__name__
                    for mutation_type in mutation_types
                )

                raise ValueError(
                    f"Reference {ref_id!r} contains mixed or unsupported "
                    f"mutation types: {mutation_type_names}"
                )

            # Validate that all positions are within the corresponding
            # residue or codon coordinate space.
            all_positions = data["covered_positions"]

            invalid_positions = sorted(
                position
                for position in all_positions
                if position < 0
                or position >= position_space_length
            )

            if invalid_positions:
                raise ValueError(
                    f"Reference {ref_id!r} contains out-of-range "
                    f"{position_unit} positions: "
                    f"{invalid_positions[:10]}"
                )

            covered_positions = len(all_positions)
            uncovered_positions = (
                position_space_length - covered_positions
            )

            coverage_percentage = (
                covered_positions
                / position_space_length
                * 100
                if position_space_length > 0
                else 0.0
            )

            successful_rows = []
            failed_rows = []

            # Apply each mutation set to the reference sequence.
            for set_index, mutation_set in data["mutation_sets"]:
                mutation_name = str(mutation_set)
                label = self.mutation_set_labels.get(
                    set_index,
                    "",
                )

                try:
                    mutated_sequence = ref_sequence.apply_mutation(
                        mutation_set
                    )

                    successful_rows.append(
                        {
                            "mutation_name": mutation_name,
                            "mutated_sequence": str(mutated_sequence),
                            "label": label,
                        }
                    )

                except Exception as error:
                    error_message = (
                        f"{type(error).__name__}: {error}"
                    )

                    failed_rows.append(
                        {
                            "mutation_name": mutation_name,
                            "label": label,
                            "error_message": error_message,
                        }
                    )

                    tqdm.write(
                        f"Warning: Could not apply mutation "
                        f"{mutation_name!r} for reference {ref_id!r}: "
                        f"{error_message}"
                    )

            # Save only successfully generated sequences to data.csv.
            data_columns = [
                "mutation_name",
                "mutated_sequence",
                "label",
            ]

            pd.DataFrame(
                successful_rows,
                columns=data_columns,
            ).to_csv(
                ref_dir / "data.csv",
                index=False,
            )

            # Save failed records separately instead of placing error text
            # inside the mutated_sequence column.
            failed_path = ref_dir / "failed_data.csv"

            if failed_rows:
                failed_columns = [
                    "mutation_name",
                    "label",
                    "error_message",
                ]

                pd.DataFrame(
                    failed_rows,
                    columns=failed_columns,
                ).to_csv(
                    failed_path,
                    index=False,
                )

            elif failed_path.exists():
                # Remove a stale failure file from an earlier run.
                failed_path.unlink()

            # Save the wild-type reference sequence.
            sequence_name = (
                ref_sequence.name
                if ref_sequence.name
                else str(ref_id)
            )

            with open(
                ref_dir / "wt.fasta",
                "w",
                encoding="utf-8",
            ) as fasta_file:
                fasta_file.write(
                    f">{sequence_name}\n{ref_sequence}\n"
                )

            # sequence_length is the actual number of sequence characters:
            # amino acids for proteins and nucleotides for DNA/RNA.
            metadata = {
                "reference_id": ref_id,
                "sequence_name": ref_sequence.name,
                "sequence_type": type(ref_sequence).__name__,
                "sequence_length": len(ref_sequence),
                "position_unit": position_unit,
                "position_space_length": position_space_length,
                "num_mutation_sets": len(data["mutation_sets"]),
                "num_saved_mutation_sets": len(successful_rows),
                "num_failed_mutation_sets": len(failed_rows),
                "total_mutations": data["total_mutations"],
                "covered_positions": covered_positions,
                "uncovered_positions": uncovered_positions,
                "coverage_percentage": coverage_percentage,
                "num_unique_labels": len(data["unique_labels"]),
                "has_unlabeled": (
                    "unlabeled" in data["unique_labels"]
                ),
                "dataset_name": self.name,
            }

            # Add sequence-type-specific length information.
            if isinstance(
                ref_sequence,
                (DNASequence, RNASequence),
            ):
                metadata.update(
                    {
                        "nucleotide_length": len(ref_sequence),
                        "codon_length": len(ref_sequence) // 3,
                    }
                )

            elif isinstance(ref_sequence, ProteinSequence):
                metadata["residue_length"] = len(ref_sequence)

            with open(
                ref_dir / "metadata.json",
                "w",
                encoding="utf-8",
            ) as metadata_file:
                json.dump(
                    metadata,
                    metadata_file,
                    indent=2,
                    default=str,
                )

    def save(
        self,
        filepath: str,
        save_type: Optional[Literal["mutcleaner", "pickle", "dataframe"]] = "mutcleaner",
    ):
        """
        Save the dataset to files.

        Parameters
        ----------
        filepath : str
            Base filepath (without extension)
        save_type : Optional[Literal["mutcleaner", "pickle", "dataframe"]], default="mutcleaner"
            Type of save format ("mutcleaner", "dataframe" or "pickle")

        For save_type="dataframe":
            - Saves mutations as {filepath}.csv
            - Saves reference sequences as {filepath}_refs.pkl
            - Saves metadata as {filepath}_meta.json

        For save_type="pickle":
            - Saves entire dataset as {filepath}.pkl

        Examples
        --------
        >>> dataset.save("my_study", "dataframe")
        >>> # Creates: my_study.csv, my_study_refs.pkl, my_study_meta.json
        """

        base_path = Path(filepath)

        if save_type == "dataframe":
            # Save mutations as CSV
            df = self.to_dataframe()
            csv_path = base_path.with_suffix(".csv")
            df.to_csv(csv_path, index=False)

            # Save reference sequences as pickle
            refs_path = base_path.with_suffix("").with_name(
                f"{base_path.name}_refs.pkl"
            )
            with open(refs_path, "wb") as f:
                pickle.dump(self.reference_sequences, f)

            # Save dataset metadata as JSON
            meta_path = base_path.with_suffix("").with_name(
                f"{base_path.name}_meta.json"
            )
            dataset_meta = {
                "name": self.name,
                "metadata": self.metadata,
                "save_type": save_type,
                "num_mutation_sets": len(self.mutation_sets),
                "num_reference_sequences": len(self.reference_sequences),
            }
            with open(meta_path, "w") as f:
                json.dump(dataset_meta, f, indent=2)

            tqdm.write(f"Dataset saved to:")
            tqdm.write(f"  Mutations: {csv_path}")
            tqdm.write(f"  References: {refs_path}")
            tqdm.write(f"  Metadata: {meta_path}")

        elif save_type == "pickle":
            # Save entire dataset as pickle
            pkl_path = base_path.with_suffix(".pkl")
            with open(pkl_path, "wb") as f:
                pickle.dump(self, f)
            tqdm.write(f"Dataset saved to: {pkl_path}")

        elif save_type == "mutcleaner":
            # Save as Mutcleaner format
            if base_path.suffix != "":
                raise ValueError(
                    f"Invalid Mutcleaner save format. Expected folder but got {base_path.suffix}."
                )
            self.save_by_reference(base_path)

        else:
            raise ValueError(
                f"Unsupported save_type: {save_type}. Use 'mutcleaner', 'dataframe' or 'pickle'"
            )

    # ====== load ======
    @classmethod
    def load_by_reference(
        cls,
        base_dir: Union[str, Path],
        dataset_name: Optional[str] = None,
        is_zero_based: bool = True,
    ) -> "MutationDataset":
        """
        Load a dataset from mutcleaner reference-based format.

        Expected directory structure::

            base_dir/
            ├── reference_id_1/
            │   ├── data.csv
            │   ├── wt.fasta
            │   └── metadata.json
            ├── reference_id_2/
            │   ├── data.csv
            │   ├── wt.fasta
            │   └── metadata.json
            └── ...

        Parameters
        ----------
        base_dir : Union[str, Path]
            Base directory containing reference folders
        dataset_name : Optional[str], default=None
            Optional name for the loaded dataset
        is_zero_based : bool, default=True
            Whether origin mutation positions are zero-based

        Returns
        -------
        MutationDataset instance
        """

        base_path = Path(base_dir)
        if not base_path.exists():
            raise FileNotFoundError(f"Base directory not found: {base_path}")
        if not base_path.is_dir():
            raise ValueError(f"Path is not a directory: {base_path}")

        tqdm.write(f"Loading dataset from: {base_path}")

        # Find all reference directories
        ref_dirs = [d for d in base_path.iterdir() if d.is_dir()]
        if not ref_dirs:
            raise ValueError(f"No reference directories found in {base_path}")

        # Create new dataset
        dataset = cls(name=dataset_name)

        total_loaded = 0
        skipped_dirs = []

        # Process each reference directory
        for ref_dir in tqdm(ref_dirs, desc="Loading dataset from mutcleaner format "):
            # Required files
            data_path = ref_dir / "data.csv"
            fasta_path = ref_dir / "wt.fasta"
            metadata_path = ref_dir / "metadata.json"

            # Check required files exist
            if not data_path.exists() or not fasta_path.exists():
                missing_files = []
                if not data_path.exists():
                    missing_files.append("data.csv")
                if not fasta_path.exists():
                    missing_files.append("wt.fasta")
                skipped_dirs.append(
                    f"{ref_dir.name} (missing: {', '.join(missing_files)})"
                )
                continue

            # Load metadata to get original reference_id and sequence_type
            original_ref_id = ref_dir.name  # Default to directory name
            sequence_type = ProteinSequence  # Default sequence type

            if metadata_path.exists():
                try:
                    with open(metadata_path, "r") as f:
                        metadata = json.load(f)
                    original_ref_id = metadata.get("reference_id", ref_dir.name)
                    sequence_type_name = metadata.get(
                        "sequence_type", "ProteinSequence"
                    )
                    sequence_type = SEQUENCE_TYPE_MAP.get(
                        sequence_type_name, ProteinSequence
                    )
                except Exception as e:
                    print(f"Warning: Could not load metadata for {ref_dir.name}: {e}")

            try:
                # Load reference sequence from FASTA
                sequences = load_sequences_from_fasta(
                    fasta_path, sequence_type, header_func=lambda x: (x, "")
                )
                if not sequences:
                    skipped_dirs.append(f"{ref_dir.name} (empty FASTA)")
                    continue

                ref_sequence = list(sequences.values())[0]  # Get first sequence
                dataset.add_reference_sequence(original_ref_id, ref_sequence)

                # Load mutation data
                df_ref = pd.read_csv(data_path)
                required_cols = ["mutation_name", "mutated_sequence", "label"]
                missing_cols = [
                    col for col in required_cols if col not in df_ref.columns
                ]

                if missing_cols:
                    skipped_dirs.append(
                        f"{ref_dir.name} (missing columns: {', '.join(missing_cols)})"
                    )
                    continue

                # Batch process mutations for this reference
                mutation_sets_added = 0
                failed_mutations = 0

                for _, row in df_ref.iterrows():
                    mutation_name = row["mutation_name"]
                    label = row["label"]

                    # Parse mutation from mutation_name
                    try:
                        mutation_set = MutationSet.from_string(
                            mutation_name, sep=",", is_zero_based=is_zero_based
                        )
                        dataset.add_mutation_set(mutation_set, original_ref_id, label)
                        mutation_sets_added += 1
                    except Exception:
                        failed_mutations += 1
                        continue

                tqdm.write(
                    f"  Loaded {original_ref_id}: {mutation_sets_added} mutation sets",
                    end="",
                )
                if failed_mutations > 0:
                    tqdm.write(f" ({failed_mutations} failed)")
                else:
                    tqdm.write("")

                total_loaded += mutation_sets_added

            except Exception as e:
                skipped_dirs.append(f"{ref_dir.name} (error: {str(e)})")
                continue

        # Report results
        if skipped_dirs:
            tqdm.write(f"Skipped {len(skipped_dirs)} directories:")
            for skip_info in skipped_dirs:
                tqdm.write(f"  - {skip_info}")

        if len(dataset) == 0:
            raise ValueError("No valid mutation sets were loaded")

        tqdm.write(
            f"Successfully loaded dataset with {len(dataset)} mutation sets from {len(dataset.reference_sequences)} references"
        )
        return dataset

    @classmethod
    def load(cls, filepath: str, load_type: Optional[str] = None) -> "MutationDataset":
        """
        Load a dataset from files.

        Parameters
        ----------
        filepath : str
            Base filepath (with or without extension)
        load_type : Optional[str], default=None
            Type of load format ("mutcleaner", "dataframe" or "pickle").
            If None, auto-detect from file extension.

        Returns
        -------
            MutationDataset instance

        Examples
        --------
        >>> # Auto-detect from extension
        >>> dataset = MutationDataset.load("my_study.csv")
        >>> dataset = MutationDataset.load("my_study.pkl")

        >>> # Explicit type
        >>> dataset = MutationDataset.load("my_study", "dataframe")
        """
        base_path = Path(filepath)

        # Auto-detect load type from extension if not specified
        if load_type is None:
            if base_path.suffix == ".csv":
                load_type = "dataframe"
            elif base_path.suffix == ".pkl":
                load_type = "pickle"
            elif base_path.suffix == "":
                load_type = "mutcleaner"
            else:
                # Try dataframe format first
                load_type = "dataframe"
                base_path = base_path.with_suffix("")  # Remove any extension

        if load_type == "dataframe":
            # Remove extension to get base path
            if base_path.suffix == ".csv":
                base_path = base_path.with_suffix("")

            csv_path = base_path.with_suffix(".csv")
            refs_path = base_path.with_suffix("").with_name(
                f"{base_path.name}_refs.pkl"
            )
            meta_path = base_path.with_suffix("").with_name(
                f"{base_path.name}_meta.json"
            )

            # Check if files exist
            if not csv_path.exists():
                raise FileNotFoundError(f"CSV file not found: {csv_path}")
            if not refs_path.exists():
                raise FileNotFoundError(f"References file not found: {refs_path}")

            # Load mutations DataFrame
            df = pd.read_csv(csv_path)

            # Load reference sequences
            with open(refs_path, "rb") as f:
                reference_sequences = pickle.load(f)

            # Load metadata if available
            dataset_name = None
            dataset_metadata = {}
            if meta_path.exists():
                with open(meta_path, "r") as f:
                    meta = json.load(f)
                    dataset_name = meta.get("name")
                    dataset_metadata = meta.get("metadata", {})

            # Create dataset using `from_dataframe`
            dataset = cls.from_dataframe(df, reference_sequences, dataset_name)
            dataset.metadata = dataset_metadata

            tqdm.write(f"Dataset loaded from:")
            tqdm.write(f"  Mutations: {csv_path}")
            tqdm.write(f"  References: {refs_path}")
            if meta_path.exists():
                tqdm.write(f"  Metadata: {meta_path}")

            return dataset

        elif load_type == "pickle":
            pkl_path = base_path.with_suffix(".pkl")
            if not pkl_path.exists():
                raise FileNotFoundError(f"Pickle file not found: {pkl_path}")

            with open(pkl_path, "rb") as f:
                dataset = pickle.load(f)

            tqdm.write(f"Dataset loaded from: {pkl_path}")
            return dataset

        elif load_type == "mutcleaner":
            return cls.load_by_reference(base_path)

        else:
            raise ValueError(
                f"Unsupported load_type: {load_type}. Use 'dataframe', 'pickle' or 'mutcleaner'"
            )

    @classmethod
    def from_dataframe(
        cls,
        df: pd.DataFrame,
        reference_sequences: Dict[str, BaseSequence],
        name: Optional[str] = None,
        specific_mutation_type: Optional[Type[BaseMutation]] = None,
    ) -> "MutationDataset":
        """
        Create a MutationDataset from a DataFrame containing mutation data.

        This method reconstructs a MutationDataset from a flattened DataFrame representation,
        typically used for loading saved mutation datasets from files. The DataFrame should
        contain mutation information with each row representing a single mutation within
        mutation sets.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame containing mutation data with the following required columns:
            - 'mutation_set_id': Identifier for grouping mutations into sets
            - 'reference_id': Identifier for the reference sequence
            - 'mutation_string': String representation of the mutation
            - 'position': Position of the mutation in the sequence
            - 'mutation_type': Type of mutation ('amino_acid', 'codon_dna', 'codon_rna')

            Optional columns include:
            - 'mutation_set_name': Name of the mutation set
            - 'label': Label associated with the mutation set
            - 'wild_amino_acid': Wild-type amino acid (for amino acid mutations)
            - 'mutant_amino_acid': Mutant amino acid (for amino acid mutations)
            - 'wild_codon': Wild-type codon (for codon mutations)
            - 'mutant_codon': Mutant codon (for codon mutations)
            - 'set_*': Columns with 'set_' prefix for mutation set metadata
            - 'mutation_*': Columns with 'mutation_' prefix for individual mutation metadata

        reference_sequences : Dict[str, BaseSequence]
            Dictionary mapping reference sequence IDs to their corresponding BaseSequence
            objects. Must contain all reference sequences referenced in the DataFrame.

        name : Optional[str], default=None
            Optional name for the created MutationDataset.

        specific_mutation_type : Optional[BaseMutation], default=None
            The type of mutations to create. If None, will infer from first mutation
            must be provided when the mutation type is neither 'amino_acid' nor any 'codon_*' type.

        Returns
        -------
        MutationDataset
            A new MutationDataset instance populated with the mutation sets and
            reference sequences from the DataFrame.

        Raises
        ------
        ValueError
            If the DataFrame is empty, missing required columns, or references
            sequences not provided in reference_sequences dict.

        Notes
        -----
        - Mutations are grouped by 'mutation_set_id' to reconstruct mutation sets
        - The method automatically determines the appropriate mutation set type
        (AminoAcidMutationSet, CodonMutationSet, or generic MutationSet) based
        on the mutation types within each set
        - Metadata is extracted from columns with 'set_' and 'mutation_' prefixes
        - Only reference sequences that are actually used in the DataFrame are
        added to the dataset

        Examples
        --------
        >>> import pandas as pd
        >>> from sequences import ProteinSequence
        >>>
        >>> # Create sample DataFrame
        >>> df = pd.DataFrame({
        ...     'mutation_set_id': ['set1', 'set1', 'set2'],
        ...     'reference_id': ['prot1', 'prot1', 'prot2'],
        ...     'mutation_string': ['A1V', 'L2P', 'G5R'],
        ...     'position': [1, 2, 5],
        ...     'mutation_type': ['amino_acid', 'amino_acid', 'amino_acid'],
        ...     'wild_amino_acid': ['A', 'L', 'G'],
        ...     'mutant_amino_acid': ['V', 'P', 'R'],
        ...     'mutation_set_name': ['variant1', 'variant1', 'variant2'],
        ...     'label': ['pathogenic', 'pathogenic', 'benign']
        ... })
        >>>
        >>> # Define reference sequences
        >>> ref_seqs = {
        ...     'prot1': ProteinSequence('ALDEFG', name='protein1'),
        ...     'prot2': ProteinSequence('MKGLRK', name='protein2')
        ... }
        >>>
        >>> # Create MutationDataset
        >>> dataset = MutationDataset.from_dataframe(df, ref_seqs, name="my_dataset")
        >>> print(len(dataset.mutation_sets))
        2
        """
        if df.empty:
            raise ValueError("DataFrame cannot be empty")

        # Validate required columns
        required_cols = [
            "mutation_set_id",
            "reference_id",
            "mutation_string",
            "position",
            "mutation_type",
        ]
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"DataFrame missing required columns: {missing_cols}")

        # Validate that all referenced sequences are provided
        df_ref_ids = set(df["reference_id"].dropna().unique())
        provided_ref_ids = set(reference_sequences.keys())
        missing_refs = df_ref_ids - provided_ref_ids
        if missing_refs:
            raise ValueError(f"Missing reference sequences for IDs: {missing_refs}")

        # Create new dataset
        dataset = cls(name=name)

        # Add reference sequences
        for ref_id, sequence in tqdm(
            reference_sequences.items(), desc="Adding reference sequences"
        ):
            if ref_id in df_ref_ids:  # Only add sequences that are actually used
                dataset.add_reference_sequence(ref_id, sequence)

        # Rocognize metadata columns
        set_metadata_cols = [
            col
            for col in df.columns
            if col.startswith("set_") and col != "set_metadata"
        ]
        mutation_metadata_cols = [
            col
            for col in df.columns
            if col.startswith("mutation_")
            and col
            not in {
                "mutation_id",
                "mutation_type",
                "mutation_string",
                "mutation_category",
            }
        ]

        # Group by mutation set to rebuild mutation sets
        grouped = df.groupby("mutation_set_id", sort=False)

        for _, group in tqdm(grouped, desc="Reconstructing mutation sets"):
            # Get mutation set info
            set_info = group.iloc[0]
            set_name = set_info.get("mutation_set_name")
            reference_id = set_info["reference_id"]
            label = set_info.get("label")

            # Extract set metadata
            set_metadata = {
                col[4:]: value
                for col in set_metadata_cols
                if pd.notna(value := set_info[col])
            }

            # Create mutations from group
            mutations = []
            columns = list(group.columns)
            for values in group.values:
                row_dict = dict(zip(columns, values))
                mutation = cls._create_mutation_from_dict(
                    row_dict, mutation_metadata_cols, specific_mutation_type
                )
                mutations.append(mutation)

            # Create appropriate mutation set type
            if mutations:
                mutation_type = type(mutations[0])

                if mutation_type == AminoAcidMutation:
                    mutation_set = AminoAcidMutationSet(
                        mutations=mutations, name=set_name, metadata=set_metadata  # type: ignore
                    )
                elif mutation_type == CodonMutation:
                    mutation_set = CodonMutationSet(
                        mutations=mutations, name=set_name, metadata=set_metadata  # type: ignore
                    )
                else:
                    mutation_set = MutationSet(
                        mutations=mutations,
                        mutation_type=mutation_type,
                        name=set_name,
                        metadata=set_metadata,
                    )

                # Add to dataset
                dataset.add_mutation_set(mutation_set, reference_id, label)

        return dataset

    @staticmethod
    def _create_mutation_from_dict(
        row_dict: Dict[str, Any],
        metadata_cols: List[str],
        specific_mutation_type: Optional[Type[BaseMutation]] = None,
    ) -> BaseMutation:
        """
        Create a mutation object from a DataFrame row dict
        Used by the from_dataframe() method

        Parameters
        ----------
        row_dict: Dict[str, Any]
            A pandas row dict containing mutation data with the following required fields:
            - 'mutation_type': Type of mutation ('amino_acid', 'codon_*', or other)

            For amino acid mutations, also requires:
            - 'wild_amino_acid': Wild-type amino acid
            - 'mutant_amino_acid': Mutant amino acid

            For codon mutations, also requires:
            - 'wild_codon': Wild-type codon
            - 'mutant_codon': Mutant codon

            For other mutations:
            - 'mutation_string': String representation of the mutation

            Optional fields:
            - Any column starting with 'mutation_' (excluding specific reserved names)
            will be added as metadata to the mutation object

        mutation_metadata_cols : List[str]
            Metadata columns to be extracted and added to the mutation object

        specific_mutation_type : Optional[Type[BaseMutation]], default=None
            Only used when the mutation type is neither 'amino_acid' nor any 'codon_*' type.

        Returns
        -------
        BaseMutation
            An instance of the appropriate mutation subclass:
            - AminoAcidMutation for 'amino_acid' type
            - CodonMutation for types starting with 'codon_'
            - Inferred mutation type for other types (parsed from mutation_string)
        """
        mutation_type = row_dict["mutation_type"]
        position = int(row_dict["position"])

        # Filter metadata
        mutation_metadata = {
            key[9:]: row_dict[key]  # Remove "mutation_" prefix
            for key in metadata_cols
            if not pd.isna(row_dict[key])
        }

        if mutation_type == "amino_acid":
            return AminoAcidMutation(
                wild_type=row_dict["wild_amino_acid"],
                position=position,
                mutant_type=row_dict["mutant_amino_acid"],
                metadata=mutation_metadata,
            )

        elif mutation_type.startswith("codon_"):
            return CodonMutation(
                wild_type=row_dict["wild_codon"],
                position=position,
                mutant_type=row_dict["mutant_codon"],
                metadata=mutation_metadata,
            )

        else:
            # FIXME: need to handle other mutation types
            # Try to parse from mutation string as fallback
            if specific_mutation_type is None:
                raise ValueError(
                    f"Unsupported mutation type: {mutation_type}, "
                    f"you must provide a specific mutation type"
                )
            mutation_string = row_dict["mutation_string"]
            try:
                return MutationSet._create_mutation(
                    mutation_string,
                    mutation_type=specific_mutation_type,
                    is_zero_based=True,
                )
            except Exception as e:
                raise ValueError(f"Cannot create mutation from row: {e}")
