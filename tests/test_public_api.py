import mutcleaner


class TestPackageWildcardImport:
    """Test package-level wildcard import and __all__ consistency"""

    def test_wildcard_import_from_package(self):
        """Test `from mutcleaner import *` imports exactly the names in __all__"""
        namespace = {}
        exec("from mutcleaner import *", namespace)
        imported_names = {name for name in namespace if name != "__builtins__"}
        assert imported_names == set(mutcleaner.__all__)

    def test_all_names_resolve_on_package(self):
        """Test every __all__ entry resolves as a mutcleaner attribute"""
        missing_names = [name for name in mutcleaner.__all__ if not hasattr(mutcleaner, name)]
        assert missing_names == []

    def test_all_has_no_duplicate_entries(self):
        """Test __all__ contains no duplicate names"""
        duplicate_names = [name for name in mutcleaner.__all__ if mutcleaner.__all__.count(name) > 1]
        assert duplicate_names == []

    def test_explicit_imports_unchanged(self):
        """Test existing explicit imports from mutcleaner still work"""
        from mutcleaner import (
            MutationDataset,
            Pipeline,
            codon_human_myoglobin_cleaner,
            protein_human_myoglobin_cleaner,
        )

        assert MutationDataset is mutcleaner.MutationDataset
        assert Pipeline is mutcleaner.Pipeline
        assert protein_human_myoglobin_cleaner is mutcleaner.protein_human_myoglobin_cleaner
        assert codon_human_myoglobin_cleaner is mutcleaner.codon_human_myoglobin_cleaner
