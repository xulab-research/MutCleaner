import mutcleaner
import mutcleaner.utils

from mutcleaner.utils import raw_data_downloader


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


class TestSubmoduleWildcardImport:
    """Test submodule wildcard imports against their __all__ lists"""

    def test_raw_data_downloader_wildcard_import(self):
        """Test `from mutcleaner.utils.raw_data_downloader import *` imports exactly __all__"""
        namespace = {}
        exec("from mutcleaner.utils.raw_data_downloader import *", namespace)
        imported_names = {name for name in namespace if name != "__builtins__"}
        assert imported_names == set(raw_data_downloader.__all__)

    def test_raw_data_downloader_all_names_resolve(self):
        """Test every downloader __all__ entry resolves on the module"""
        missing_names = [name for name in raw_data_downloader.__all__ if not hasattr(raw_data_downloader, name)]
        assert missing_names == []

    def test_utils_wildcard_import(self):
        """Test `from mutcleaner.utils import *` imports exactly its __all__"""
        namespace = {}
        exec("from mutcleaner.utils import *", namespace)
        imported_names = {name for name in namespace if name != "__builtins__"}
        assert imported_names == set(mutcleaner.utils.__all__)
