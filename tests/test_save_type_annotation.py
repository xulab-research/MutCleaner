import typing

from mutcleaner.core import dataset as dataset_module
from mutcleaner.core.dataset import MutationDataset


class TestSaveTypeAnnotation:
    """Test the save_type annotation of MutationDataset.save."""

    def test_save_type_annotation_allows_pickle(self):
        """The save_type annotation must list the pickle format."""
        type_globals = {**vars(dataset_module), "Literal": typing.Literal, "Optional": typing.Optional}
        save_type_hint = typing.get_type_hints(MutationDataset.save, globalns=type_globals)["save_type"]
        allowed_values = typing.get_args(typing.get_args(save_type_hint)[0])

        assert "pickle" in allowed_values
        assert " " not in allowed_values
