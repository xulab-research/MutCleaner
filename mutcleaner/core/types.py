# mutcleaner/core/types.py
"""Type variables used across mutcleaner.

Attributes
----------
SequenceType : typing.TypeVar
    Type variable bound to :class:`~mutcleaner.core.sequence.BaseSequence`.
MutationType : typing.TypeVar
    Bound to :class:`~mutcleaner.core.mutation.BaseMutation`.
MutationSetType : typing.TypeVar
    Bound to :class:`~mutcleaner.core.sequence.MutationSet`.
CleanerConfigType : typing.TypeVar
    Bound to :class:`~mutcleaner.cleaners.base_config.BaseCleanerConfig`.
"""

from __future__ import annotations

from typing import TypeVar, TYPE_CHECKING

if TYPE_CHECKING:
    from .sequence import BaseSequence, MutationSet
    from .mutation import BaseMutation
    from ..cleaners.base_config import BaseCleanerConfig

SequenceType = TypeVar("SequenceType", bound="BaseSequence")
MutationType = TypeVar("MutationType", bound="BaseMutation")
MutationSetType = TypeVar("MutationSetType", bound="MutationSet")
CleanerConfigType = TypeVar("CleanerConfigType", bound="BaseCleanerConfig")
