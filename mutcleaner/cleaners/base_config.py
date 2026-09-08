# mutcleaner/cleaners/base_config.py
from __future__ import annotations

import json
from abc import ABC, abstractmethod
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any, Dict, List, Type, Union

    from ..core.types import CleanerConfigType

__all__ = ["BaseCleanerConfig"]


def __dir__() -> List[str]:
    return __all__


@dataclass
class BaseCleanerConfig(ABC):
    """Base configuration class for all dataset cleaners

    This abstract base class provides common configuration functionality
    that can be inherited by specific cleaner configurations.

    Attributes
    ----------
    pipeline_name : str
        Name of the cleaning pipeline
    num_workers : int
        Default number of worker processes
    validate_config : bool
        Whether to validate configuration before use
    """

    # Common configuration options
    pipeline_name: str = field(default="base_cleaner")
    num_workers: int = field(default=16)
    validate_config: bool = field(default=True)

    def __post_init__(self):
        """Post-initialization validation"""
        if self.validate_config:
            self.validate()

    @abstractmethod
    def validate(self) -> None:
        """Validate the configuration

        This method should be implemented by subclasses to perform
        specific validation logic.

        Raises
        ------
        ValueError
            If configuration is invalid
        """
        # Common validations
        if self.num_workers < 1:
            raise ValueError(f"num_workers must be at least 1, got {self.num_workers}")

    @classmethod
    def from_dict(cls: Type[CleanerConfigType], config_dict: Dict[str, Any]) -> CleanerConfigType:
        """Create configuration object from dictionary

        Parameters
        ----------
        config_dict : Dict[str, Any]
            Dictionary containing configuration parameters

        Returns
        -------
        BaseCleanerConfig
            Configuration object
        """
        return cls(**config_dict)

    @classmethod
    def from_json(cls: Type[CleanerConfigType], json_path: Union[str, Path]) -> CleanerConfigType:
        """Load configuration from JSON file

        Parameters
        ----------
        json_path : Union[str, Path]
            Path to JSON configuration file

        Returns
        -------
        BaseCleanerConfig
            Configuration object

        Raises
        ------
        FileNotFoundError
            If configuration file does not exist
        """
        json_path = Path(json_path)
        if not json_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {json_path}")

        with open(json_path, "r") as f:
            config_dict = json.load(f)

        return cls.from_dict(config_dict)

    def to_dict(self, exclude_callables: bool = True) -> Dict[str, Any]:
        """Convert configuration to dictionary

        Parameters
        ----------
        exclude_callables : bool, optional
            Whether to exclude callable objects (functions, lambdas), by default True

        Returns
        -------
        Dict[str, Any]
            Dictionary representation of the configuration
        """
        data = asdict(self)

        if exclude_callables:
            # Remove any callable values that can't be serialized
            data = {k: v for k, v in data.items() if not callable(v) and not (isinstance(v, dict) and any(callable(vv) for vv in v.values()))}

        return data

    def to_json(self, json_path: Union[str, Path], **json_kwargs) -> None:
        """Save configuration to JSON file

        Parameters
        ----------
        json_path : Union[str, Path]
            Path where to save the JSON file
        **json_kwargs
            Additional arguments passed to json.dump
        """
        json_path = Path(json_path)
        json_path.parent.mkdir(parents=True, exist_ok=True)

        config_dict = self.to_dict(exclude_callables=True)

        with open(json_path, "w") as f:
            json.dump(config_dict, f, indent=2, **json_kwargs)

    def merge(self: CleanerConfigType, partial_config: Dict[str, Any]) -> CleanerConfigType:
        """Merge partial configuration with current configuration

        Parameters
        ----------
        partial_config : Dict[str, Any]
            Dictionary containing configuration values to update

        Returns
        -------
        BaseCleanerConfig
            New configuration object with merged values
        """
        current_dict = asdict(self)

        # Deep merge for nested dictionaries
        def deep_merge(base: dict, update: dict) -> dict:
            result = base.copy()
            for key, value in update.items():
                if key in result and isinstance(result[key], dict) and isinstance(value, dict):
                    result[key] = deep_merge(result[key], value)
                else:
                    result[key] = value
            return result

        merged_dict = deep_merge(current_dict, partial_config)
        return self.__class__.from_dict(merged_dict)

    def get_summary(self) -> str:
        """Get a human-readable summary of the configuration

        Returns
        -------
        str
            String summary of the configuration
        """
        lines = [f"{self.__class__.__name__} Configuration:"]
        for key, value in self.to_dict().items():
            if isinstance(value, dict):
                lines.append(f"  {key}:")
                for k, v in value.items():
                    lines.append(f"    {k}: {v}")
            else:
                lines.append(f"  {key}: {value}")
        return "\n".join(lines)
