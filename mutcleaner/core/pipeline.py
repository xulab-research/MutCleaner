# mutcleaner/core/pipeline.py
from __future__ import annotations

from typing import (
    Any,
    Callable,
    NamedTuple,
    ParamSpec,
    Tuple,
    TypeVar,
    overload,
    TYPE_CHECKING,
)
from copy import deepcopy
import logging
from functools import wraps
import time
import pickle
from dataclasses import dataclass

if TYPE_CHECKING:
    from typing import (
        List,
        Dict,
        Optional,
        Union,
    )

__all__ = ["Pipeline", "create_pipeline", "multiout_step", "pipeline_step"]


def __dir__() -> List[str]:
    return __all__


class MultiOutput(NamedTuple):
    """Container for functions that return multiple outputs"""

    main: Any
    side: Dict[str, Any] = {}


@dataclass
class PipelineOutput:
    """Structured output from pipeline steps"""

    data: Any
    artifacts: Dict[str, Any]

    def __getitem__(self, key: str):
        """Allow dictionary-style access to artifacts"""
        return self.artifacts.get(key)


@dataclass
class DelayedStep:
    """Represents a delayed step that hasn't been executed yet"""

    name: str
    function: Callable
    args: tuple
    kwargs: dict

    def to_pipeline_step(self) -> "PipelineStep":
        """Convert to a PipelineStep for execution"""
        return PipelineStep(self.name, self.function, *self.args, **self.kwargs)


class PipelineStep:
    """Represents a single step in the pipeline"""

    def __init__(self, name: str, function: Callable, *args, **kwargs):
        self.name = name
        self.function = function
        self.args = args
        self.kwargs = kwargs
        self.result = None
        self.side_outputs = {}  # Store side outputs
        self.execution_time = None
        self.success = False
        self.error = None

    def execute(self, data: Any) -> Any:
        """Execute this pipeline step"""
        start_time = time.time()

        try:
            # Apply function with args and kwargs - always pass actual data
            result = self.function(data, *self.args, **self.kwargs)

            # Handle different output types based on step type
            step_type = getattr(self.function, "_step_type", "unknown")

            if isinstance(result, MultiOutput):
                # Result from @multiout_step decorated function
                self.result = result.main
                self.side_outputs = result.side
                final_result = result.main

            elif step_type == "multi_output":
                # This shouldn't happen if @multiout_step is working correctly
                raise RuntimeError(f"Function {self.function.__name__} is marked as multi_output " f"but didn't return MultiOutput. This indicates a decorator bug.")

            elif step_type == "single_output":
                # Function decorated with @pipeline_step - treat any result as single value
                self.result = result
                final_result = result

            else:
                # Undecorated function - treat as single output with warning
                import warnings

                warnings.warn(
                    f"Function {getattr(self.function, '__name__', 'unknown')} is not " f"decorated with @pipeline_step or @multiout_step. " f"Consider adding @pipeline_step for better pipeline integration.",
                    UserWarning,
                )
                self.result = result
                final_result = result

            self.success = True
            return final_result

        except Exception as e:
            self.success = False
            self.error = e
            raise
        finally:
            self.execution_time = time.time() - start_time

    def get_step_info(self) -> Dict[str, Any]:
        """Get detailed information about this step"""
        return {
            "name": self.name,
            "step_type": getattr(self.function, "_step_type", "unknown"),
            "expected_outputs": getattr(self.function, "_expected_output_count", 1),
            "output_names": getattr(self.function, "_output_names", ["main"]),
            "success": self.success,
            "execution_time": self.execution_time,
            "has_side_outputs": bool(self.side_outputs),
            "side_output_keys": (list(self.side_outputs.keys()) if self.side_outputs else []),
            "error": str(self.error) if self.error else None,
        }


class Pipeline:
    """Pipeline for processing data with pandas-style method chaining"""

    def __init__(self, data: Any = None, name: Optional[str] = None, logging_level: str = "INFO"):
        self.name = name or "Pipeline"
        self._data = data
        self._artifacts: Dict[str, Any] = {}
        self.steps: List[PipelineStep] = []
        self.delayed_steps: List[DelayedStep] = []
        self.results: List[Any] = []

        # Setup logging
        self.logger = logging.getLogger(self.name)
        self.logger.setLevel(getattr(logging, logging_level))

        # Add handler if logger doesn't have one
        if not self.logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)

    @property
    def data(self) -> Any:
        """
        Always return the actual data, never PipelineOutput.

        This ensures consistent user experience - pipeline.data can always be
        used with methods like .copy(), .append(), etc.
        """
        if isinstance(self._data, PipelineOutput):
            return self._data.data
        return self._data

    @property
    def artifacts(self) -> Dict[str, Any]:
        """
        Always return the artifacts dictionary.

        This provides direct access to all stored artifacts from pipeline steps.
        """
        if isinstance(self._data, PipelineOutput):
            return self._data.artifacts
        return self._artifacts

    @property
    def structured_data(self) -> PipelineOutput:
        """
        Return PipelineOutput object with both data and artifacts.

        Use this when you need the complete pipeline state for serialization,
        passing to other systems, or when working with structured data flows.
        """
        if isinstance(self._data, PipelineOutput):
            return self._data
        return PipelineOutput(data=self._data, artifacts=self._artifacts.copy())

    @property
    def has_pending_steps(self) -> bool:
        """Check if there are delayed steps waiting to be executed"""
        return len(self.delayed_steps) > 0

    def then(self, func: Callable, *args, **kwargs) -> "Pipeline":
        """Apply a function to the current data (pandas.pipe style)"""
        # Check if there are pending delayed steps
        if self.delayed_steps:
            import warnings

            warnings.warn(
                f"Pipeline has {len(self.delayed_steps)} pending delayed steps. " f"Using then() will execute immediately without running delayed steps first. " f"Consider using execute() to run delayed steps first, or use delayed_then() " f"to add this step to the delayed queue.",
                UserWarning,
                stacklevel=2,
            )

        # Use custom step name if available from decorator
        if hasattr(func, "_step_name"):
            step_name = func._step_name
        else:
            step_name = getattr(func, "__name__", str(func))

        if self._data is None:
            raise ValueError("No data to process. Initialize pipeline with data.")

        # Validate if function is marked as pipeline step
        if hasattr(func, "_is_pipeline_step") and func._is_pipeline_step:
            self.logger.debug(f"Executing pipeline step: {step_name}")

        # Create and execute step
        step = PipelineStep(step_name, func, *args, **kwargs)
        self.steps.append(step)

        self.logger.info(f"Executing step: {step_name}")

        try:
            # Execute step - always pass actual data to function
            result = step.execute(self.data)

            # Update internal data
            self._data = result

            # Store side outputs in artifacts
            if step.side_outputs:
                for key, value in step.side_outputs.items():
                    artifact_key = f"{step_name}.{key}" if key else step_name
                    self._artifacts[artifact_key] = value

            # Store result for history
            self.results.append(self.data)

            self.logger.info(f"Step '{step_name}' completed in {step.execution_time:.3f}s")

            # Log side outputs if any
            if step.side_outputs:
                self.logger.info(f"Step '{step_name}' produced {len(step.side_outputs)} side outputs")

        except Exception as e:
            self.logger.error(f"Step '{step_name}' failed: {str(e)}")
            raise RuntimeError(f"Pipeline failed at step '{step_name}': {str(e)}") from e

        return self

    def delayed_then(self, func: Callable, *args, **kwargs) -> "Pipeline":
        """Add a function to the delayed execution queue without running it immediately"""
        # Use custom step name if available from decorator
        if hasattr(func, "_step_name"):
            step_name = func._step_name
        else:
            step_name = getattr(func, "__name__", str(func))

        # Create delayed step
        delayed_step = DelayedStep(step_name, func, args, kwargs)
        self.delayed_steps.append(delayed_step)

        self.logger.debug(f"Added delayed step: {step_name}")

        return self

    def add_delayed_step(self, func: Callable, index: Optional[int] = None, *args, **kwargs) -> "Pipeline":
        """
        Add a delayed step before a specific position in the delayed execution queue.

        Performs a similar action to the `list.insert()` method.

        Parameters
        ----------
        func : Callable
            Function to add as delayed step
        index : Optional[int]
            Position to insert the step. If None, appends to the end.
            Supports negative indexing.
        *args, **kwargs
            Arguments to pass to the function

        Returns
        -------
        Pipeline
            Self for method chaining

        Examples
        --------
        >>> # Add step at the beginning
        >>> pipeline.add_delayed_step(func1, 0)

        >>> # Add step at the end (same as delayed_then)
        >>> pipeline.add_delayed_step(func2)

        >>> # Insert step at position 2
        >>> pipeline.add_delayed_step(func3, 2)

        >>> # Insert step before the last one
        >>> pipeline.add_delayed_step(func4, -1)
        """
        # Use custom step name if available from decorator
        if hasattr(func, "_step_name"):
            step_name = func._step_name
        else:
            step_name = getattr(func, "__name__", str(func))

        # Create delayed step
        delayed_step = DelayedStep(step_name, func, args, kwargs)

        if index is None:
            # Append to the end (same as delayed_then)
            self.delayed_steps.append(delayed_step)
            self.logger.debug(f"Added delayed step '{step_name}' at end (position {len(self.delayed_steps)-1})")
        else:
            # Insert at specific position
            if index < 0:
                # Handle negative indexing
                actual_index = len(self.delayed_steps) + index
            else:
                actual_index = index

            # Validate index
            if actual_index < 0:
                actual_index = 0
            elif actual_index > len(self.delayed_steps):
                actual_index = len(self.delayed_steps)

            self.delayed_steps.insert(actual_index, delayed_step)
            self.logger.debug(f"Inserted delayed step '{step_name}' at position {actual_index}")

        return self

    def remove_delayed_step(self, index_or_name: Union[int, str]) -> "Pipeline":
        """
        Remove a delayed step at the specified index.

        Parameters
        ----------
        index : int
            Index of the delayed step to remove

        Returns
        -------
        Pepline
            Self for method chaining

        Raises
        ------
        ValueError
            If no delayed step is found with the specified index or name
        """
        if isinstance(index_or_name, int):
            index = index_or_name
        elif isinstance(index_or_name, str):
            # Find index by name
            index = next(
                (i for i, step in enumerate(self.delayed_steps) if step.name == index_or_name),
                None,
            )
            if index is None:
                self.logger.debug(f"Cannot remove delayed step with name '{index_or_name}'. No such step found.")
                raise ValueError(f"Cannot remove delayed step with name '{index_or_name}'. No such step found.")
        else:
            raise TypeError(f"Expect int or str for type(index_or_name), got {type(index_or_name)}")

        if index >= len(self.delayed_steps):
            self.logger.debug(f"Cannot remove delayed step at index {index}. Index out of range.")
            raise ValueError(f"Cannot remove delayed step at index {index}. Index out of range.")

        self.logger.debug(f"Removed delayed step at position {index}")
        del self.delayed_steps[index]
        return self

    def execute(self, steps: Optional[Union[int, List[int]]] = None) -> "Pipeline":
        """
        Execute delayed steps.

        Parameters
        ----------
        steps : Optional[Union[int, List[int]]]
            Which delayed steps to execute:
            - None: execute all delayed steps
            - int: execute the first N delayed steps
            - List[int]: execute specific delayed steps by index

        Returns
        -------
        Pipeline
            Self for method chaining
        """
        if not self.delayed_steps:
            self.logger.info("No delayed steps to execute")
            return self

        if self._data is None:
            raise ValueError("No data to process. Initialize pipeline with data.")

        # Determine which steps to execute
        if steps is None:
            # Execute all delayed steps
            steps_to_execute = self.delayed_steps.copy()
            self.delayed_steps = []
        elif isinstance(steps, int):
            # Execute first N steps
            steps_to_execute = self.delayed_steps[:steps]
            self.delayed_steps = self.delayed_steps[steps:]
        elif isinstance(steps, list):
            # Execute specific steps by index
            steps_to_execute = []
            remaining_steps = []
            for i, delayed_step in enumerate(self.delayed_steps):
                if i in steps:
                    steps_to_execute.append(delayed_step)
                else:
                    remaining_steps.append(delayed_step)
            self.delayed_steps = remaining_steps
        else:
            raise ValueError("steps parameter must be None, int, or List[int]")

        self.logger.info(f"Executing {len(steps_to_execute)} delayed steps")

        # Execute the selected steps
        for delayed_step in steps_to_execute:
            # Convert delayed step to pipeline step and execute
            step = delayed_step.to_pipeline_step()
            self.steps.append(step)

            self.logger.info(f"Executing delayed step: {step.name}")

            try:
                # Execute step - always pass actual data to function
                result = step.execute(self.data)

                # Update internal data
                self._data = result

                # Store side outputs in artifacts
                if step.side_outputs:
                    for key, value in step.side_outputs.items():
                        artifact_key = f"{step.name}.{key}" if key else step.name
                        self._artifacts[artifact_key] = value

                # Store result for history
                self.results.append(self.data)

                self.logger.info(f"Delayed step '{step.name}' completed in {step.execution_time:.3f}s")

                # Log side outputs if any
                if step.side_outputs:
                    self.logger.info(f"Delayed step '{step.name}' produced {len(step.side_outputs)} side outputs")

            except Exception as e:
                self.logger.error(f"Delayed step '{step.name}' failed: {str(e)}")
                raise RuntimeError(f"Pipeline failed at delayed step '{step.name}': {str(e)}") from e

        return self

    def get_delayed_steps_info(self) -> List[Dict[str, Any]]:
        """Get information about delayed steps"""
        return [
            {
                "index": i,
                "name": step.name,
                "function": step.function.__name__,
                "args_count": len(step.args),
                "kwargs_keys": list(step.kwargs.keys()),
                "step_type": getattr(step.function, "_step_type", "unknown"),
                "is_pipeline_step": getattr(step.function, "_is_pipeline_step", False),
            }
            for i, step in enumerate(self.delayed_steps)
        ]

    def apply(self, func: Callable, *args, **kwargs) -> "Pipeline":
        """Apply function and return new Pipeline (functional style)"""
        new_pipeline = self.copy()
        return new_pipeline.then(func, *args, **kwargs)

    def assign(self, **kwargs) -> "Pipeline":
        """Add attributes or computed values to data"""

        def assign_values(data):
            # Handle different data types
            if hasattr(data, "__dict__"):
                for key, value in kwargs.items():
                    setattr(data, key, value)
            elif isinstance(data, dict):
                data.update(kwargs)
            else:
                # For immutable types, wrap in a container
                return {"data": data, **kwargs}
            return data

        return self.then(assign_values)

    def filter(self, condition: Callable) -> "Pipeline":
        """Filter data based on condition"""

        def filter_data(data):
            if hasattr(data, "__iter__") and not isinstance(data, (str, bytes)):
                return type(data)(item for item in data if condition(item))
            elif condition(data):
                return data
            else:
                raise ValueError("Data does not meet filter condition")

        return self.then(filter_data)

    def transform(self, transformer: Callable, *args, **kwargs) -> "Pipeline":
        """Alias of `then`, used to define format transformations."""
        return self.then(transformer, *args, **kwargs)

    def validate(self, validator: Callable, error_msg: str = "Validation failed") -> "Pipeline":
        """Validate data and raise error if invalid"""

        def validate_data(data):
            if not validator(data):
                raise ValueError(error_msg)
            return data

        return self.then(validate_data)

    def peek(self, func: Optional[Callable] = None, prefix: str = "") -> "Pipeline":
        """Inspect data without modifying it (for debugging)"""

        def peek_data(data):
            if func:
                func(data)
            else:
                msg = f"{prefix}Pipeline data: {repr(data)}" if prefix else f"Pipeline data: {repr(data)}"
                self.logger.debug(msg)
            return data

        return self.then(peek_data)

    def store(self, name: str, extractor: Optional[Callable] = None) -> "Pipeline":
        """Store current data or extracted value as artifact"""

        def store_data(data):
            if extractor:
                self._artifacts[name] = extractor(data)
            else:
                self._artifacts[name] = deepcopy(data)
            return data

        return self.then(store_data)

    def copy(self) -> "Pipeline":
        """Create a deep copy of this pipeline"""
        new_pipeline = Pipeline(
            deepcopy(self.data),  # Always copy actual data
            f"{self.name}_copy",
            logging_level=logging.getLevelName(self.logger.level),
        )
        # Copy artifacts and delayed steps
        new_pipeline._artifacts = deepcopy(self.artifacts)
        new_pipeline.delayed_steps = deepcopy(self.delayed_steps)
        return new_pipeline

    def get_data(self) -> Any:
        """
        Get current data (same as .data property).

        Kept for backward compatibility.
        """
        return self.data

    def get_artifact(self, name: str) -> Any:
        """Get a specific artifact by name"""
        if name in self.artifacts:
            return self.artifacts[name]
        else:
            raise KeyError(f"Artifact '{name}' not found. Available: {list(self.artifacts.keys())}")

    def get_all_artifacts(self) -> Dict[str, Any]:
        """Get all stored artifacts"""
        return self.artifacts.copy()

    def get_step_result(self, step_index: Union[int, str]) -> Any:
        """Get result from a specific step by index or name"""
        if isinstance(step_index, int):
            if 0 <= step_index < len(self.results):
                return self.results[step_index]
            else:
                raise IndexError(f"Step index {step_index} out of range")
        else:
            for i, step in enumerate(self.steps):
                if step.name == step_index:
                    if i < len(self.results):
                        return self.results[i]
                    else:
                        raise ValueError(f"Step '{step_index}' has not completed execution")
            raise ValueError(f"Step '{step_index}' not found")

    def get_execution_summary(self) -> Dict[str, Any]:
        """Get summary of pipeline execution"""
        summary = {
            "pipeline_name": self.name,
            "total_steps": len(self.steps),
            "delayed_steps": len(self.delayed_steps),
            "successful_steps": sum(1 for s in self.steps if s.success),
            "failed_steps": sum(1 for s in self.steps if s.error is not None),
            "total_execution_time": sum(s.execution_time or 0 for s in self.steps),
            "artifacts_count": len(self.artifacts),
            "steps": [],
            "delayed_steps_info": self.get_delayed_steps_info(),
        }

        for i, step in enumerate(self.steps):
            step_info = {
                "index": i,
                "name": step.name,
                "success": step.success,
                "execution_time": step.execution_time,
                "error": str(step.error) if step.error else None,
                "side_outputs": (list(step.side_outputs.keys()) if step.side_outputs else []),
            }
            summary["steps"].append(step_info)

        return summary

    def visualize_pipeline(self) -> str:
        """Generate a text visualization of the pipeline"""
        lines = [f"Pipeline: {self.name}", "=" * 40]

        # Show executed steps
        for i, step in enumerate(self.steps):
            status = "✓" if step.success else "✗" if step.error else "○"
            time_str = f"({step.execution_time:.3f}s)" if step.execution_time else "(pending)"

            # Check if it's a decorated pipeline step
            if hasattr(step.function, "_is_pipeline_step"):
                lines.append(f"{status} Step {i+1}: {step.name} {time_str} [validated]")
            else:
                lines.append(f"{status} Step {i+1}: {step.name} {time_str}")

            # Add description if available
            if hasattr(step.function, "_step_description") and step.function._step_description:
                lines.append(f"   └─ {step.function._step_description.strip()}")

            # Show side outputs
            if step.side_outputs:
                for key in step.side_outputs:
                    lines.append(f"   └─ side output: {key}")

        # Show delayed steps
        if self.delayed_steps:
            lines.append("\nDelayed Steps:")
            lines.append("-" * 20)
            for i, delayed_step in enumerate(self.delayed_steps):
                is_pipeline_step = getattr(delayed_step.function, "_is_pipeline_step", False)
                status_str = "[validated]" if is_pipeline_step else ""
                lines.append(f"⏸ Delayed {i+1}: {delayed_step.name} {status_str}")

                # Add description if available
                if hasattr(delayed_step.function, "_step_description") and delayed_step.function._step_description:
                    lines.append(f"   └─ {delayed_step.function._step_description.strip()}")

        lines.append("=" * 40)
        lines.append(f"Current data type: {type(self.data).__name__}")
        lines.append(f"Artifacts stored: {len(self.artifacts)}")
        lines.append(f"Delayed steps: {len(self.delayed_steps)}")

        return "\n".join(lines)

    def save(self, filepath: str, format: str = "pickle") -> "Pipeline":
        """Save current data to file"""
        data_to_save = self.data

        if format == "pickle":
            with open(filepath, "wb") as f:
                pickle.dump(data_to_save, f)
        elif format == "json":
            import json

            with open(filepath, "w") as f:
                json.dump(data_to_save, f)
        else:
            raise ValueError(f"Unsupported format: {format}")

        self.logger.info(f"Pipeline data saved to {filepath}")
        return self

    def save_artifacts(self, filepath: str, format: str = "pickle") -> "Pipeline":
        """Save all artifacts to file"""
        if format == "pickle":
            with open(filepath, "wb") as f:
                pickle.dump(self.artifacts, f)
        elif format == "json":
            import json

            with open(filepath, "w") as f:
                json.dump(self.artifacts, f)
        else:
            raise ValueError(f"Unsupported format: {format}")

        self.logger.info(f"Pipeline artifacts saved to {filepath}")
        return self

    def save_structured_data(self, filepath: str, format: str = "pickle") -> "Pipeline":
        """Save structured data (data + artifacts) to file"""
        structured = self.structured_data

        if format == "pickle":
            with open(filepath, "wb") as f:
                pickle.dump(structured, f)
        elif format == "json":
            import json

            # Convert to dict for JSON serialization
            data_dict = {"data": structured.data, "artifacts": structured.artifacts}
            with open(filepath, "w") as f:
                json.dump(data_dict, f)
        else:
            raise ValueError(f"Unsupported format: {format}")

        self.logger.info(f"Pipeline structured data saved to {filepath}")
        return self

    @classmethod
    def load(cls, filepath: str, format: str = "pickle", name: Optional[str] = None) -> "Pipeline":
        """Load data from file and create new pipeline"""
        if format == "pickle":
            with open(filepath, "rb") as f:
                data = pickle.load(f)
        elif format == "json":
            import json

            with open(filepath, "r") as f:
                data = json.load(f)
        else:
            raise ValueError(f"Unsupported format: {format}")

        return cls(data, name or f"Pipeline_from_{filepath}")

    @classmethod
    def load_structured_data(cls, filepath: str, format: str = "pickle", name: Optional[str] = None) -> "Pipeline":
        """Load structured data from file and create new pipeline"""
        if format == "pickle":
            with open(filepath, "rb") as f:
                structured = pickle.load(f)

            if isinstance(structured, PipelineOutput):
                pipeline = cls(structured.data, name)
                pipeline._artifacts = structured.artifacts
                return pipeline
            else:
                return cls(structured, name)

        elif format == "json":
            import json

            with open(filepath, "r") as f:
                data_dict = json.load(f)

            if isinstance(data_dict, dict) and "data" in data_dict and "artifacts" in data_dict:
                pipeline = cls(data_dict["data"], name)
                pipeline._artifacts = data_dict["artifacts"]
                return pipeline
            else:
                return cls(data_dict, name)
        else:
            raise ValueError(f"Unsupported format: {format}")

    def __str__(self) -> str:
        success_count = sum(1 for s in self.steps if s.success)
        artifacts_str = f", {len(self.artifacts)} artifacts" if self.artifacts else ""
        delayed_str = f", {len(self.delayed_steps)} delayed" if self.delayed_steps else ""
        return f"Pipeline('{self.name}'): {success_count}/{len(self.steps)} steps executed{artifacts_str}{delayed_str}"

    def __repr__(self) -> str:
        return f"<Pipeline name='{self.name}' steps={len(self.steps)} delayed={len(self.delayed_steps)} data_type={type(self.data).__name__} artifacts={len(self.artifacts)}>"


def create_pipeline(data: Any, name: Optional[str] = None, **kwargs) -> Pipeline:
    """Create a new pipeline with initial data"""
    return Pipeline(data, name, **kwargs)


P = ParamSpec("P")
R = TypeVar("R")


@overload
def pipeline_step(_func: Callable[P, R]) -> Callable[P, R]: ...
@overload
def pipeline_step(*, name: Optional[str] = ...) -> Callable[[Callable[P, R]], Callable[P, R]]: ...


def pipeline_step(_func: Optional[Callable[P, R]] = None, *, name: Optional[str] = None) -> Union[Callable[P, R], Callable[[Callable[P, R]], Callable[P, R]]]:
    """
    Decorator for single-output pipeline functions.

    Use this for functions that return a single value (including tuples as single values).
    For multiple outputs, use @multiout_step instead.

    Parameters
    ----------
    name : Optional[str]
        Custom name for the step. If None, uses function name.

    Examples
    --------
    >>> @pipeline_step
    ... def process(data):
    ...     return processed_data  # Single output

    >>> @pipeline_step("process_data")
    ... def process(data):
    ...     return processed_data  # Single output

    >>> @pipeline_step()
    ... def get_coordinates():
    ...     return (10, 20)  # Single tuple output
    """

    def _decorate(func: Callable[P, R]) -> Callable[P, R]:
        step_name = name or func.__name__

        @wraps(func)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
            # place for logging/metrics using `step_name`
            return func(*args, **kwargs)

        # Add consistent metadata attributes
        setattr(wrapper, "_is_pipeline_step", True)
        setattr(wrapper, "_step_type", "single_output")
        setattr(wrapper, "_step_name", step_name)
        setattr(wrapper, "_step_description", func.__doc__)
        setattr(wrapper, "_expected_output_count", 1)
        setattr(wrapper, "_output_names", ["main"])
        return wrapper

    # bare usage: @pipeline_step
    if _func is not None:
        return _decorate(_func)
    # factory usage: @pipeline_step(name="...")
    return _decorate


def multiout_step(
    **outputs: str,
) -> Callable[[Callable[P, Tuple[Any, ...] | Any]], Callable[P, MultiOutput]]:
    """
    Decorator factory for multi-output pipeline functions.

    Parameters
    ----------
    **outputs : str
        Named outputs. Use 'main' to specify which output is the main data flow.
        If 'main' is not specified, the first return value is treated as main.

    Note
    ----
    - Bare usage is NOT allowed: you must call it with parentheses, e.g.
    @multiout_step(stats="statistics", plot="visualization")
    - At least one output must be declared; otherwise use @pipeline_step.

    Examples
    --------
    >>> @multiout_step(stats="statistics", plot="visualization")
    ... def analyze_data(data):
    ...     # return (main, stats, plot)
    ...     return processed, stats_dict, plot_obj
    """
    if not outputs:
        raise TypeError("@multiout_step requires at least one named output; " "for single-output functions, use @pipeline_step.")

    has_explicit_main = "main" in outputs
    side_output_items = [(k, v) for k, v in outputs.items() if k != "main"]
    side_output_names = [v for _, v in side_output_items]
    # If 'main' is provided, expected_count == number of declared outputs.
    # Otherwise, default: first item is main, rest are side outputs.
    expected_count = len(outputs) if has_explicit_main else len(side_output_items) + 1

    def _decorate(func: Callable[P, Tuple[Any, ...] | Any]) -> Callable[P, MultiOutput]:
        @wraps(func)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> MultiOutput:
            results = func(*args, **kwargs)

            # When multiple outputs are declared, a tuple of that length is required.
            if not isinstance(results, tuple):
                if expected_count > 1:
                    raise ValueError(f"Function {func.__name__} decorated with @multiout_step " f"expected {expected_count} return values but got 1 non-tuple value. " f"For single outputs, use @pipeline_step instead.")
                # This branch should be rare because we require at least one side OR explicit main,
                # but keep it for completeness when expected_count == 1.
                return MultiOutput(main=results, side={})

            if len(results) != expected_count:
                raise ValueError(f"Function {func.__name__} decorated with @multiout_step " f"expected {expected_count} return values but got {len(results)}. " f"Declared outputs: {list(outputs.keys())}")

            # Map tuple -> (main, side)
            if has_explicit_main:
                main_index = list(outputs.keys()).index("main")
                main = results[main_index]
                side: Dict[str, Any] = {}
                for i, (_, side_name) in enumerate(side_output_items):
                    # shift index if the side comes after the main slot
                    actual_index = i if i < main_index else i + 1
                    side[side_name] = results[actual_index]
            else:
                main = results[0]
                side = {side_name: results[i + 1] for i, (_, side_name) in enumerate(side_output_items)}

            return MultiOutput(main=main, side=side)

        # Attach metadata for consistency
        setattr(wrapper, "_is_pipeline_step", True)
        setattr(wrapper, "_step_type", "multi_output")
        setattr(wrapper, "_step_name", func.__name__)
        setattr(wrapper, "_step_description", func.__doc__)
        setattr(wrapper, "_expected_output_count", expected_count)
        setattr(wrapper, "_output_names", ["main"] + side_output_names)
        setattr(wrapper, "_declared_outputs", outputs)

        return wrapper

    return _decorate
