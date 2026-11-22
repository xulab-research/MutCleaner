import pytest

from mutcleaner.core.pipeline import Pipeline, pipeline_step


def failing_step(data):
    raise ValueError("intentional failure")


def append_two_step(data):
    return data + [2]


def append_three_step(data):
    return data + [3]


def double_step(data):
    return [value * 2 for value in data]


def add_one_step(data):
    return [value + 1 for value in data]


class TestDelayedPipelineFailureState:
    """Delayed steps and step results must stay consistent after a failed step"""

    def test_failed_step_keeps_unexecuted_steps_in_delayed_queue(self):
        pipeline = Pipeline([1])
        pipeline.delayed_then(failing_step).delayed_then(append_two_step)

        with pytest.raises(RuntimeError, match="Pipeline failed at delayed step 'failing_step'"):
            pipeline.execute()

        assert [step.name for step in pipeline.delayed_steps] == ["append_two_step"]
        assert pipeline.has_pending_steps is True
        assert [(step.name, step.success) for step in pipeline.steps] == [("failing_step", False)]

    def test_failed_step_is_consumed_from_delayed_queue(self):
        pipeline = Pipeline([1])
        pipeline.delayed_then(failing_step).delayed_then(append_two_step)

        with pytest.raises(RuntimeError):
            pipeline.execute()

        assert "failing_step" not in [step.name for step in pipeline.delayed_steps]
        assert len(pipeline.steps) == 1

    def test_resume_after_failure_runs_queued_step(self):
        pipeline = Pipeline([1])
        pipeline.delayed_then(failing_step).delayed_then(append_two_step)

        with pytest.raises(RuntimeError):
            pipeline.execute()

        pipeline.execute()

        assert pipeline.data == [1, 2]
        assert [(step.name, step.success) for step in pipeline.steps] == [("failing_step", False), ("append_two_step", True)]
        assert pipeline.results == [None, [1, 2]]
        assert len(pipeline.delayed_steps) == 0

    def test_step_result_does_not_mix_after_failure_and_resume(self):
        pipeline = Pipeline([1])
        pipeline.delayed_then(failing_step).delayed_then(append_two_step)

        with pytest.raises(RuntimeError):
            pipeline.execute()
        pipeline.execute()

        assert pipeline.get_step_result(0) is None
        assert pipeline.get_step_result(1) == [1, 2]
        assert pipeline.get_step_result("failing_step") is None
        assert pipeline.get_step_result("append_two_step") == [1, 2]

    def test_results_stay_aligned_with_steps_on_failure(self):
        pipeline = Pipeline([1])
        pipeline.delayed_then(failing_step).delayed_then(append_two_step)

        with pytest.raises(RuntimeError):
            pipeline.execute()

        assert len(pipeline.results) == len(pipeline.steps) == 1
        assert pipeline.results[0] is None

    def test_multiple_remaining_steps_resume_in_order(self):
        pipeline = Pipeline([1])
        pipeline.delayed_then(failing_step).delayed_then(append_two_step).delayed_then(append_three_step)

        with pytest.raises(RuntimeError):
            pipeline.execute()

        assert [step.name for step in pipeline.delayed_steps] == ["append_two_step", "append_three_step"]

        pipeline.execute()

        assert pipeline.data == [1, 2, 3]
        assert pipeline.results == [None, [1, 2], [1, 2, 3]]
        assert [step.success for step in pipeline.steps] == [False, True, True]

    def test_execution_summary_after_failure_and_resume(self):
        pipeline = Pipeline([1])
        pipeline.delayed_then(failing_step).delayed_then(append_two_step)

        with pytest.raises(RuntimeError):
            pipeline.execute()

        summary = pipeline.get_execution_summary()
        assert summary["total_steps"] == 1
        assert summary["successful_steps"] == 0
        assert summary["failed_steps"] == 1
        assert summary["delayed_steps"] == 1
        assert summary["steps"][0]["name"] == "failing_step"
        assert summary["steps"][0]["success"] is False

        pipeline.execute()

        summary = pipeline.get_execution_summary()
        assert summary["total_steps"] == 2
        assert summary["successful_steps"] == 1
        assert summary["failed_steps"] == 1
        assert summary["delayed_steps"] == 0
        assert [step_info["success"] for step_info in summary["steps"]] == [False, True]


class TestPipelineSuccessPathsUnchanged:
    """Immediate and delayed success paths keep their previous behavior"""

    def test_immediate_execution_chain_unchanged(self):
        pipeline = Pipeline([1, 2, 3])
        pipeline.then(double_step).then(add_one_step)

        assert pipeline.data == [3, 5, 7]
        assert pipeline.results == [[2, 4, 6], [3, 5, 7]]
        assert pipeline.get_step_result(0) == [2, 4, 6]
        assert pipeline.get_step_result(1) == [3, 5, 7]
        assert pipeline.get_step_result("double_step") == [2, 4, 6]
        assert pipeline.get_step_result("add_one_step") == [3, 5, 7]

    def test_delayed_execution_success_unchanged(self):
        pipeline = Pipeline([1, 2, 3])
        pipeline.delayed_then(double_step).delayed_then(add_one_step)
        pipeline.execute()

        assert pipeline.data == [3, 5, 7]
        assert len(pipeline.delayed_steps) == 0
        assert pipeline.results == [[2, 4, 6], [3, 5, 7]]
        assert pipeline.get_step_result(0) == [2, 4, 6]
        assert pipeline.get_step_result(1) == [3, 5, 7]
        assert pipeline.get_step_result("double_step") == [2, 4, 6]
        assert pipeline.get_step_result("add_one_step") == [3, 5, 7]

    def test_execute_selected_indices_success_unchanged(self):
        pipeline = Pipeline([1, 2, 3])
        pipeline.delayed_then(double_step).delayed_then(add_one_step).delayed_then(append_two_step)

        pipeline.execute([0, 2])

        assert pipeline.data == [2, 4, 6, 2]
        assert [step.name for step in pipeline.delayed_steps] == ["add_one_step"]
        assert pipeline.get_step_result(0) == [2, 4, 6]
        assert pipeline.get_step_result(1) == [2, 4, 6, 2]

    def test_partial_execute_then_resume_success_unchanged(self):
        pipeline = Pipeline([1, 2, 3])
        pipeline.delayed_then(double_step).delayed_then(add_one_step)

        pipeline.execute(1)
        assert pipeline.data == [2, 4, 6]
        assert [step.name for step in pipeline.delayed_steps] == ["add_one_step"]

        pipeline.execute()
        assert pipeline.data == [3, 5, 7]
        assert len(pipeline.delayed_steps) == 0
        assert pipeline.results == [[2, 4, 6], [3, 5, 7]]


class TestFailedImmediateStepResult:
    """A failed immediate step returns None instead of borrowing a later result"""

    def test_failed_then_step_result_is_none(self):
        pipeline = Pipeline([1])
        pipeline.then(double_step)

        with pytest.raises(RuntimeError, match="Pipeline failed at step 'failing_step'"):
            pipeline.then(failing_step)

        pipeline.then(add_one_step)

        assert pipeline.data == [3]
        assert pipeline.results == [[2], None, [3]]
        assert pipeline.get_step_result(0) == [2]
        assert pipeline.get_step_result(1) is None
        assert pipeline.get_step_result(2) == [3]
        assert pipeline.get_step_result("failing_step") is None
        assert pipeline.get_step_result("add_one_step") == [3]


@pipeline_step(name="decorated_step")
def decorated_step(data):
    return data


def test_name_lookup_uses_matching_step_after_failure():
    pipeline = Pipeline([1])

    @pipeline_step(name="first_step")
    def first_step(data):
        return [2]

    pipeline.then(first_step)
    with pytest.raises(RuntimeError):
        pipeline.then(failing_step)
    pipeline.then(decorated_step)

    assert pipeline.get_step_result("first_step") == [2]
    assert pipeline.get_step_result("failing_step") is None
    assert pipeline.get_step_result("decorated_step") == [2]
