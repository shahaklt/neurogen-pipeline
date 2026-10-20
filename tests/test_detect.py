import numpy as np
import pandas as pd
import pytest

from pipeline.detect import profile_expression_matrix, validate_metadata


def test_raw_counts_detected():
    rng = np.random.default_rng(0)
    expr = pd.DataFrame(
        rng.integers(0, 5000, size=(50, 8)),
        columns=[f"s{i}" for i in range(8)],
    )
    profile = profile_expression_matrix(expr)
    assert profile.data_type == "raw_counts"
    assert profile.deg_method == "deseq2"


def test_normalized_data_detected():
    rng = np.random.default_rng(0)
    expr = pd.DataFrame(
        rng.normal(0, 1, size=(50, 8)),
        columns=[f"s{i}" for i in range(8)],
    )
    profile = profile_expression_matrix(expr)
    assert profile.data_type == "normalized"
    assert profile.deg_method == "limma"


def test_log_scale_data_flagged_as_normalized():
    rng = np.random.default_rng(0)
    expr = pd.DataFrame(
        rng.uniform(2, 14, size=(50, 8)),
        columns=[f"s{i}" for i in range(8)],
    )
    profile = profile_expression_matrix(expr)
    assert profile.data_type == "normalized"
    assert profile.looks_log_scale
    assert profile.deg_method == "limma"


def test_missing_values_ignored_in_detection():
    rng = np.random.default_rng(0)
    values = rng.integers(0, 5000, size=(50, 8)).astype(float)
    values[0, 0] = np.nan
    expr = pd.DataFrame(values, columns=[f"s{i}" for i in range(8)])
    profile = profile_expression_matrix(expr)
    assert profile.data_type == "raw_counts"


def test_validate_metadata_requires_sample_overlap():
    meta = pd.DataFrame({"sample": list("abcde"), "diagnosis": ["OUD", "OUD", "control", "control", "control"]})
    with pytest.raises(ValueError, match="only 2 samples overlap"):
        validate_metadata(meta, "diagnosis", ["a", "b", "x", "y", "z"])


def test_validate_metadata_requires_sample_column():
    meta = pd.DataFrame({"id": ["a", "b"], "diagnosis": ["OUD", "control"]})
    with pytest.raises(ValueError, match="sample"):
        validate_metadata(meta, "diagnosis", ["a", "b"])


def test_validate_metadata_requires_two_groups():
    meta = pd.DataFrame({"sample": ["a", "b", "c", "d"], "diagnosis": ["OUD"] * 4})
    with pytest.raises(ValueError, match="fewer than 2 groups"):
        validate_metadata(meta, "diagnosis", ["a", "b", "c", "d"])
