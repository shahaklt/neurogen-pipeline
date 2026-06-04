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


def test_validate_metadata_requires_sample_column():
    meta = pd.DataFrame({"id": ["a", "b"], "diagnosis": ["OUD", "control"]})
    with pytest.raises(ValueError, match="sample"):
        validate_metadata(meta, "diagnosis", ["a", "b"])


def test_validate_metadata_requires_two_groups():
    meta = pd.DataFrame({"sample": ["a", "b", "c", "d"], "diagnosis": ["OUD"] * 4})
    with pytest.raises(ValueError, match="fewer than 2 groups"):
        validate_metadata(meta, "diagnosis", ["a", "b", "c", "d"])

# wip note: improve covariate auto-detection from metadata columns

# wip note: correct mouse org.mm.eg.db enrichment example

# wip note: clean up cli deg subcommand

# wip note: wire up module_hubs.csv writer

# wip note: improve library-size caveat on raw-count correlation

# wip note: fix cli network subcommand

# wip note: refactor unit tests for detect.py edge cases

# wip note: introduce vst() variance-stabilizing transform for network input

# wip note: document limma normalized-data deg path

# wip note: document batch-correction caveat note

# wip note: note unit tests for detect.py edge cases
