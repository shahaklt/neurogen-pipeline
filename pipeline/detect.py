"""
detect.py

The "morphing" part of the pipeline: look at whatever expression matrix
got handed to us and decide how it should be treated downstream, instead
of making the user tell us. This is what lets one CLI handle raw RNA-seq
counts, already-normalized microarray intensities, or god knows what else
someone exports from a processing tool.

Heuristics are deliberately simple and conservative - when in doubt we
fall back to the safer/more general path (limma) rather than guessing
DESeq2 and blowing up on non-integer input.
"""
from dataclasses import dataclass
import numpy as np
import pandas as pd


@dataclass
class DatasetProfile:
    data_type: str          # "raw_counts" | "normalized"
    deg_method: str          # "deseq2" | "limma"
    n_genes: int
    n_samples: int
    has_negative_values: bool
    looks_log_scale: bool
    notes: list


def profile_expression_matrix(expr: pd.DataFrame) -> DatasetProfile:
    notes = []
    values = expr.to_numpy(dtype=float, copy=False)
    non_na = values[~np.isnan(values)]

    has_negative = bool((non_na < 0).any())
    # raw counts are integers (or very close to it, e.g. kallisto/salmon
    # estimated counts) and non-negative
    frac_non_integer = np.mean(np.abs(non_na - np.round(non_na)) > 1e-6) if len(non_na) else 1.0
    is_integer_like = frac_non_integer < 0.01 and not has_negative

    # log-scale data (log2 CPM, microarray intensities) rarely exceeds ~20-25
    # while raw counts routinely run into the thousands
    looks_log_scale = bool(non_na.size and np.percentile(non_na, 99) < 30 and not is_integer_like)

    if is_integer_like:
        data_type = "raw_counts"
        deg_method = "deseq2"
        notes.append("values are non-negative integers -> treating as raw RNA-seq counts, using DESeq2")
    else:
        data_type = "normalized"
        deg_method = "limma"
        reason = "contains negative values" if has_negative else "non-integer values"
        notes.append(f"{reason} -> treating as pre-normalized data (microarray/logCPM/etc), using limma")

    if looks_log_scale:
        notes.append("value range is consistent with log-scale data")

    return DatasetProfile(
        data_type=data_type,
        deg_method=deg_method,
        n_genes=expr.shape[0],
        n_samples=expr.shape[1],
        has_negative_values=has_negative,
        looks_log_scale=looks_log_scale,
        notes=notes,
    )


def validate_metadata(meta: pd.DataFrame, condition_col: str, expr_columns) -> list:
    """Sanity-check metadata against the expression matrix, return list of warnings."""
    warnings = []
    if "sample" not in meta.columns:
        raise ValueError('metadata file must have a "sample" column matching the expression matrix column names')
    if condition_col not in meta.columns:
        raise ValueError(f'condition column "{condition_col}" not found in metadata columns: {list(meta.columns)}')

    overlap = set(meta["sample"]) & set(expr_columns)
    if len(overlap) < 4:
        raise ValueError(
            f"only {len(overlap)} samples overlap between expression matrix and metadata - "
            "check that sample ids match exactly"
        )
    if len(overlap) < len(expr_columns):
        warnings.append(f"{len(expr_columns) - len(overlap)} expression columns have no metadata row and will be dropped")

    n_levels = meta.loc[meta["sample"].isin(overlap), condition_col].nunique()
    if n_levels < 2:
        raise ValueError(f'condition column "{condition_col}" has fewer than 2 groups among overlapping samples')
    if n_levels > 2:
        warnings.append(
            f'condition column "{condition_col}" has {n_levels} levels - DESeq2/limma will use '
            "the reference-level contrast; consider a dedicated multi-group design if that's not what you want"
        )
    return warnings

# wip note: introduce brain_region as alternate --condition example

# wip note: add rscript subprocess runner in utils.py

# wip note: correct multiscale module detection

# wip note: introduce cli deg subcommand
