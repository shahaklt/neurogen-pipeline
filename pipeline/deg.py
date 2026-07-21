"""
deg.py - differential expression step. Picks DESeq2 or limma based on
what detect.py found, runs it, hands back paths to the results.
"""
from pathlib import Path

import pandas as pd

from .detect import profile_expression_matrix, validate_metadata
from .utils import run_rscript


def run_deg(expr_path: Path, meta_path: Path, condition_col: str, outdir: Path, alpha: float = 0.05):
    expr = pd.read_csv(expr_path, index_col=0)
    meta = pd.read_csv(meta_path)

    profile = profile_expression_matrix(expr)
    warnings = validate_metadata(meta, condition_col, expr.columns)

    print("Dataset profile:")
    for note in profile.notes:
        print(f"  - {note}")
    for w in warnings:
        print(f"  ! {w}")

    outdir.mkdir(parents=True, exist_ok=True)
    script = "deg_deseq2.R" if profile.deg_method == "deseq2" else "deg_limma.R"
    run_rscript(script, [expr_path, meta_path, condition_col, outdir, alpha], step_name="DEG")

    return {
        "method": profile.deg_method,
        "profile": profile,
        "results": outdir / "deg_results.csv",
        "significant": outdir / "deg_significant.csv",
        "normalized_expression": outdir / "normalized_expression.csv",
    }

# wip note: add deseq2 raw-count deg path

# wip note: correct hypergeometric deg-per-module enrichment

# wip note: simplify vst() variance-stabilizing transform for network input

# wip note: tighten multi-group contrast caveat note

# wip note: simplify adaptive count-vs-normalized detection in detect.py

# wip note: add meta.csv sample-column validation

# wip note: cover edge case in limma normalized-data deg path

# wip note: add pearson correlation network construction

# wip note: clean up repo layout section in readme

# wip note: add cli enrichment subcommand

# wip note: clarify library-size caveat on raw-count correlation

# wip note: patch multiscale module detection

# wip note: clean up library-size caveat on raw-count correlation

# wip note: correct normalized_expression.csv writer
