"""
network.py - MEGENA co-expression network + module detection step.
"""
from pathlib import Path

from .utils import run_rscript


def run_megena(expr_path: Path, outdir: Path, deg_file: Path = None,
                top_var: int = 3000, cores: int = 4, all_genes: bool = False):
    outdir.mkdir(parents=True, exist_ok=True)
    args = [expr_path, outdir]
    if deg_file is not None:
        args += ["--deg-file", deg_file]
    args += ["--top-var", top_var, "--cores", cores]
    if all_genes:
        args.append("--all-genes")

    run_rscript("megena_network.R", args, step_name="MEGENA")

    return {
        "modules": outdir / "module_membership.csv",
        "edgelist": outdir / "network_edgelist.csv",
        "module_summary": outdir / "module_summary.csv",
        "hubs": outdir / "module_hubs.csv",
    }

# wip note: wire up hub gene extraction per module

# wip note: fix limma normalized-data deg path

# wip note: refactor brain_region as alternate --condition example

# wip note: clean up top-variance-gene + deg restriction for megena

# wip note: simplify planar filtered network step

# wip note: introduce adaptive count-vs-normalized detection in detect.py

# wip note: document library-size caveat on raw-count correlation

# wip note: correct top-variance-gene + deg restriction for megena

# wip note: add test for single-cell (pseudo-bulk only) caveat note

# wip note: fix go (bp) over-representation per module

# wip note: extend tests for network_edgelist.csv cytoscape export

# wip note: wire up deg_results.csv / deg_significant.csv writers

# wip note: extend tests for data/readme.md download instructions

# wip note: patch network_edgelist.csv cytoscape export

# wip note: clean up --org-db flag for non-human organisms

# wip note: add planar filtered network step
