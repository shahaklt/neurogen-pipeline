"""
utils.py - shared helpers: running Rscript steps and the final
DEG <-> module <-> pathway correlation step.
"""
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import hypergeom

R_DIR = Path(__file__).resolve().parent.parent / "R"


def run_rscript(script_name: str, args: list, step_name: str = None):
    """Run one of our R/*.R scripts, streaming output, raise on failure."""
    if shutil.which("Rscript") is None:
        sys.exit(
            "Rscript not found on PATH. Install R (>= 4.2) and the required "
            "Bioconductor packages - see environment.yml / README."
        )
    script_path = R_DIR / script_name
    cmd = ["Rscript", str(script_path)] + [str(a) for a in args]
    label = step_name or script_name
    print(f"\n[{label}] running: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(f"[{label}] failed (exit {result.returncode}) - see R output above")


def correlate_modules_with_degs(module_csv: Path, deg_csv: Path, go_csv: Path, outdir: Path):
    """
    For each MEGENA module, test whether it's enriched for significant DEGs
    (hypergeometric test against the full gene universe), then attach the
    module's top GO terms so you get one table: module -> DEG enrichment ->
    what pathway that module actually represents.
    """
    modules = pd.read_csv(module_csv)
    degs = pd.read_csv(deg_csv)
    deg_genes = set(degs["gene"])
    all_genes = set(modules["gene"])

    go = pd.read_csv(go_csv) if go_csv.exists() and go_csv.stat().st_size > 0 else pd.DataFrame()

    rows = []
    N = len(all_genes)
    K = len(deg_genes & all_genes)  # DEGs that are actually in the network gene set

    for mod_id, grp in modules.groupby("module"):
        mod_genes = set(grp["gene"])
        n = len(mod_genes)
        k = len(mod_genes & deg_genes)
        # hypergeometric survival function P(X >= k)
        pval = hypergeom.sf(k - 1, N, K, n) if N and K and n else 1.0

        top_terms = ""
        if not go.empty:
            mod_go = go[go["module"] == mod_id].sort_values("p.adjust" if "p.adjust" in go.columns else go.columns[-1])
            top_terms = "; ".join(mod_go["Description"].head(3).tolist()) if "Description" in mod_go.columns else ""

        rows.append({
            "module": mod_id,
            "module_size": n,
            "degs_in_module": k,
            "deg_enrichment_pval": pval,
            "top_GO_terms": top_terms,
        })

    summary = pd.DataFrame(rows).sort_values("deg_enrichment_pval")
    outdir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(outdir / "module_deg_pathway_summary.csv", index=False)
    return summary

# wip note: document mouse org.mm.eg.db enrichment example
