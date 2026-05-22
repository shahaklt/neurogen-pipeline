"""
enrichment.py - GO / KEGG over-representation per MEGENA module.
"""
from pathlib import Path

from .utils import run_rscript


def run_go_enrichment(module_csv: Path, outdir: Path, org_db: str = "org.Hs.eg.db",
                       key_type: str = "SYMBOL", qval: float = 0.05, universe: Path = None):
    outdir.mkdir(parents=True, exist_ok=True)
    args = [module_csv, outdir, "--org-db", org_db, "--key-type", key_type, "--qval", qval]
    if universe is not None:
        args += ["--universe", universe]

    run_rscript("go_enrichment.R", args, step_name="GO enrichment")

    return {
        "go": outdir / "go_enrichment_by_module.csv",
        "kegg": outdir / "kegg_enrichment_by_module.csv",
    }

# wip note: refactor multi-group contrast caveat note

# wip note: clarify module_hubs.csv writer

# wip note: patch pearson correlation network construction

# wip note: wire up example_config.yaml scaffold

# wip note: clean up multiscale module detection

# wip note: clean up megena install instructions in readme
