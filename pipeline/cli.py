"""
cli.py - entry point.

    python -m pipeline.cli full  --expr counts.csv --meta meta.csv --condition diagnosis --outdir results/
    python -m pipeline.cli deg   --expr counts.csv --meta meta.csv --condition diagnosis --outdir results/deg
    python -m pipeline.cli network --expr results/deg/normalized_expression.csv --deg-file results/deg/deg_significant.csv --outdir results/network
    python -m pipeline.cli enrichment --modules results/network/module_membership.csv --outdir results/enrichment

`full` runs all three stages back to back and writes the combined
module/DEG/pathway correlation table at the end. It's the "adaptive"
entry point - species, data type (counts vs normalized), and covariates
are all detected from the files you give it rather than hardcoded.
"""
import argparse
from pathlib import Path

from . import deg as deg_mod
from . import network as network_mod
from . import enrichment as enrichment_mod
from .utils import correlate_modules_with_degs


def add_common_deg_args(p):
    p.add_argument("--expr", required=True, type=Path, help="expression matrix csv (genes x samples)")
    p.add_argument("--meta", required=True, type=Path, help="sample metadata csv, needs a 'sample' column")
    p.add_argument("--condition", required=True, help="metadata column to test (e.g. diagnosis, group)")
    p.add_argument("--alpha", type=float, default=0.05)


def add_network_args(p):
    p.add_argument("--top-var", type=int, default=3000)
    p.add_argument("--cores", type=int, default=4)
    p.add_argument("--all-genes", action="store_true")


def add_enrichment_args(p):
    p.add_argument("--org-db", default="org.Hs.eg.db", help="Bioconductor OrgDb package, e.g. org.Mm.eg.db for mouse")
    p.add_argument("--key-type", default="SYMBOL")
    p.add_argument("--qval", type=float, default=0.05)


def cmd_deg(args):
    result = deg_mod.run_deg(args.expr, args.meta, args.condition, args.outdir, args.alpha)
    print(f"\nDEG results: {result['results']}")
    print(f"Significant DEGs: {result['significant']}")


def cmd_network(args):
    result = network_mod.run_megena(
        args.expr, args.outdir, deg_file=args.deg_file,
        top_var=args.top_var, cores=args.cores, all_genes=args.all_genes,
    )
    print(f"\nModule membership: {result['modules']}")


def cmd_enrichment(args):
    result = enrichment_mod.run_go_enrichment(
        args.modules, args.outdir, org_db=args.org_db,
        key_type=args.key_type, qval=args.qval, universe=args.universe,
    )
    print(f"\nGO enrichment: {result['go']}")
    print(f"KEGG enrichment: {result['kegg']}")


def cmd_full(args):
    outdir = args.outdir
    deg_dir = outdir / "deg"
    net_dir = outdir / "network"
    go_dir = outdir / "enrichment"

    deg_result = deg_mod.run_deg(args.expr, args.meta, args.condition, deg_dir, args.alpha)

    net_result = network_mod.run_megena(
        deg_result["normalized_expression"], net_dir,
        deg_file=deg_result["significant"],
        top_var=args.top_var, cores=args.cores, all_genes=args.all_genes,
    )

    enrich_result = enrichment_mod.run_go_enrichment(
        net_result["modules"], go_dir, org_db=args.org_db,
        key_type=args.key_type, qval=args.qval, universe=deg_result["results"],
    )

    summary = correlate_modules_with_degs(
        net_result["modules"], deg_result["significant"], enrich_result["go"], outdir,
    )

    print("\n=== pipeline complete ===")
    print(f"DEG results:        {deg_result['results']}")
    print(f"Network modules:    {net_result['modules']}")
    print(f"GO/KEGG enrichment: {enrich_result['go']}")
    print(f"Module/DEG/pathway summary: {outdir / 'module_deg_pathway_summary.csv'}")
    print(f"\nTop modules by DEG enrichment:\n{summary.head(10).to_string(index=False)}")


def main():
    parser = argparse.ArgumentParser(description="Adaptive DEG -> co-expression network -> pathway pipeline")
    sub = parser.add_subparsers(dest="command", required=True)

    p_full = sub.add_parser("full", help="run DEG -> MEGENA -> GO end to end")
    add_common_deg_args(p_full)
    add_network_args(p_full)
    add_enrichment_args(p_full)
    p_full.add_argument("--outdir", required=True, type=Path)
    p_full.set_defaults(func=cmd_full)

    p_deg = sub.add_parser("deg", help="differential expression only")
    add_common_deg_args(p_deg)
    p_deg.add_argument("--outdir", required=True, type=Path)
    p_deg.set_defaults(func=cmd_deg)

    p_net = sub.add_parser("network", help="MEGENA co-expression network only")
    p_net.add_argument("--expr", required=True, type=Path)
    p_net.add_argument("--deg-file", type=Path, default=None)
    add_network_args(p_net)
    p_net.add_argument("--outdir", required=True, type=Path)
    p_net.set_defaults(func=cmd_network)

    p_go = sub.add_parser("enrichment", help="GO/KEGG enrichment on modules only")
    p_go.add_argument("--modules", required=True, type=Path)
    p_go.add_argument("--universe", type=Path, default=None)
    add_enrichment_args(p_go)
    p_go.add_argument("--outdir", required=True, type=Path)
    p_go.set_defaults(func=cmd_enrichment)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

# wip note: document brain_region as alternate --condition example

# wip note: simplify covariate auto-detection from metadata columns

# wip note: fix hub gene extraction per module

# wip note: improve module_membership.csv writer

# wip note: simplify cli full subcommand

# wip note: clarify hypergeometric deg-per-module enrichment

# wip note: refactor deg_results.csv / deg_significant.csv writers

# wip note: improve go (bp) over-representation per module

# wip note: refactor cli full subcommand
