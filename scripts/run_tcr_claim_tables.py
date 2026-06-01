#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import tcrbio as tb


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run TCR-CLAIM tables from cell-level TCR metadata.")
    parser.add_argument("--input", required=True, help="CSV/TSV cell metadata with ct_strict and ct_vgene.")
    parser.add_argument("--out", required=True, help="Output directory.")
    parser.add_argument("--tissue-col", default="tissue_type")
    return parser.parse_args()


def read_input(path: str) -> pd.DataFrame:
    return tb.read_table(path)


def main() -> None:
    args = parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    tcr = read_input(args.input)
    if "ct_strict" not in tcr.columns or "ct_vgene" not in tcr.columns:
        tcr = tb.define_clones(tcr)

    cell_table = tb.cell_tcr_table(tcr)
    strict_table = tb.strict_clone_table(tcr)
    risk = tb.collapse_risk(tcr)
    agreement = tb.clone_count_agreement(tcr)
    sharing = tb.tissue_sharing(tcr, tissue_col=args.tissue_col)
    claims = tb.claim_checker(tcr=tcr, risk=risk, sharing=sharing)
    panel = tb.design_xenium_panel_from_candidates()

    cell_table.to_csv(out / "cell_tcr_table.csv", index=False)
    strict_table.to_csv(out / "strict_clone_table.csv", index=False)
    risk.to_csv(out / "vgene_group_table.csv", index=False)
    agreement.to_csv(out / "clone_count_agreement.csv", index=False)
    sharing.to_csv(out / "sharing_table.csv", index=False)
    claims.to_csv(out / "claim_table.csv", index=False)
    panel.to_csv(out / "xenium_panel_roadmap.csv", index=False)
    tb.report_clone_cards(claims, output=out / "clone_cards.html")
    print(f"TCR-CLAIM outputs written to {out.resolve()}")


if __name__ == "__main__":
    main()
