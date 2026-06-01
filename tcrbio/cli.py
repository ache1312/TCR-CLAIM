from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from . import (
    cell_tcr_table,
    claim_checker,
    clonal_diversity,
    clone_count_agreement,
    collapse_risk,
    define_clones,
    design_xenium_panel_from_candidates,
    prioritize_candidates,
    qc_summary,
    read_table,
    report_clone_cards,
    report_markdown_summary,
    strict_clone_table,
    strict_vs_relaxed_diversity,
    tissue_sharing,
)


AGREEMENT_KEYS = [
    "dataset_id",
    "cancer_type",
    "donor_id",
    "sample_id",
    "tissue_type",
    "cell_class",
    "analysis_set",
]
SHARING_KEYS = ["dataset_id", "cancer_type", "donor_id", "cell_class", "tissue_a", "tissue_b", "analysis_set"]


def run_tables_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Run TCR-CLAIM tables from cell-level TCR metadata.")
    parser.add_argument("--input", required=True, help="CSV/TSV cell metadata with ct_strict and ct_vgene.")
    parser.add_argument("--out", required=True, help="Output directory.")
    parser.add_argument("--tissue-col", default="tissue_type")
    parser.add_argument("--all-cells", action="store_true", help="Disable the primary CD4/CD8 paired-TCR filter.")
    args = parser.parse_args(argv)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    tcr = read_table(args.input)
    if "ct_strict" not in tcr.columns or "ct_vgene" not in tcr.columns:
        tcr = define_clones(tcr)

    primary_only = not args.all_cells
    risk = collapse_risk(tcr, primary_only=primary_only)
    sharing = tissue_sharing(tcr, tissue_col=args.tissue_col, primary_only=primary_only)
    qc = qc_summary(tcr)
    strict_diversity = clonal_diversity(tcr, clone_key="ct_strict", primary_only=primary_only)
    relaxed_diversity = clonal_diversity(tcr, clone_key="ct_vgene", primary_only=primary_only)
    diversity_comparison = strict_vs_relaxed_diversity(tcr, primary_only=primary_only)
    candidates = prioritize_candidates(tcr, risk=risk, primary_only=primary_only)
    claims = claim_checker(tcr=tcr, risk=risk, sharing=sharing, candidates=candidates)

    cell_tcr_table(tcr).to_csv(out / "cell_tcr_table.csv", index=False)
    qc.to_csv(out / "qc_summary.csv", index=False)
    strict_clone_table(tcr, primary_only=primary_only).to_csv(out / "strict_clone_table.csv", index=False)
    risk.to_csv(out / "vgene_group_table.csv", index=False)
    clone_count_agreement(tcr, primary_only=primary_only).to_csv(out / "clone_count_agreement.csv", index=False)
    strict_diversity.to_csv(out / "strict_clonal_diversity.csv", index=False)
    relaxed_diversity.to_csv(out / "relaxed_clonal_diversity.csv", index=False)
    diversity_comparison.to_csv(out / "strict_vs_relaxed_diversity.csv", index=False)
    sharing.to_csv(out / "sharing_table.csv", index=False)
    candidates.to_csv(out / "candidate_table.csv", index=False)
    claims.to_csv(out / "claim_table.csv", index=False)
    design_xenium_panel_from_candidates().to_csv(out / "xenium_panel_roadmap.csv", index=False)
    report_clone_cards(claims, output=out / "clone_cards.html")
    report_markdown_summary(
        output=out / "tcr_claim_report.md",
        qc=qc,
        diversity=diversity_comparison,
        candidates=candidates,
        sharing=sharing,
        claims=claims,
        risk=risk,
    )
    print(f"TCR-CLAIM outputs written to {out.resolve()}")


def validate_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Validate Python TCR-CLAIM metrics against R benchmark outputs.")
    parser.add_argument("--results-dir", required=True, help="Benchmark result directory with cell_metadata_with_tcr.csv.")
    parser.add_argument("--out", required=True, help="Output CSV validation summary path.")
    parser.add_argument("--atol", type=float, default=1e-9)
    parser.add_argument("--rtol", type=float, default=1e-9)
    args = parser.parse_args(argv)

    summary = validate_result_dir(Path(args.results_dir), atol=args.atol, rtol=args.rtol)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out, index=False)
    failures = int((summary["status"] == "fail").sum()) if not summary.empty else 0
    print(f"Wrote validation summary to {out.resolve()}")
    print(f"Validation failures: {failures}")
    if failures:
        raise SystemExit(1)


def validate_batch_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Validate multiple R benchmark result directories.")
    parser.add_argument("--results", required=True, help="Comma-separated result directories.")
    parser.add_argument("--out", required=True, help="Output CSV validation summary path.")
    parser.add_argument("--atol", type=float, default=1e-9)
    parser.add_argument("--rtol", type=float, default=1e-9)
    args = parser.parse_args(argv)

    rows = []
    for result_dir in _split_csv_arg(args.results):
        result_path = Path(result_dir)
        if not (result_path / "cell_metadata_with_tcr.csv").exists():
            rows.append(
                {
                    "result_id": result_path.name,
                    "table": "_input",
                    "column": "cell_metadata_with_tcr.csv",
                    "n_expected_rows": 0,
                    "n_observed_rows": 0,
                    "n_merged_rows": 0,
                    "n_mismatches": 1,
                    "status": "fail",
                }
            )
            continue
        summary = validate_result_dir(result_path, atol=args.atol, rtol=args.rtol)
        summary.insert(0, "result_id", result_path.name)
        rows.extend(summary.to_dict("records"))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    batch = pd.DataFrame(rows)
    batch.to_csv(out, index=False)
    failures = int((batch["status"] == "fail").sum()) if not batch.empty else 0
    print(f"Wrote batch validation summary to {out.resolve()}")
    print(f"Validation failures: {failures}")
    if failures:
        raise SystemExit(1)


def validate_result_dir(results_dir: Path, *, atol: float = 1e-9, rtol: float = 1e-9) -> pd.DataFrame:
    cell_path = results_dir / "cell_metadata_with_tcr.csv"
    if not cell_path.exists():
        raise SystemExit(f"Missing {cell_path}")

    tcr = read_table(cell_path)
    rows: list[dict[str, object]] = []

    agreement_path = results_dir / "clone_count_agreement.csv"
    if agreement_path.exists():
        rows.extend(
            compare_tables(
                read_table(agreement_path),
                clone_count_agreement(tcr),
                keys=AGREEMENT_KEYS,
                table_name="clone_count_agreement",
                atol=atol,
                rtol=rtol,
            )
        )

    sharing_path = results_dir / "tissue_sharing.csv"
    if sharing_path.exists():
        rows.extend(
            compare_tables(
                read_table(sharing_path),
                tissue_sharing(tcr),
                keys=SHARING_KEYS,
                table_name="tissue_sharing",
                atol=atol,
                rtol=rtol,
            )
        )
    return pd.DataFrame(rows)


def _split_csv_arg(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def compare_tables(
    expected: pd.DataFrame,
    observed: pd.DataFrame,
    *,
    keys: list[str],
    table_name: str,
    atol: float,
    rtol: float,
) -> list[dict[str, object]]:
    common = [column for column in expected.columns if column in observed.columns]
    keys = [key for key in keys if key in common]
    merged = expected[common].merge(
        observed[common],
        on=keys,
        how="outer",
        suffixes=("_expected", "_observed"),
        indicator=True,
    )
    rows = [
        {
            "table": table_name,
            "column": "_merge",
            "n_expected_rows": len(expected),
            "n_observed_rows": len(observed),
            "n_merged_rows": len(merged),
            "n_mismatches": int((merged["_merge"] != "both").sum()),
            "status": "pass" if (merged["_merge"] == "both").all() else "fail",
        }
    ]
    for column in common:
        if column in keys:
            continue
        expected_col = f"{column}_expected"
        observed_col = f"{column}_observed"
        if expected_col not in merged.columns or observed_col not in merged.columns:
            continue
        left = merged[expected_col]
        right = merged[observed_col]
        if pd.api.types.is_numeric_dtype(left) or pd.api.types.is_numeric_dtype(right):
            equal = np.isclose(
                pd.to_numeric(left, errors="coerce"),
                pd.to_numeric(right, errors="coerce"),
                equal_nan=True,
                atol=atol,
                rtol=rtol,
            )
        else:
            equal = left.fillna("<NA>").astype(str) == right.fillna("<NA>").astype(str)
        rows.append(
            {
                "table": table_name,
                "column": column,
                "n_expected_rows": len(expected),
                "n_observed_rows": len(observed),
                "n_merged_rows": len(merged),
                "n_mismatches": int((~equal).sum()),
                "status": "pass" if bool(equal.all()) else "fail",
            }
        )
    return rows
