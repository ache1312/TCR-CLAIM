from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from .io import read_table
from .metrics import clone_count_agreement, tissue_sharing
from .pipeline import discover_result_dirs, run_tcr_claim_batch, run_tcr_claim_pipeline


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
    parser.add_argument("--phenotype-scores", default=None, help="Comma-separated phenotype score columns to summarize per candidate.")
    parser.add_argument("--include-xenium-cdr3", action="store_true", help="Add top candidate CDR3 sequences to the Xenium roadmap and export FASTA.")
    parser.add_argument("--max-cdr3-targets", type=int, default=20, help="Maximum CDR3 targets for optional Xenium advanced design.")
    parser.add_argument("--all-cells", action="store_true", help="Disable the primary CD4/CD8 paired-TCR filter.")
    args = parser.parse_args(argv)

    phenotype_scores = _split_csv_arg(args.phenotype_scores) if args.phenotype_scores else None
    run_tcr_claim_pipeline(
        args.input,
        args.out,
        tissue_col=args.tissue_col,
        phenotype_scores=phenotype_scores,
        include_xenium_cdr3=args.include_xenium_cdr3,
        max_cdr3_targets=args.max_cdr3_targets,
        primary_only=not args.all_cells,
    )
    print(f"TCR-CLAIM outputs written to {Path(args.out).resolve()}")


def run_batch_main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Run TCR-CLAIM tables across multiple benchmark result directories.")
    parser.add_argument("--results-root", help="Root directory to search for cell_metadata_with_tcr.csv files.")
    parser.add_argument("--results", help="Comma-separated result directories.")
    parser.add_argument("--out", required=True, help="Output root directory for batch outputs.")
    parser.add_argument("--tissue-col", default="tissue_type")
    parser.add_argument("--phenotype-scores", default=None, help="Comma-separated phenotype score columns to summarize per candidate.")
    parser.add_argument("--include-xenium-cdr3", action="store_true", help="Add top candidate CDR3 sequences to each Xenium roadmap and export FASTA.")
    parser.add_argument("--max-cdr3-targets", type=int, default=20)
    parser.add_argument("--all-cells", action="store_true", help="Disable the primary CD4/CD8 paired-TCR filter.")
    parser.add_argument("--fail-on-error", action="store_true", help="Exit non-zero if any dataset fails.")
    args = parser.parse_args(argv)

    result_dirs = []
    if args.results:
        result_dirs.extend(Path(result) for result in _split_csv_arg(args.results))
    if args.results_root:
        result_dirs.extend(discover_result_dirs(args.results_root))
    if not result_dirs:
        raise SystemExit("Provide --results or --results-root.")

    phenotype_scores = _split_csv_arg(args.phenotype_scores) if args.phenotype_scores else None
    summary = run_tcr_claim_batch(
        result_dirs,
        args.out,
        tissue_col=args.tissue_col,
        phenotype_scores=phenotype_scores,
        include_xenium_cdr3=args.include_xenium_cdr3,
        max_cdr3_targets=args.max_cdr3_targets,
        primary_only=not args.all_cells,
        continue_on_error=True,
    )
    failures = int((summary["status"] == "fail").sum()) if not summary.empty and "status" in summary.columns else 0
    print(f"TCR-CLAIM batch outputs written to {Path(args.out).resolve()}")
    print(f"Datasets processed: {len(summary)}")
    print(f"Dataset failures: {failures}")
    if failures and args.fail_on_error:
        raise SystemExit(1)


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
