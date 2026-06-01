#!/usr/bin/env python
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run TCR-CLAIM tables on benchmark result directories.")
    parser.add_argument("--results", required=True, help="Comma-separated result directories.")
    parser.add_argument("--out", required=True, help="Output root directory.")
    parser.add_argument("--all-cells", action="store_true", help="Disable primary CD4/CD8 filtering.")
    return parser.parse_args()


def split_csv_arg(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def main() -> None:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    runner = repo_root / "scripts" / "run_tcr_claim_tables.py"
    out_root = Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)

    for result_dir in split_csv_arg(args.results):
        result_path = Path(result_dir)
        cell_path = result_path / "cell_metadata_with_tcr.csv"
        if not cell_path.exists():
            print(f"Skipping {result_path}: missing cell_metadata_with_tcr.csv", file=sys.stderr)
            continue
        out_dir = out_root / result_path.name
        cmd = [
            sys.executable,
            str(runner),
            "--input",
            str(cell_path),
            "--out",
            str(out_dir),
        ]
        if args.all_cells:
            cmd.append("--all-cells")
        subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
