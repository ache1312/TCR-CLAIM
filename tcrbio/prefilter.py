from __future__ import annotations

from pathlib import Path

import pandas as pd


def prefilter_cell_metadata(
    input_path: str | Path,
    output_path: str | Path,
    *,
    chunksize: int = 250_000,
    strict_key: str = "ct_strict",
    relaxed_key: str = "ct_vgene",
    cell_class_col: str = "cell_class",
    allowed_cell_classes: tuple[str, ...] = ("CD4", "CD8"),
    require_paired_tcr: bool = True,
) -> dict[str, object]:
    """Write a filtered cell metadata table without loading the full input at once."""

    input_path = Path(input_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    stats_path = output_path.with_suffix(".prefilter_summary.csv")

    header = pd.read_csv(input_path, nrows=0, low_memory=False)
    has_strict = strict_key in header.columns
    has_relaxed = relaxed_key in header.columns
    has_cell_class = cell_class_col in header.columns

    n_rows = 0
    n_paired = 0
    n_class_allowed = 0
    n_written = 0
    wrote_header = False
    for chunk in pd.read_csv(input_path, chunksize=chunksize, low_memory=False):
        n_rows += len(chunk)
        keep = pd.Series(True, index=chunk.index)
        if require_paired_tcr:
            strict = _present(chunk[strict_key]) if has_strict else pd.Series(False, index=chunk.index)
            relaxed = _present(chunk[relaxed_key]) if has_relaxed else pd.Series(False, index=chunk.index)
            paired = strict & relaxed
            n_paired += int(paired.sum())
            keep &= paired
        if allowed_cell_classes and has_cell_class:
            allowed = chunk[cell_class_col].fillna("").astype(str).isin(allowed_cell_classes)
            n_class_allowed += int(allowed.sum())
            keep &= allowed
        filtered = chunk[keep]
        n_written += len(filtered)
        filtered.to_csv(output_path, mode="a", index=False, header=not wrote_header)
        wrote_header = True

    if not wrote_header:
        header.to_csv(output_path, index=False)

    summary = {
        "input_path": str(input_path),
        "output_path": str(output_path),
        "n_input_rows": n_rows,
        "n_paired_tcr": n_paired,
        "n_allowed_cell_class": n_class_allowed if has_cell_class else pd.NA,
        "n_output_rows": n_written,
        "retained_fraction": n_written / n_rows if n_rows else float("nan"),
        "require_paired_tcr": require_paired_tcr,
        "allowed_cell_classes": ",".join(allowed_cell_classes),
        "has_ct_strict": has_strict,
        "has_ct_vgene": has_relaxed,
        "has_cell_class": has_cell_class,
    }
    pd.DataFrame([summary]).to_csv(stats_path, index=False)
    return summary


def _present(series: pd.Series) -> pd.Series:
    return series.notna() & series.astype(str).str.strip().ne("")
