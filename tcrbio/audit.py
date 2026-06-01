from __future__ import annotations

from pathlib import Path

import pandas as pd

from .pipeline import discover_result_dirs


KEY_COLUMNS = ["cell_class", "ct_strict", "ct_vgene"]


def audit_result_dirs(
    result_dirs: list[str | Path],
    out_dir: str | Path,
    *,
    chunksize: int = 250_000,
    allowed_cell_classes: tuple[str, ...] = ("CD4", "CD8"),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Audit benchmark result directories without loading large tables at once."""

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    cell_class_rows = []
    for result_dir in result_dirs:
        result_path = Path(result_dir)
        input_path = result_path / "cell_metadata_with_tcr.csv"
        if not input_path.exists():
            rows.append(
                {
                    "result_id": result_path.name,
                    "status": "missing_input",
                    "input_path": str(input_path),
                }
            )
            continue
        summary, cell_counts = audit_cell_metadata(
            input_path,
            result_id=result_path.name,
            chunksize=chunksize,
            allowed_cell_classes=allowed_cell_classes,
        )
        rows.append(summary)
        cell_class_rows.extend(cell_counts.to_dict("records"))

    summary_df = pd.DataFrame(rows)
    cell_class_df = pd.DataFrame(cell_class_rows)
    summary_df.to_csv(out / "dataset_audit_summary.csv", index=False)
    cell_class_df.to_csv(out / "dataset_cell_class_counts.csv", index=False)
    write_audit_report(summary_df, cell_class_df, out / "dataset_audit_report.md")
    return summary_df, cell_class_df


def audit_results_root(
    results_root: str | Path,
    out_dir: str | Path,
    *,
    chunksize: int = 250_000,
    allowed_cell_classes: tuple[str, ...] = ("CD4", "CD8"),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Discover and audit all result directories below a root."""

    return audit_result_dirs(
        discover_result_dirs(results_root),
        out_dir,
        chunksize=chunksize,
        allowed_cell_classes=allowed_cell_classes,
    )


def audit_cell_metadata(
    input_path: str | Path,
    *,
    result_id: str,
    chunksize: int = 250_000,
    allowed_cell_classes: tuple[str, ...] = ("CD4", "CD8"),
) -> tuple[dict[str, object], pd.DataFrame]:
    """Audit one cell_metadata_with_tcr.csv table in chunks."""

    input_path = Path(input_path)
    header = pd.read_csv(input_path, nrows=0, low_memory=False)
    columns = list(header.columns)
    usecols = [column for column in KEY_COLUMNS if column in columns]
    cell_counts: dict[str, int] = {}
    n_rows = 0
    n_strict = 0
    n_vgene = 0
    n_paired = 0
    n_primary = 0
    for chunk in pd.read_csv(input_path, usecols=usecols or None, chunksize=chunksize, low_memory=False):
        n_rows += len(chunk)
        strict = _present(chunk["ct_strict"]) if "ct_strict" in chunk.columns else pd.Series(False, index=chunk.index)
        vgene = _present(chunk["ct_vgene"]) if "ct_vgene" in chunk.columns else pd.Series(False, index=chunk.index)
        paired = strict & vgene
        n_strict += int(strict.sum())
        n_vgene += int(vgene.sum())
        n_paired += int(paired.sum())
        if "cell_class" in chunk.columns:
            classes = chunk["cell_class"].fillna("<NA>").astype(str)
            counts = classes.value_counts(dropna=False)
            for label, count in counts.items():
                cell_counts[label] = cell_counts.get(label, 0) + int(count)
            primary = paired & classes.isin(allowed_cell_classes)
        else:
            primary = paired
        n_primary += int(primary.sum())

    status = "ok"
    reason = ""
    if n_rows == 0:
        status = "empty_input"
        reason = "No rows in cell metadata table."
    elif n_paired == 0:
        status = "no_paired_tcr"
        reason = "No rows with both ct_strict and ct_vgene."
    elif n_primary == 0:
        status = "zero_primary_universe"
        reason = "Rows have TCR keys, but none pass the CD4/CD8 primary universe."

    summary = {
        "result_id": result_id,
        "status": status,
        "reason": reason,
        "input_path": str(input_path),
        "n_rows": n_rows,
        "n_columns": len(columns),
        "has_cell_class": "cell_class" in columns,
        "has_ct_strict": "ct_strict" in columns,
        "has_ct_vgene": "ct_vgene" in columns,
        "n_ct_strict": n_strict,
        "n_ct_vgene": n_vgene,
        "n_paired_tcr": n_paired,
        "n_primary_cd4_cd8_paired": n_primary,
        "paired_tcr_fraction": _ratio(n_paired, n_rows),
        "primary_cd4_cd8_paired_fraction": _ratio(n_primary, n_rows),
    }
    cell_class_df = (
        pd.DataFrame(
            [
                {"result_id": result_id, "cell_class": label, "n_cells": count}
                for label, count in sorted(cell_counts.items(), key=lambda item: item[1], reverse=True)
            ]
        )
        if cell_counts
        else pd.DataFrame(columns=["result_id", "cell_class", "n_cells"])
    )
    return summary, cell_class_df


def write_audit_report(summary: pd.DataFrame, cell_class_counts: pd.DataFrame, output: str | Path) -> Path:
    """Write a Markdown report for dataset inclusion/exclusion review."""

    output = Path(output)
    lines = [
        "# TCR-CLAIM Dataset Audit",
        "",
        "This audit explains which benchmark datasets enter the primary CD4/CD8 paired-TCR universe and which require special handling.",
        "",
        "## Status Counts",
        "",
        _markdown_table(summary["status"].value_counts(dropna=False).rename_axis("status").reset_index(name="n_datasets"))
        if not summary.empty and "status" in summary.columns
        else "_No datasets._",
        "",
        "## Dataset Summary",
        "",
        _markdown_table(
            summary[
                [
                    column
                    for column in [
                        "result_id",
                        "status",
                        "n_rows",
                        "n_paired_tcr",
                        "n_primary_cd4_cd8_paired",
                        "paired_tcr_fraction",
                        "primary_cd4_cd8_paired_fraction",
                        "reason",
                    ]
                    if column in summary.columns
                ]
            ]
        ),
    ]
    if not cell_class_counts.empty:
        top = (
            cell_class_counts.sort_values(["result_id", "n_cells"], ascending=[True, False])
            .groupby("result_id", as_index=False)
            .head(8)
        )
        lines.extend(["", "## Top Cell-Class Labels", "", _markdown_table(top)])
    output.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    return output


def _present(series: pd.Series) -> pd.Series:
    return series.notna() & series.astype(str).str.strip().ne("")


def _ratio(numerator: int, denominator: int) -> float:
    if denominator == 0:
        return float("nan")
    return numerator / denominator


def _markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    out = pd.DataFrame(df).copy()
    for column in out.columns:
        if pd.api.types.is_float_dtype(out[column]):
            out[column] = out[column].map(lambda value: "" if pd.isna(value) else f"{value:.4f}")
        else:
            out[column] = out[column].map(lambda value: "" if pd.isna(value) else str(value))
    header = "| " + " | ".join(out.columns) + " |"
    separator = "| " + " | ".join(["---"] * len(out.columns)) + " |"
    rows = ["| " + " | ".join(row) + " |" for row in out.astype(str).to_numpy()]
    return "\n".join([header, separator, *rows])
