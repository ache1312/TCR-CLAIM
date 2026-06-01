from __future__ import annotations

import pandas as pd


DEFAULT_QC_CONTEXT = ["dataset_id", "cancer_type", "donor_id", "sample_id", "tissue_type"]


def qc_summary(
    tcr: pd.DataFrame,
    *,
    groupby: list[str] | None = None,
    strict_key: str = "ct_strict",
    relaxed_key: str = "ct_vgene",
    cell_class_col: str = "cell_class",
) -> pd.DataFrame:
    """Summarize TCR-CLAIM eligibility and chain ambiguity per context."""

    groupby = DEFAULT_QC_CONTEXT if groupby is None else groupby
    df = pd.DataFrame(tcr).copy()
    groupby = [column for column in groupby if column in df.columns]
    if "has_multi_tra" not in df.columns:
        df["has_multi_tra"] = pd.NA
    if "has_multi_trb" not in df.columns:
        df["has_multi_trb"] = pd.NA
    if cell_class_col not in df.columns:
        df[cell_class_col] = pd.NA

    paired = df[strict_key].notna() & df[relaxed_key].notna()
    primary = paired & df[cell_class_col].isin(["CD4", "CD8"])

    df = df.assign(_paired_tcr=paired, _primary_cd4_cd8_paired=primary)
    grouped = df.groupby(groupby, dropna=False) if groupby else [((), df)]
    rows = []
    for key, group in grouped:
        context = dict(zip(groupby, key if isinstance(key, tuple) else (key,), strict=False)) if groupby else {}
        rows.append(
            {
                **context,
                "n_cells": int(len(group)),
                "n_cd4": int((group[cell_class_col] == "CD4").sum()),
                "n_cd8": int((group[cell_class_col] == "CD8").sum()),
                "n_paired_tcr": int(group["_paired_tcr"].sum()),
                "n_primary_cd4_cd8_paired": int(group["_primary_cd4_cd8_paired"].sum()),
                "paired_tcr_fraction": float(group["_paired_tcr"].mean()) if len(group) else pd.NA,
                "primary_cd4_cd8_paired_fraction": float(group["_primary_cd4_cd8_paired"].mean()) if len(group) else pd.NA,
                "multi_tra_rate": _safe_bool_mean(group["has_multi_tra"]),
                "multi_trb_rate": _safe_bool_mean(group["has_multi_trb"]),
            }
        )
    return pd.DataFrame(rows)


def _safe_bool_mean(series: pd.Series):
    valid = series.dropna()
    if valid.empty:
        return pd.NA
    return float((valid == True).mean())
