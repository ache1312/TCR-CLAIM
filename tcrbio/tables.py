from __future__ import annotations

import pandas as pd

from .filters import primary_tcr_cells


DEFAULT_CONTEXT = ["dataset_id", "cancer_type", "donor_id", "sample_id", "tissue_type", "cell_class"]


def strict_clone_table(
    tcr: pd.DataFrame,
    *,
    clone_key: str = "ct_strict",
    relaxed_key: str = "ct_vgene",
    groupby: list[str] | None = None,
    primary_only: bool = True,
) -> pd.DataFrame:
    """Summarize one row per strict clone and biological context."""

    groupby = DEFAULT_CONTEXT if groupby is None else groupby
    df = pd.DataFrame(tcr).copy()
    if primary_only:
        df = primary_tcr_cells(df, strict_key=clone_key, relaxed_key=relaxed_key)
    groupby = [column for column in groupby if column in df.columns]
    required = groupby + [clone_key]
    df = df[df[clone_key].notna()]

    aggregations = {"n_cells": (clone_key, "size")}
    if relaxed_key in df.columns:
        aggregations["n_relaxed_groups"] = (relaxed_key, "nunique")
        aggregations["relaxed_group"] = (relaxed_key, lambda x: sorted(set(x.dropna()))[0] if x.dropna().any() else pd.NA)
    for column in ["trav", "traj", "tra_cdr3_nt", "trbv", "trbj", "trb_cdr3_nt"]:
        if column in df.columns:
            aggregations[column] = (column, lambda x: x.dropna().iloc[0] if len(x.dropna()) else pd.NA)

    out = df.groupby(required, dropna=False).agg(**aggregations).reset_index()
    out["evidence_level"] = "confirmed_clone"
    return out


def cell_tcr_table(tcr: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of the canonical per-cell TCR table."""

    return pd.DataFrame(tcr).copy()
