from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd


DEFAULT_CONTEXT = ["dataset_id", "cancer_type", "donor_id", "sample_id", "tissue_type", "cell_class"]
DEFAULT_SHARING_CONTEXT = ["dataset_id", "cancer_type", "donor_id", "cell_class"]


def _present(df: pd.DataFrame, columns: list[str]) -> list[str]:
    return [column for column in columns if column in df.columns]


def _iter_groups(df: pd.DataFrame, columns: list[str]):
    columns = _present(df, columns)
    if not columns:
        yield {}, df
        return
    for key, group in df.groupby(columns, dropna=False):
        if len(columns) == 1:
            key = (key,)
        yield dict(zip(columns, key, strict=False)), group


def _entropy(counts: pd.Series) -> float:
    values = counts[counts > 0].to_numpy(dtype=float)
    if len(values) <= 1:
        return 0.0
    p = values / values.sum()
    return float(-(p * np.log(p)).sum())


def collapse_risk(
    tcr: pd.DataFrame,
    *,
    relaxed_key: str = "ct_vgene",
    strict_key: str = "ct_strict",
    groupby: list[str] | None = None,
    low_risk_fraction: float = 0.90,
    medium_risk_fraction: float = 0.70,
) -> pd.DataFrame:
    """Quantify how many strict CDR3 clones are collapsed by each V-gene group."""

    groupby = DEFAULT_CONTEXT if groupby is None else groupby
    df = pd.DataFrame(tcr).copy()
    keep = df[relaxed_key].notna() & df[strict_key].notna()
    df = df[keep]
    rows = []

    for context, group in _iter_groups(df, groupby + [relaxed_key]):
        counts = group[strict_key].value_counts(dropna=True)
        total_cells = int(counts.sum())
        n_strict = int(len(counts))
        dominant_clone = counts.index[0] if n_strict else pd.NA
        dominant_cells = int(counts.iloc[0]) if n_strict else 0
        dominant_fraction = dominant_cells / total_cells if total_cells else np.nan
        risk_label = (
            "not_evaluable"
            if not total_cells
            else "low"
            if n_strict == 1 or dominant_fraction >= low_risk_fraction
            else "medium"
            if dominant_fraction >= medium_risk_fraction
            else "high"
        )
        rows.append(
            {
                **context,
                "n_cells": total_cells,
                "n_strict_clones": n_strict,
                "dominant_strict_clone": dominant_clone,
                "dominant_strict_cells": dominant_cells,
                "dominant_fraction": dominant_fraction,
                "cdr3_entropy": _entropy(counts),
                "is_dominated": bool(total_cells and dominant_fraction >= low_risk_fraction),
                "is_mixed": bool(n_strict > 1),
                "risk_label": risk_label,
                "evidence_level": "dominant_vgene_group" if risk_label == "low" else "mixed_vgene_group",
            }
        )

    return pd.DataFrame(rows)


def _count_clones(df: pd.DataFrame, clone_key: str, min_clone_size: int) -> int:
    counts = df[clone_key].dropna().value_counts()
    return int((counts >= min_clone_size).sum())


def clone_count_agreement(
    tcr: pd.DataFrame,
    *,
    strict_key: str = "ct_strict",
    relaxed_key: str = "ct_vgene",
    groupby: list[str] | None = None,
    thresholds: list[int] | None = None,
    include_top10: bool = True,
) -> pd.DataFrame:
    """Compare strict clone counts with relaxed V-gene group counts."""

    groupby = DEFAULT_CONTEXT if groupby is None else groupby
    thresholds = [1, 2, 5, 10] if thresholds is None else thresholds
    df = pd.DataFrame(tcr).copy()
    df = df[df[strict_key].notna() & df[relaxed_key].notna()]
    rows = []

    for context, group in _iter_groups(df, groupby):
        for threshold in thresholds:
            strict_n = _count_clones(group, strict_key, threshold)
            relaxed_n = _count_clones(group, relaxed_key, threshold)
            rows.append(
                _agreement_row(context, f"size_ge{threshold}", threshold, len(group), strict_n, relaxed_n)
            )
        if include_top10:
            top = group[strict_key].value_counts().head(10).index
            top_group = group[group[strict_key].isin(top)]
            rows.append(
                _agreement_row(
                    context,
                    "top10_strict",
                    pd.NA,
                    len(top_group),
                    top_group[strict_key].nunique(),
                    top_group[relaxed_key].nunique(),
                )
            )
    return pd.DataFrame(rows)


def _agreement_row(context, analysis_set, min_clone_size, n_cells, strict_n, relaxed_n):
    absolute = int(relaxed_n - strict_n)
    relative = absolute / strict_n if strict_n else np.nan
    ratio = relaxed_n / strict_n if strict_n else np.nan
    return {
        **context,
        "analysis_set": analysis_set,
        "min_clone_size": min_clone_size,
        "n_cells": int(n_cells),
        "strict_clone_count": int(strict_n),
        "vgene_clone_count": int(relaxed_n),
        "absolute_difference": absolute,
        "relative_difference": relative,
        "clone_count_ratio": ratio,
    }


def _clone_presence(pair_df: pd.DataFrame, clone_key: str, tissue_col: str, min_clone_size: int) -> pd.DataFrame:
    sizes = pair_df[clone_key].dropna().value_counts().rename("clone_size")
    presence = pair_df[pair_df[clone_key].notna()][[tissue_col, clone_key]].drop_duplicates()
    presence = presence.rename(columns={clone_key: "clone"})
    presence = presence.merge(sizes, left_on="clone", right_index=True, how="left")
    presence = presence[presence["clone_size"] >= min_clone_size]
    if presence.empty:
        return pd.DataFrame(columns=["clone", "n_tissues", "clone_size", "shared"])
    out = (
        presence.groupby("clone", dropna=False)
        .agg(n_tissues=(tissue_col, "nunique"), clone_size=("clone_size", "first"))
        .reset_index()
    )
    out["shared"] = out["n_tissues"] >= 2
    return out


def _shared_with_clone_column(presence: pd.DataFrame, clone_key: str) -> pd.DataFrame:
    if presence.empty:
        return pd.DataFrame(columns=[clone_key, "n_tissues", "clone_size", "shared"])
    shared = presence[presence["shared"]].copy()
    if "clone" not in shared.columns:
        return pd.DataFrame(columns=[clone_key, "n_tissues", "clone_size", "shared"])
    return shared.rename(columns={"clone": clone_key})


def tissue_sharing(
    tcr: pd.DataFrame,
    *,
    strict_key: str = "ct_strict",
    relaxed_key: str = "ct_vgene",
    tissue_col: str = "tissue_type",
    groupby: list[str] | None = None,
    thresholds: list[int] | None = None,
) -> pd.DataFrame:
    """Separate strict sharing from relaxed apparent-only sharing."""

    groupby = DEFAULT_SHARING_CONTEXT if groupby is None else groupby
    thresholds = [1, 2, 5, 10] if thresholds is None else thresholds
    df = pd.DataFrame(tcr).copy()
    df = df[df[strict_key].notna() & df[relaxed_key].notna() & df[tissue_col].notna()]
    rows = []

    for context, group in _iter_groups(df, groupby):
        tissues = sorted(group[tissue_col].dropna().unique())
        for tissue_a, tissue_b in combinations(tissues, 2):
            pair_df = group[group[tissue_col].isin([tissue_a, tissue_b])].copy()
            strict_map = pair_df[[strict_key, relaxed_key]].dropna().drop_duplicates()
            for threshold in thresholds:
                strict_presence = _clone_presence(pair_df, strict_key, tissue_col, threshold)
                relaxed_presence = _clone_presence(pair_df, relaxed_key, tissue_col, threshold)
                strict_shared = _shared_with_clone_column(strict_presence, strict_key)
                relaxed_shared = _shared_with_clone_column(relaxed_presence, relaxed_key)

                strict_shared_vgenes = strict_shared[[strict_key]].merge(strict_map, on=strict_key, how="left")
                strict_shared_vgenes = set(strict_shared_vgenes[relaxed_key].dropna())
                relaxed_shared_keys = set(relaxed_shared[relaxed_key].dropna())
                strict_recovered = strict_shared[[strict_key]].merge(strict_map, on=strict_key, how="left")
                strict_recovered_n = int(strict_recovered[relaxed_key].isin(relaxed_shared_keys).sum())
                backed_n = int(relaxed_shared[relaxed_key].isin(strict_shared_vgenes).sum())
                strict_shared_n = int(len(strict_shared))
                relaxed_shared_n = int(len(relaxed_shared))

                rows.append(
                    {
                        **context,
                        "tissue_a": tissue_a,
                        "tissue_b": tissue_b,
                        "analysis_set": f"size_ge{threshold}",
                        "min_clone_size": threshold,
                        "strict_union_clones": int(len(strict_presence)),
                        "vgene_union_groups": int(len(relaxed_presence)),
                        "strict_shared_clones": strict_shared_n,
                        "vgene_shared_groups": relaxed_shared_n,
                        "vgene_shared_strict_backed": backed_n,
                        "vgene_shared_apparent_only": relaxed_shared_n - backed_n,
                        "strict_shared_recovered_by_vgene": strict_recovered_n,
                        "vgene_sharing_precision_like": backed_n / relaxed_shared_n if relaxed_shared_n else np.nan,
                        "strict_sharing_recall_by_vgene": strict_recovered_n / strict_shared_n
                        if strict_shared_n
                        else np.nan,
                    }
                )
    return pd.DataFrame(rows)
