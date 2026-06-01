from __future__ import annotations

import numpy as np
import pandas as pd

from .filters import primary_tcr_cells


DEFAULT_DIVERSITY_CONTEXT = ["dataset_id", "cancer_type", "tissue_type", "cell_class"]


def clonal_diversity(
    tcr: pd.DataFrame,
    *,
    clone_key: str,
    groupby: list[str] | None = None,
    strict_key: str = "ct_strict",
    relaxed_key: str = "ct_vgene",
    primary_only: bool = True,
) -> pd.DataFrame:
    """Compute repertoire diversity metrics for a clone definition."""

    groupby = DEFAULT_DIVERSITY_CONTEXT if groupby is None else groupby
    df = pd.DataFrame(tcr).copy()
    if primary_only:
        df = primary_tcr_cells(df, strict_key=strict_key, relaxed_key=relaxed_key)
    df = df[df[clone_key].notna()]
    groupby = [column for column in groupby if column in df.columns]
    grouped = df.groupby(groupby, dropna=False) if groupby else [((), df)]

    rows = []
    for key, group in grouped:
        context = dict(zip(groupby, key if isinstance(key, tuple) else (key,), strict=False)) if groupby else {}
        counts = group[clone_key].value_counts()
        rows.append({**context, **_diversity_from_counts(counts)})
    return pd.DataFrame(rows)


def strict_vs_relaxed_diversity(
    tcr: pd.DataFrame,
    *,
    strict_key: str = "ct_strict",
    relaxed_key: str = "ct_vgene",
    groupby: list[str] | None = None,
    primary_only: bool = True,
) -> pd.DataFrame:
    """Compare diversity metrics between strict clones and relaxed V-gene groups."""

    groupby = DEFAULT_DIVERSITY_CONTEXT if groupby is None else groupby
    strict = clonal_diversity(
        tcr,
        clone_key=strict_key,
        groupby=groupby,
        strict_key=strict_key,
        relaxed_key=relaxed_key,
        primary_only=primary_only,
    ).add_suffix("_strict")
    relaxed = clonal_diversity(
        tcr,
        clone_key=relaxed_key,
        groupby=groupby,
        strict_key=strict_key,
        relaxed_key=relaxed_key,
        primary_only=primary_only,
    ).add_suffix("_relaxed")
    if strict.empty or relaxed.empty:
        return pd.DataFrame()

    left_keys = [f"{column}_strict" for column in groupby if f"{column}_strict" in strict.columns]
    right_keys = [f"{column}_relaxed" for column in groupby if f"{column}_relaxed" in relaxed.columns]
    merged = strict.merge(relaxed, left_on=left_keys, right_on=right_keys, how="outer")
    for left_key, right_key, original in zip(left_keys, right_keys, groupby, strict=False):
        merged[original] = merged[left_key].combine_first(merged[right_key])

    metrics = ["richness", "effective_shannon", "inverse_simpson", "clonality", "max_clone_fraction", "gini"]
    for metric in metrics:
        s_col = f"{metric}_strict"
        r_col = f"{metric}_relaxed"
        if s_col not in merged.columns or r_col not in merged.columns:
            continue
        merged[f"{metric}_difference"] = merged[r_col] - merged[s_col]
        merged[f"{metric}_ratio"] = merged[r_col] / merged[s_col]
        merged.loc[merged[s_col] == 0, f"{metric}_ratio"] = np.nan
        merged[f"{metric}_relative_difference"] = merged[f"{metric}_difference"] / merged[s_col]
        merged.loc[merged[s_col] == 0, f"{metric}_relative_difference"] = np.nan

    ordered = [column for column in groupby if column in merged.columns]
    remaining = [column for column in merged.columns if column not in ordered and not column.endswith("_strict") and not column.endswith("_relaxed")]
    kept_strict = [column for column in merged.columns if column.endswith("_strict") and column.removesuffix("_strict") not in groupby]
    kept_relaxed = [column for column in merged.columns if column.endswith("_relaxed") and column.removesuffix("_relaxed") not in groupby]
    return merged[ordered + kept_relaxed + kept_strict + remaining]


def _diversity_from_counts(counts: pd.Series) -> dict[str, float | int]:
    values = counts[counts > 0].to_numpy(dtype=float)
    total = float(values.sum())
    richness = int(len(values))
    if total == 0 or richness == 0:
        return {
            "n_cells": 0,
            "richness": 0,
            "singleton_clones": 0,
            "shannon": 0.0,
            "effective_shannon": 0.0,
            "simpson_diversity": 0.0,
            "inverse_simpson": 0.0,
            "pielou_evenness": 0.0,
            "clonality": 0.0,
            "max_clone_fraction": np.nan,
            "gini": 0.0,
        }
    p = values / total
    shannon = float(-(p * np.log(p)).sum()) if richness > 1 else 0.0
    simpson = float(1.0 - np.square(p).sum())
    inverse_simpson = float(1.0 / np.square(p).sum()) if np.square(p).sum() else np.nan
    pielou = float(shannon / np.log(richness)) if richness > 1 else 1.0
    return {
        "n_cells": int(total),
        "richness": richness,
        "singleton_clones": int((values == 1).sum()),
        "shannon": shannon,
        "effective_shannon": float(np.exp(shannon)),
        "simpson_diversity": simpson,
        "inverse_simpson": inverse_simpson,
        "pielou_evenness": pielou,
        "clonality": float(1.0 - pielou),
        "max_clone_fraction": float(values.max() / total),
        "gini": _gini(values),
    }


def _gini(values: np.ndarray) -> float:
    values = np.sort(values[values >= 0])
    n = len(values)
    if n <= 1 or values.sum() == 0:
        return 0.0
    return float((2 * np.sum(np.arange(1, n + 1) * values) / (n * values.sum())) - ((n + 1) / n))
