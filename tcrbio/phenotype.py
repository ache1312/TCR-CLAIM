from __future__ import annotations

import numpy as np
import pandas as pd


def clone_phenotype_association(
    tcr: pd.DataFrame,
    *,
    clone_key: str = "ct_strict",
    groupby: list[str] | None = None,
    scores: list[str] | None = None,
    thresholds: list[int] | None = None,
) -> pd.DataFrame:
    """Summarize phenotype score shifts in expanded strict clones vs singletons."""

    groupby = ["dataset_id", "cancer_type", "tissue_type", "cell_class"] if groupby is None else groupby
    scores = [] if scores is None else scores
    thresholds = [2, 5, 10] if thresholds is None else thresholds
    df = pd.DataFrame(tcr).copy()
    df = df[df[clone_key].notna()]
    groupby = [column for column in groupby if column in df.columns]
    rows = []

    if not scores:
        return pd.DataFrame()

    clone_sizes = df.groupby(groupby + [clone_key], dropna=False).size().rename("clone_size").reset_index()
    df = df.merge(clone_sizes, on=groupby + [clone_key], how="left")

    grouped = df.groupby(groupby, dropna=False) if groupby else [((), df)]
    for key, group in grouped:
        context = dict(zip(groupby, key if isinstance(key, tuple) else (key,), strict=False)) if groupby else {}
        singleton = group[group["clone_size"] <= 1]
        for threshold in thresholds:
            expanded = group[group["clone_size"] >= threshold]
            for score in scores:
                if score not in group.columns:
                    continue
                mean_expanded = expanded[score].mean()
                mean_singleton = singleton[score].mean()
                rows.append(
                    {
                        **context,
                        "score": score,
                        "expansion_threshold": threshold,
                        "n_expanded": int(expanded.shape[0]),
                        "n_singleton": int(singleton.shape[0]),
                        "mean_expanded": mean_expanded,
                        "mean_singleton": mean_singleton,
                        "delta_mean": mean_expanded - mean_singleton
                        if not np.isnan(mean_expanded) and not np.isnan(mean_singleton)
                        else np.nan,
                    }
                )
    return pd.DataFrame(rows)
