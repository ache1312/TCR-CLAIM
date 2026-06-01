from __future__ import annotations

import numpy as np
import pandas as pd

from .filters import primary_tcr_cells


DEFAULT_CANDIDATE_CONTEXT = ["dataset_id", "cancer_type", "donor_id", "sample_id", "tissue_type", "cell_class"]

PHENOTYPE_SCORE_ALIASES = {
    "cytotoxic": [
        "cytotoxic",
        "cytotoxicity",
        "cytotoxic_score",
        "cytotoxicity_score",
        "score_cytotoxic",
        "gzm_prf_score",
    ],
    "exhausted": [
        "exhausted",
        "exhaustion",
        "exhaustion_score",
        "exhausted_score",
        "score_exhaustion",
        "t_cell_exhaustion",
    ],
    "proliferative": [
        "proliferative",
        "proliferation",
        "proliferation_score",
        "cycling",
        "cycling_score",
        "mki67_score",
    ],
    "progenitor_memory": [
        "progenitor",
        "progenitor_score",
        "memory",
        "memory_score",
        "naive_memory",
        "tcf7_score",
    ],
    "treg": [
        "treg",
        "treg_score",
        "suppressive",
        "suppressive_score",
        "foxp3_score",
    ],
    "ifn_response": [
        "ifn",
        "ifn_score",
        "interferon",
        "interferon_score",
        "ifn_response",
        "ifn_response_score",
    ],
}


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


def infer_phenotype_score_columns(tcr: pd.DataFrame, scores: list[str] | None = None) -> list[str]:
    """Detect known numeric phenotype score columns without treating all numeric metadata as phenotype."""

    df = pd.DataFrame(tcr)
    if scores is not None:
        return [score for score in scores if score in df.columns and pd.api.types.is_numeric_dtype(df[score])]

    normalized_to_column = {_normalize_name(column): column for column in df.columns}
    detected = []
    for aliases in PHENOTYPE_SCORE_ALIASES.values():
        for alias in aliases:
            column = normalized_to_column.get(_normalize_name(alias))
            if column is not None and pd.api.types.is_numeric_dtype(df[column]) and column not in detected:
                detected.append(column)
    return detected


def candidate_phenotype_table(
    tcr: pd.DataFrame,
    candidates: pd.DataFrame,
    *,
    scores: list[str] | None = None,
    strict_key: str = "ct_strict",
    relaxed_key: str = "ct_vgene",
    context_cols: list[str] | None = None,
    primary_only: bool = True,
    min_delta: float = 0.1,
) -> pd.DataFrame:
    """Summarize available phenotype scores for prioritized strict clone candidates.

    The output is descriptive. It reports candidate-vs-context differences and
    intentionally does not infer antigen specificity, tumor reactivity, or
    functional activity.
    """

    context_cols = DEFAULT_CANDIDATE_CONTEXT if context_cols is None else context_cols
    df = pd.DataFrame(tcr).copy()
    candidate_df = pd.DataFrame(candidates).copy()
    if primary_only:
        df = primary_tcr_cells(df, strict_key=strict_key, relaxed_key=relaxed_key)
    if candidate_df.empty:
        return pd.DataFrame()

    score_cols = infer_phenotype_score_columns(df, scores=scores)
    rows = []
    for _, candidate in candidate_df.iterrows():
        candidate_id = candidate.get("candidate_id", candidate.get(strict_key))
        context = _candidate_context(candidate, context_cols)
        context_mask = _context_mask(df, context)
        context_df = df[context_mask]
        clone_df = context_df[context_df[strict_key] == candidate_id] if strict_key in context_df.columns else context_df.iloc[0:0]
        reference_df = context_df[context_df[strict_key] != candidate_id] if strict_key in context_df.columns else context_df
        base = {
            **context,
            "candidate_rank": candidate.get("candidate_rank", pd.NA),
            "candidate_id": candidate_id,
            "n_candidate_cells": int(len(clone_df)),
            "n_reference_cells": int(len(reference_df)),
            "rank_score": candidate.get("rank_score", pd.NA),
            "collapse_risk_label": candidate.get("risk_label", pd.NA),
            "collapse_dominant_fraction": candidate.get("dominant_fraction", pd.NA),
        }
        if not score_cols:
            rows.append(
                {
                    **base,
                    "phenotype_state": pd.NA,
                    "score_column": pd.NA,
                    "candidate_mean": pd.NA,
                    "reference_mean": pd.NA,
                    "delta_mean": pd.NA,
                    "candidate_median": pd.NA,
                    "reference_median": pd.NA,
                    "direction": "not_evaluable",
                    "phenotype_evidence_status": "no_scores_available",
                    "allowed_interpretation": "Clone can be prioritized by TCR evidence, but no phenotype score was available.",
                    "not_allowed_interpretation": "Phenotype-associated clone.",
                }
            )
            continue

        for score in score_cols:
            candidate_values = pd.to_numeric(clone_df[score], errors="coerce").dropna()
            reference_values = pd.to_numeric(reference_df[score], errors="coerce").dropna()
            candidate_mean = candidate_values.mean() if not candidate_values.empty else np.nan
            reference_mean = reference_values.mean() if not reference_values.empty else np.nan
            delta = candidate_mean - reference_mean if not np.isnan(candidate_mean) and not np.isnan(reference_mean) else np.nan
            state = phenotype_state_for_score(score)
            rows.append(
                {
                    **base,
                    "phenotype_state": state,
                    "score_column": score,
                    "candidate_mean": candidate_mean,
                    "reference_mean": reference_mean,
                    "delta_mean": delta,
                    "candidate_median": candidate_values.median() if not candidate_values.empty else np.nan,
                    "reference_median": reference_values.median() if not reference_values.empty else np.nan,
                    "direction": _direction(delta, min_delta=min_delta),
                    "phenotype_evidence_status": "scored",
                    "allowed_interpretation": f"Candidate has descriptive {state} score context.",
                    "not_allowed_interpretation": "Antigen specificity or tumor reactivity.",
                }
            )
    return pd.DataFrame(rows)


def phenotype_state_for_score(score: str) -> str:
    normalized = _normalize_name(score)
    for state, aliases in PHENOTYPE_SCORE_ALIASES.items():
        if normalized in {_normalize_name(alias) for alias in aliases}:
            return state
    return normalized


def _candidate_context(candidate: pd.Series, context_cols: list[str]) -> dict[str, object]:
    context = {}
    for column in context_cols:
        if column in candidate.index and pd.notna(candidate[column]):
            context[column] = candidate[column]
    return context


def _context_mask(df: pd.DataFrame, context: dict[str, object]) -> pd.Series:
    mask = pd.Series(True, index=df.index)
    for column, value in context.items():
        if column in df.columns:
            mask &= df[column] == value
    return mask


def _direction(delta: float, *, min_delta: float) -> str:
    if np.isnan(delta):
        return "not_evaluable"
    if delta >= min_delta:
        return "enriched"
    if delta <= -min_delta:
        return "depleted"
    return "no_clear_shift"


def _normalize_name(value: str) -> str:
    return str(value).strip().lower().replace("-", "_").replace(" ", "_")
