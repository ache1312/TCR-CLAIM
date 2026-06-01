from __future__ import annotations

import numpy as np
import pandas as pd

from .filters import primary_tcr_cells
from .metrics import collapse_risk
from .tables import strict_clone_table


def prioritize_candidates(
    tcr: pd.DataFrame,
    *,
    risk: pd.DataFrame | None = None,
    strict_key: str = "ct_strict",
    relaxed_key: str = "ct_vgene",
    primary_only: bool = True,
    focus: str | None = "tumor_cd8",
    top_n: int | None = 100,
) -> pd.DataFrame:
    """Rank strict clone candidates with collapse-risk context.

    The ranking is intentionally heuristic: it is a triage score for biological
    review, not a statistical probability of tumor reactivity.
    """

    df = pd.DataFrame(tcr).copy()
    if primary_only:
        df = primary_tcr_cells(df, strict_key=strict_key, relaxed_key=relaxed_key)
    if risk is None:
        risk = collapse_risk(df, strict_key=strict_key, relaxed_key=relaxed_key, primary_only=False)

    clones = strict_clone_table(
        df,
        clone_key=strict_key,
        relaxed_key=relaxed_key,
        primary_only=False,
    )
    if clones.empty:
        return clones

    risk_cols = [
        "dataset_id",
        "cancer_type",
        "donor_id",
        "sample_id",
        "tissue_type",
        "cell_class",
        relaxed_key,
        "dominant_fraction",
        "n_strict_clones",
        "risk_label",
        "evidence_level",
    ]
    risk_df = pd.DataFrame(risk)
    risk_cols = [column for column in risk_cols if column in risk_df.columns]
    left_keys = [column for column in ["dataset_id", "cancer_type", "donor_id", "sample_id", "tissue_type", "cell_class"] if column in clones.columns and column in risk_df.columns]
    if "relaxed_group" in clones.columns and relaxed_key in risk_df.columns:
        clones = clones.rename(columns={"relaxed_group": relaxed_key})
        left_keys = left_keys + [relaxed_key]

    out = clones.merge(risk_df[risk_cols], on=left_keys, how="left", suffixes=("", "_risk"))
    out["candidate_id"] = out[strict_key]
    out["candidate_type"] = "strict_clone"
    out["focus"] = focus
    out["tumor_context"] = _is_tumor_context(out["tissue_type"]) if "tissue_type" in out.columns else False
    out["is_cd8"] = out["cell_class"].eq("CD8") if "cell_class" in out.columns else False
    out["collapse_risk_penalty"] = out["risk_label"].map({"low": 0.0, "medium": 10.0, "high": 25.0}).fillna(15.0)
    out["dominance_bonus"] = 20.0 * pd.to_numeric(out.get("dominant_fraction"), errors="coerce").fillna(0.0)
    out["focus_bonus"] = 0.0
    if focus == "tumor_cd8":
        out.loc[out["tumor_context"] & out["is_cd8"], "focus_bonus"] = 50.0
    elif focus == "tumor":
        out.loc[out["tumor_context"], "focus_bonus"] = 35.0

    out["rank_score"] = (
        np.log1p(pd.to_numeric(out["n_cells"], errors="coerce").fillna(0.0)) * 20.0
        + out["dominance_bonus"]
        + out["focus_bonus"]
        - out["collapse_risk_penalty"]
    )
    out["allowed_claim"] = "Expanded or recurrent strict paired-CDR3 TCR clone candidate."
    out["not_allowed_claim"] = "Antigen-specific or tumor-reactive clone without orthogonal validation."
    out = out.sort_values(["rank_score", "n_cells"], ascending=[False, False]).reset_index(drop=True)
    out["candidate_rank"] = np.arange(1, len(out) + 1)
    if top_n is not None:
        out = out.head(top_n)
    return out


def _is_tumor_context(series: pd.Series) -> pd.Series:
    tumor_values = {
        "tumor",
        "primary_tumor",
        "metastasis",
        "metastatic_lymph_node",
        "tumor_resident_memory",
        "tumor_recirculating",
        "tumor_cd8_sorted",
        "tumor_pre_treatment",
        "tumor_post_treatment",
    }
    return pd.Series(series).isin(tumor_values)
