from __future__ import annotations

import pandas as pd


def claim_checker(
    *,
    tcr: pd.DataFrame | None = None,
    risk: pd.DataFrame | None = None,
    sharing: pd.DataFrame | None = None,
    strict_key: str = "ct_strict",
    relaxed_key: str = "ct_vgene",
) -> pd.DataFrame:
    """Create explicit allowed and not-allowed biological claims."""

    rows = []

    if tcr is not None and strict_key in tcr.columns:
        clone_counts = pd.DataFrame(tcr)[strict_key].dropna().value_counts()
        for clone_id, n_cells in clone_counts.items():
            rows.append(
                {
                    "entity_id": clone_id,
                    "entity_type": "strict_clone",
                    "evidence_level": "confirmed_clone",
                    "n_cells": int(n_cells),
                    "allowed_claim": "Operational paired CDR3-defined TCR clonotype.",
                    "not_allowed_claim": "Antigen-specific or tumor-reactive clone.",
                    "reason": "No antigen validation, TCR reconstruction, pMHC, or functional assay was provided.",
                }
            )

    if risk is not None and not risk.empty:
        risk_df = pd.DataFrame(risk)
        for _, row in risk_df.iterrows():
            group_id = row.get(relaxed_key, row.get("ct_vgene", pd.NA))
            if row.get("risk_label") == "low":
                allowed = "Candidate TCR V-gene group dominated by one strict paired CDR3 clonotype."
                not_allowed = "Confirmed clonotype unless the dominant strict CDR3 is explicitly used."
            else:
                allowed = "Ambiguous relaxed TCR group requiring strict CDR3 confirmation."
                not_allowed = "Single clonotype or confirmed shared clonotype."
            rows.append(
                {
                    "entity_id": group_id,
                    "entity_type": "relaxed_vgene_group",
                    "evidence_level": row.get("evidence_level"),
                    "n_cells": int(row.get("n_cells", 0)),
                    "risk_label": row.get("risk_label"),
                    "dominant_fraction": row.get("dominant_fraction"),
                    "allowed_claim": allowed,
                    "not_allowed_claim": not_allowed,
                    "reason": "TRAV-TRBV can collapse multiple paired CDR3 clonotypes into one relaxed group.",
                }
            )

    if sharing is not None and not sharing.empty:
        sharing_df = pd.DataFrame(sharing)
        for _, row in sharing_df.iterrows():
            pair = f"{row.get('tissue_a')}|{row.get('tissue_b')}|{row.get('analysis_set')}"
            if row.get("vgene_shared_strict_backed", 0) > 0:
                rows.append(
                    {
                        "entity_id": pair,
                        "entity_type": "tissue_sharing",
                        "evidence_level": "strict_sharing",
                        "allowed_claim": "At least one shared relaxed V-gene group is backed by a strict shared CDR3 clonotype.",
                        "not_allowed_claim": "Antigen-specific spatial or functional sharing.",
                        "reason": "Strict paired CDR3 supports clonal sharing, but antigen specificity was not tested.",
                    }
                )
            if row.get("vgene_shared_apparent_only", 0) > 0:
                rows.append(
                    {
                        "entity_id": pair,
                        "entity_type": "tissue_sharing",
                        "evidence_level": "apparent_sharing",
                        "allowed_claim": "Relaxed TRAV-TRBV sharing candidate.",
                        "not_allowed_claim": "Confirmed shared clonotype across tissues.",
                        "reason": "Shared TRAV-TRBV group lacks strict paired CDR3 sharing support.",
                    }
                )

    return pd.DataFrame(rows)
