from __future__ import annotations

import pandas as pd


def _safe_component(series: pd.Series) -> pd.Series:
    return series.astype(object).where(series.notna(), "None").astype(str)


def define_clones(
    tcr: pd.DataFrame,
    *,
    strict_policy: str = "paired_cdr3_nt_vjc",
    relaxed_policy: str = "trav_trbv",
) -> pd.DataFrame:
    """Add strict paired CDR3 and relaxed TRAV-TRBV clone keys."""

    if strict_policy not in {"paired_cdr3_nt_vjc", "paired_nt_vjc"}:
        raise ValueError(f"Unsupported strict_policy: {strict_policy}")
    if relaxed_policy != "trav_trbv":
        raise ValueError(f"Unsupported relaxed_policy: {relaxed_policy}")

    out = pd.DataFrame(tcr).copy()
    needed = [
        "trav",
        "traj",
        "trac",
        "tra_cdr3_nt",
        "trbv",
        "trbd",
        "trbj",
        "trbc",
        "trb_cdr3_nt",
    ]
    for column in needed:
        if column not in out.columns:
            out[column] = pd.NA

    strict_ready = out[["trav", "trbv", "tra_cdr3_nt", "trb_cdr3_nt"]].notna().all(axis=1)
    relaxed_ready = out[["trav", "trbv"]].notna().all(axis=1)

    strict = (
        _safe_component(out["trav"])
        + "."
        + _safe_component(out["traj"])
        + "."
        + _safe_component(out["trac"])
        + ";"
        + _safe_component(out["tra_cdr3_nt"])
        + "_"
        + _safe_component(out["trbv"])
        + "."
        + _safe_component(out["trbd"])
        + "."
        + _safe_component(out["trbj"])
        + "."
        + _safe_component(out["trbc"])
        + ";"
        + _safe_component(out["trb_cdr3_nt"])
    )
    relaxed = _safe_component(out["trav"]) + "_" + _safe_component(out["trbv"])

    out["ct_strict"] = strict.where(strict_ready, pd.NA)
    out["ct_vgene"] = relaxed.where(relaxed_ready, pd.NA)
    out["has_paired_tcr"] = strict_ready
    out["ct_strict_policy"] = strict_policy
    out["ct_strict_level"] = "confirmed_clone"
    out["ct_vgene_policy"] = relaxed_policy
    out["ct_vgene_level"] = "candidate_group"
    return out
