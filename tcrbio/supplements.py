from __future__ import annotations

from pathlib import Path

import pandas as pd


PER_DATASET_TABLES = {
    "collapse_risk": "vgene_group_table.csv",
    "diversity_compression": "strict_vs_relaxed_diversity.csv",
    "sharing_apparent_vs_strict": "sharing_table.csv",
    "candidate_index": "candidate_table.csv",
    "claim_inventory": "claim_table.csv",
}


def create_supplement_tables(batch_root: str | Path, out_dir: str | Path | None = None) -> dict[str, pd.DataFrame]:
    """Create publication-oriented supplemental tables from a batch output root."""

    root = Path(batch_root)
    out = root / "supplement_tables" if out_dir is None else Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    outputs = {
        "supp_dataset_qc": _dataset_qc(root),
        "supp_collapse_risk": _concat_per_dataset(root, PER_DATASET_TABLES["collapse_risk"]),
        "supp_diversity_compression": _concat_per_dataset(root, PER_DATASET_TABLES["diversity_compression"]),
        "supp_sharing_apparent_vs_strict": _concat_per_dataset(root, PER_DATASET_TABLES["sharing_apparent_vs_strict"]),
        "supp_candidate_index": _concat_per_dataset(root, PER_DATASET_TABLES["candidate_index"]),
        "supp_claim_inventory": _claim_inventory(root),
    }
    for name, table in outputs.items():
        table.to_csv(out / f"{name}.csv", index=False)
    return outputs


def _dataset_qc(root: Path) -> pd.DataFrame:
    summary = _read_optional(root / "batch_run_summary.csv")
    if summary.empty:
        return pd.DataFrame()
    preferred = [
        "result_id",
        "status",
        "n_contexts",
        "n_cells",
        "n_paired_tcr",
        "n_primary_cd4_cd8_paired",
        "primary_paired_fraction",
        "n_strict_clones",
        "n_vgene_groups",
        "n_candidates",
        "mean_richness_ratio",
        "mean_effective_shannon_ratio",
        "strict_shared_clones_total",
        "relaxed_shared_groups_total",
        "apparent_relaxed_only_total",
        "candidate_phenotype_scored_rows",
        "candidate_phenotype_not_evaluable_rows",
    ]
    columns = [column for column in preferred if column in summary.columns]
    remaining = [column for column in summary.columns if column.startswith("risk_groups_")]
    return summary[columns + remaining].copy()


def _claim_inventory(root: Path) -> pd.DataFrame:
    claims = _concat_per_dataset(root, PER_DATASET_TABLES["claim_inventory"])
    if claims.empty or not {"result_id", "entity_type", "evidence_level"}.issubset(claims.columns):
        return claims
    return (
        claims.groupby(["result_id", "entity_type", "evidence_level"], dropna=False)
        .size()
        .rename("n_claims")
        .reset_index()
    )


def _concat_per_dataset(root: Path, filename: str) -> pd.DataFrame:
    frames = []
    per_dataset = root / "per_dataset"
    if not per_dataset.exists():
        return pd.DataFrame()
    for dataset_dir in sorted(path for path in per_dataset.iterdir() if path.is_dir()):
        table = _read_optional(dataset_dir / filename)
        if table.empty:
            continue
        table = table.copy()
        table.insert(0, "result_id", dataset_dir.name)
        frames.append(table)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def _read_optional(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()
