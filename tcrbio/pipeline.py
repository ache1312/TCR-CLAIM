from __future__ import annotations

from pathlib import Path

import pandas as pd

from .claims import claim_checker
from .definitions import define_clones
from .diversity import clonal_diversity, strict_vs_relaxed_diversity
from .io import read_table
from .metrics import clone_count_agreement, collapse_risk, tissue_sharing
from .phenotype import candidate_phenotype_table
from .prioritize import prioritize_candidates
from .qc import qc_summary
from .reporting import report_batch_summary, report_candidate_cards, report_clone_cards, report_markdown_summary
from .tables import cell_tcr_table, strict_clone_table
from .xenium import design_xenium_panel_from_candidates, export_cdr3_fasta_for_xenium


def run_tcr_claim_pipeline(
    input_path: str | Path,
    out_dir: str | Path,
    *,
    tissue_col: str = "tissue_type",
    phenotype_scores: list[str] | None = None,
    include_xenium_cdr3: bool = False,
    max_cdr3_targets: int = 20,
    primary_only: bool = True,
) -> dict[str, pd.DataFrame]:
    """Run one complete TCR-CLAIM table/report workflow."""

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    tcr = read_table(input_path)
    if "ct_strict" not in tcr.columns or "ct_vgene" not in tcr.columns:
        tcr = define_clones(tcr)

    risk = collapse_risk(tcr, primary_only=primary_only)
    sharing = tissue_sharing(tcr, tissue_col=tissue_col, primary_only=primary_only)
    qc = qc_summary(tcr)
    strict_diversity = clonal_diversity(tcr, clone_key="ct_strict", primary_only=primary_only)
    relaxed_diversity = clonal_diversity(tcr, clone_key="ct_vgene", primary_only=primary_only)
    diversity_comparison = strict_vs_relaxed_diversity(tcr, primary_only=primary_only)
    candidates = prioritize_candidates(tcr, risk=risk, primary_only=primary_only)
    candidate_phenotypes = candidate_phenotype_table(
        tcr,
        candidates,
        scores=phenotype_scores,
        primary_only=primary_only,
    )
    claims = claim_checker(tcr=tcr, risk=risk, sharing=sharing, candidates=candidates)
    panel = design_xenium_panel_from_candidates(
        candidates=candidates,
        candidate_phenotypes=candidate_phenotypes,
        include_top_cdr3=include_xenium_cdr3,
        max_cdr3_targets=max_cdr3_targets,
    )

    outputs = {
        "cell_tcr_table": cell_tcr_table(tcr),
        "qc_summary": qc,
        "strict_clone_table": strict_clone_table(tcr, primary_only=primary_only),
        "vgene_group_table": risk,
        "clone_count_agreement": clone_count_agreement(tcr, primary_only=primary_only),
        "strict_clonal_diversity": strict_diversity,
        "relaxed_clonal_diversity": relaxed_diversity,
        "strict_vs_relaxed_diversity": diversity_comparison,
        "sharing_table": sharing,
        "candidate_table": candidates,
        "candidate_phenotype_table": candidate_phenotypes,
        "claim_table": claims,
        "xenium_panel_roadmap": panel,
    }
    for name, table in outputs.items():
        table.to_csv(out / f"{name}.csv", index=False)

    if include_xenium_cdr3:
        export_cdr3_fasta_for_xenium(candidates, out / "xenium_cdr3_targets.fasta", max_targets=max_cdr3_targets)
    report_clone_cards(claims, output=out / "clone_cards.html")
    report_candidate_cards(
        candidates=candidates,
        candidate_phenotypes=candidate_phenotypes,
        output=out / "candidate_cards.html",
    )
    report_markdown_summary(
        output=out / "tcr_claim_report.md",
        qc=qc,
        diversity=diversity_comparison,
        candidates=candidates,
        candidate_phenotypes=candidate_phenotypes,
        sharing=sharing,
        claims=claims,
        risk=risk,
    )
    return outputs


def discover_result_dirs(results_root: str | Path) -> list[Path]:
    """Find benchmark result directories containing cell_metadata_with_tcr.csv."""

    root = Path(results_root)
    if (root / "cell_metadata_with_tcr.csv").exists():
        return [root]
    return sorted({path.parent for path in root.rglob("cell_metadata_with_tcr.csv")})


def run_tcr_claim_batch(
    result_dirs: list[str | Path],
    out_root: str | Path,
    *,
    tissue_col: str = "tissue_type",
    phenotype_scores: list[str] | None = None,
    include_xenium_cdr3: bool = False,
    max_cdr3_targets: int = 20,
    primary_only: bool = True,
    continue_on_error: bool = True,
) -> pd.DataFrame:
    """Run TCR-CLAIM across multiple benchmark result directories."""

    out = Path(out_root)
    per_dataset = out / "per_dataset"
    per_dataset.mkdir(parents=True, exist_ok=True)
    rows = []
    for result_dir in result_dirs:
        result_path = Path(result_dir)
        result_id = result_path.name
        input_path = result_path / "cell_metadata_with_tcr.csv"
        output_dir = per_dataset / _safe_path_name(result_id)
        if not input_path.exists():
            rows.append(_failed_summary(result_id, input_path, output_dir, "missing cell_metadata_with_tcr.csv"))
            if not continue_on_error:
                break
            continue
        try:
            run_tcr_claim_pipeline(
                input_path,
                output_dir,
                tissue_col=tissue_col,
                phenotype_scores=phenotype_scores,
                include_xenium_cdr3=include_xenium_cdr3,
                max_cdr3_targets=max_cdr3_targets,
                primary_only=primary_only,
            )
            rows.append(summarize_tcr_claim_outputs(output_dir, result_id=result_id, input_path=input_path))
        except Exception as exc:
            rows.append(_failed_summary(result_id, input_path, output_dir, str(exc)))
            if not continue_on_error:
                raise

    summary = pd.DataFrame(rows)
    summary.to_csv(out / "batch_run_summary.csv", index=False)
    report_batch_summary(summary, output=out / "batch_report.md")
    return summary


def summarize_tcr_claim_outputs(
    out_dir: str | Path,
    *,
    result_id: str | None = None,
    input_path: str | Path | None = None,
) -> dict[str, object]:
    """Summarize one TCR-CLAIM output directory into a publication-oriented row."""

    out = Path(out_dir)
    qc = _read_optional(out / "qc_summary.csv")
    strict = _read_optional(out / "strict_clone_table.csv")
    risk = _read_optional(out / "vgene_group_table.csv")
    diversity = _read_optional(out / "strict_vs_relaxed_diversity.csv")
    sharing = _read_optional(out / "sharing_table.csv")
    candidates = _read_optional(out / "candidate_table.csv")
    phenotypes = _read_optional(out / "candidate_phenotype_table.csv")
    claims = _read_optional(out / "claim_table.csv")

    top = candidates.sort_values("candidate_rank").head(1) if "candidate_rank" in candidates.columns else candidates.head(1)
    row = {
        "result_id": result_id if result_id is not None else out.name,
        "status": "pass",
        "input_path": str(input_path) if input_path is not None else pd.NA,
        "output_dir": str(out),
        "n_contexts": int(len(qc)) if not qc.empty else 0,
        "n_cells": int(_sum(qc, "n_cells")),
        "n_paired_tcr": int(_sum(qc, "n_paired_tcr")),
        "n_primary_cd4_cd8_paired": int(_sum(qc, "n_primary_cd4_cd8_paired")),
        "primary_paired_fraction": _ratio(_sum(qc, "n_primary_cd4_cd8_paired"), _sum(qc, "n_cells")),
        "n_strict_clones": int(strict["ct_strict"].nunique()) if "ct_strict" in strict.columns else 0,
        "n_vgene_groups": int(risk["ct_vgene"].nunique()) if "ct_vgene" in risk.columns else 0,
        "n_candidates": int(len(candidates)),
        "top_candidate_n_cells": int(top["n_cells"].iloc[0]) if not top.empty and "n_cells" in top.columns else 0,
        "top_candidate_risk_label": top["risk_label"].iloc[0] if not top.empty and "risk_label" in top.columns else pd.NA,
        "mean_richness_ratio": _mean(diversity, "richness_ratio"),
        "mean_effective_shannon_ratio": _mean(diversity, "effective_shannon_ratio"),
        "strict_shared_clones_total": int(_sum(sharing, "strict_shared_clones")),
        "relaxed_shared_groups_total": int(_sum(sharing, "vgene_shared_groups")),
        "apparent_relaxed_only_total": int(_sum(sharing, "vgene_shared_apparent_only")),
        "candidate_phenotype_scored_rows": int((phenotypes.get("phenotype_evidence_status", pd.Series(dtype=object)) == "scored").sum()),
        "candidate_phenotype_not_evaluable_rows": int((phenotypes.get("phenotype_evidence_status", pd.Series(dtype=object)) == "no_scores_available").sum()),
        "n_claims": int(len(claims)),
    }
    if "risk_label" in risk.columns:
        for label, value in risk["risk_label"].value_counts(dropna=False).items():
            row[f"risk_groups_{label}"] = int(value)
    return row


def _read_optional(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _sum(df: pd.DataFrame, column: str) -> float:
    if column not in df.columns:
        return 0.0
    return float(pd.to_numeric(df[column], errors="coerce").fillna(0).sum())


def _mean(df: pd.DataFrame, column: str) -> float:
    if column not in df.columns:
        return float("nan")
    return float(pd.to_numeric(df[column], errors="coerce").mean())


def _ratio(numerator: float, denominator: float) -> float:
    if denominator == 0:
        return float("nan")
    return numerator / denominator


def _failed_summary(result_id: str, input_path: Path, output_dir: Path, reason: str) -> dict[str, object]:
    return {
        "result_id": result_id,
        "status": "fail",
        "input_path": str(input_path),
        "output_dir": str(output_dir),
        "failure_reason": reason,
    }


def _safe_path_name(value: str) -> str:
    return "".join(char if char.isalnum() or char in {"-", "_", "."} else "_" for char in value)
