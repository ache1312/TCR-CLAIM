from __future__ import annotations

from pathlib import Path
from html import escape

import pandas as pd


def report_clone_cards(claims: pd.DataFrame, output: str | Path | None = None) -> str:
    """Render claim rows as a simple HTML clone-card report."""

    df = pd.DataFrame(claims).copy()
    cards = []
    for _, row in df.iterrows():
        title = escape(str(row.get("entity_id", "candidate")))
        evidence = escape(str(row.get("evidence_level", "not_evaluable")))
        allowed = escape(str(row.get("allowed_claim", "")))
        not_allowed = escape(str(row.get("not_allowed_claim", "")))
        reason = escape(str(row.get("reason", "")))
        cards.append(
            f"""
<section class="card">
  <h2>{title}</h2>
  <p><strong>Evidence level:</strong> {evidence}</p>
  <p><strong>Allowed claim:</strong> {allowed}</p>
  <p><strong>Not allowed:</strong> {not_allowed}</p>
  <p><strong>Reason:</strong> {reason}</p>
</section>"""
        )

    html = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>TCR candidate clone cards</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 32px; color: #202124; }}
    .card {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 16px; margin: 16px 0; }}
    h1 {{ margin-bottom: 8px; }}
    h2 {{ margin-top: 0; font-size: 18px; }}
  </style>
</head>
<body>
  <h1>TCR candidate clone cards</h1>
  <p>These cards separate supported biological claims from unsupported claims.</p>
  {''.join(cards)}
</body>
</html>"""

    if output is not None:
        Path(output).write_text(html, encoding="utf-8")
    return html


def report_markdown_summary(
    *,
    output: str | Path | None = None,
    qc: pd.DataFrame | None = None,
    diversity: pd.DataFrame | None = None,
    candidates: pd.DataFrame | None = None,
    candidate_phenotypes: pd.DataFrame | None = None,
    sharing: pd.DataFrame | None = None,
    claims: pd.DataFrame | None = None,
    risk: pd.DataFrame | None = None,
    top_n: int = 10,
) -> str:
    """Render a compact Markdown summary for one TCR-CLAIM run."""

    sections = [
        "# TCR-CLAIM Run Summary",
        "",
        "TCR-CLAIM separates strict paired-CDR3 clone evidence from relaxed TRAV-TRBV candidate groups, collapse risk, apparent sharing, and unsupported biological claims.",
        "",
        "## Interpretation Boundary",
        "",
        "- Supported: operational paired-CDR3 clonotypes, relaxed candidate groups, collapse risk, strict/apparent sharing, and candidates for downstream biological review.",
        "- Not supported by TCR-CLAIM alone: antigen specificity, tumor reactivity, functional killing, or spatial validation.",
    ]

    qc_df = pd.DataFrame() if qc is None else pd.DataFrame(qc)
    if not qc_df.empty:
        total_cells = _sum_numeric(qc_df, "n_cells")
        paired = _sum_numeric(qc_df, "n_paired_tcr")
        primary = _sum_numeric(qc_df, "n_primary_cd4_cd8_paired")
        sections.extend(
            [
                "",
                "## QC",
                "",
                f"- Contexts: {len(qc_df)}",
                f"- Cells: {_fmt_int(total_cells)}",
                f"- Paired TCR cells: {_fmt_int(paired)} ({_fmt_fraction(paired, total_cells)})",
                f"- Primary CD4/CD8 paired cells: {_fmt_int(primary)} ({_fmt_fraction(primary, total_cells)})",
            ]
        )

    diversity_df = pd.DataFrame() if diversity is None else pd.DataFrame(diversity)
    if not diversity_df.empty:
        sections.extend(["", "## Strict Versus Relaxed Diversity", ""])
        for column, label in [
            ("richness_ratio", "Mean relaxed/strict richness ratio"),
            ("richness_relative_difference", "Mean relative richness difference"),
            ("effective_shannon_ratio", "Mean relaxed/strict effective Shannon ratio"),
        ]:
            if column in diversity_df.columns:
                value = pd.to_numeric(diversity_df[column], errors="coerce").dropna()
                if not value.empty:
                    sections.append(f"- {label}: {value.mean():.3f}")

    risk_df = pd.DataFrame() if risk is None else pd.DataFrame(risk)
    if not risk_df.empty and "risk_label" in risk_df.columns:
        risk_counts = risk_df["risk_label"].value_counts(dropna=False).rename_axis("risk_label").reset_index(name="n_groups")
        sections.extend(["", "## Collapse Risk", "", _markdown_table(risk_counts)])

    sharing_df = pd.DataFrame() if sharing is None else pd.DataFrame(sharing)
    if not sharing_df.empty:
        sections.extend(
            [
                "",
                "## Tissue Sharing",
                "",
                f"- Strict shared clones across reported tissue pairs: {_fmt_int(_sum_numeric(sharing_df, 'strict_shared_clones'))}",
                f"- Relaxed shared TRAV-TRBV groups: {_fmt_int(_sum_numeric(sharing_df, 'vgene_shared_groups'))}",
                f"- Apparent relaxed-only shared groups: {_fmt_int(_sum_numeric(sharing_df, 'vgene_shared_apparent_only'))}",
            ]
        )

    candidate_df = pd.DataFrame() if candidates is None else pd.DataFrame(candidates)
    if not candidate_df.empty:
        columns = [
            column
            for column in ["candidate_rank", "candidate_id", "n_cells", "tissue_type", "cell_class", "risk_label", "rank_score"]
            if column in candidate_df.columns
        ]
        top = candidate_df[columns].head(top_n).copy()
        if "candidate_id" in top.columns:
            top["candidate_id"] = top["candidate_id"].map(_shorten)
        sections.extend(["", "## Prioritized Candidates", "", _markdown_table(top)])

    phenotype_df = pd.DataFrame() if candidate_phenotypes is None else pd.DataFrame(candidate_phenotypes)
    if not phenotype_df.empty:
        sections.extend(["", "## Candidate Phenotypes", ""])
        status_counts = phenotype_df["phenotype_evidence_status"].value_counts(dropna=False).rename_axis("status").reset_index(name="n_rows")
        sections.append(_markdown_table(status_counts))
        enriched = phenotype_df[phenotype_df.get("direction", pd.Series(index=phenotype_df.index, dtype=object)) == "enriched"]
        if not enriched.empty:
            columns = [
                column
                for column in ["candidate_rank", "candidate_id", "phenotype_state", "score_column", "delta_mean"]
                if column in enriched.columns
            ]
            top_enriched = enriched.sort_values("delta_mean", ascending=False)[columns].head(top_n).copy()
            if "candidate_id" in top_enriched.columns:
                top_enriched["candidate_id"] = top_enriched["candidate_id"].map(_shorten)
            sections.extend(["", "Top enriched phenotype shifts:", "", _markdown_table(top_enriched)])

    claim_df = pd.DataFrame() if claims is None else pd.DataFrame(claims)
    if not claim_df.empty and {"entity_type", "evidence_level"}.issubset(claim_df.columns):
        claim_counts = (
            claim_df.groupby(["entity_type", "evidence_level"], dropna=False)
            .size()
            .rename("n_claims")
            .reset_index()
        )
        sections.extend(["", "## Claim Inventory", "", _markdown_table(claim_counts)])

    markdown = "\n".join(sections).rstrip() + "\n"
    if output is not None:
        Path(output).write_text(markdown, encoding="utf-8")
    return markdown


def _sum_numeric(df: pd.DataFrame, column: str) -> float:
    if column not in df.columns:
        return 0.0
    return float(pd.to_numeric(df[column], errors="coerce").fillna(0).sum())


def _fmt_int(value: float) -> str:
    return f"{int(value):,}"


def _fmt_fraction(numerator: float, denominator: float) -> str:
    if denominator == 0:
        return "not evaluable"
    return f"{numerator / denominator:.3f}"


def _shorten(value: object, max_len: int = 48) -> str:
    text = str(value)
    if len(text) <= max_len:
        return text
    return f"{text[: max_len - 3]}..."


def _markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    out = pd.DataFrame(df).copy()
    for column in out.columns:
        if pd.api.types.is_float_dtype(out[column]):
            out[column] = out[column].map(lambda value: "" if pd.isna(value) else f"{value:.3f}")
        else:
            out[column] = out[column].map(lambda value: "" if pd.isna(value) else str(value))
    header = "| " + " | ".join(out.columns) + " |"
    separator = "| " + " | ".join(["---"] * len(out.columns)) + " |"
    rows = ["| " + " | ".join(row) + " |" for row in out.astype(str).to_numpy()]
    return "\n".join([header, separator, *rows])
