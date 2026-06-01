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


def report_candidate_cards(
    *,
    candidates: pd.DataFrame,
    candidate_phenotypes: pd.DataFrame | None = None,
    output: str | Path | None = None,
    top_n: int = 25,
) -> str:
    """Render richer candidate-centered HTML cards for biological review."""

    candidate_df = pd.DataFrame(candidates).copy()
    phenotype_df = pd.DataFrame() if candidate_phenotypes is None else pd.DataFrame(candidate_phenotypes)
    if "candidate_rank" in candidate_df.columns:
        candidate_df = candidate_df.sort_values("candidate_rank")
    cards = []
    for _, candidate in candidate_df.head(top_n).iterrows():
        candidate_id = candidate.get("candidate_id", candidate.get("ct_strict", "candidate"))
        phenotype_rows = _phenotype_rows_for_candidate(phenotype_df, candidate_id)
        cards.append(_candidate_card_html(candidate, phenotype_rows))

    html = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>TCR-CLAIM candidate cards</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 32px; color: #202124; line-height: 1.45; }}
    .card {{ border: 1px solid #d0d7de; border-radius: 6px; padding: 16px; margin: 16px 0; }}
    .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 8px 18px; }}
    .label {{ color: #57606a; font-size: 12px; text-transform: uppercase; }}
    .value {{ font-weight: 600; overflow-wrap: anywhere; }}
    h1 {{ margin-bottom: 8px; }}
    h2 {{ margin-top: 0; font-size: 18px; overflow-wrap: anywhere; }}
    table {{ border-collapse: collapse; width: 100%; margin-top: 8px; }}
    th, td {{ border: 1px solid #d0d7de; padding: 6px; text-align: left; font-size: 13px; }}
    th {{ background: #f6f8fa; }}
  </style>
</head>
<body>
  <h1>TCR-CLAIM candidate cards</h1>
  <p>These cards summarize candidate evidence for biological review. They do not claim antigen specificity, tumor reactivity, functional killing, or spatial validation.</p>
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


def report_batch_summary(summary: pd.DataFrame, output: str | Path | None = None, top_n: int = 20) -> str:
    """Render a compact Markdown summary for a multi-dataset TCR-CLAIM run."""

    df = pd.DataFrame(summary).copy()
    sections = [
        "# TCR-CLAIM Batch Summary",
        "",
        "This report summarizes multi-dataset TCR-CLAIM outputs. It is intended for benchmark tracking, supplemental tables, and manuscript triage.",
        "",
        "## Interpretation Boundary",
        "",
        "- Batch metrics compare candidate evidence across datasets.",
        "- They do not establish antigen specificity, tumor reactivity, functional killing, or spatial validation.",
    ]
    if df.empty:
        sections.extend(["", "No datasets were processed."])
    else:
        passed = int((df.get("status", pd.Series(dtype=object)) == "pass").sum())
        failed = int((df.get("status", pd.Series(dtype=object)) == "fail").sum())
        sections.extend(
            [
                "",
                "## Run Status",
                "",
                f"- Datasets: {len(df)}",
                f"- Passed: {passed}",
                f"- Failed: {failed}",
                f"- Cells: {_fmt_int(_sum_numeric(df, 'n_cells'))}",
                f"- Primary CD4/CD8 paired cells: {_fmt_int(_sum_numeric(df, 'n_primary_cd4_cd8_paired'))}",
                f"- Candidates: {_fmt_int(_sum_numeric(df, 'n_candidates'))}",
                f"- Apparent relaxed-only sharing calls: {_fmt_int(_sum_numeric(df, 'apparent_relaxed_only_total'))}",
            ]
        )
        table_cols = [
            column
            for column in [
                "result_id",
                "status",
                "n_cells",
                "n_primary_cd4_cd8_paired",
                "n_candidates",
                "mean_richness_ratio",
                "apparent_relaxed_only_total",
                "candidate_phenotype_scored_rows",
            ]
            if column in df.columns
        ]
        sections.extend(["", "## Dataset Index", "", _markdown_table(df[table_cols].head(top_n))])
        failed_df = df[df.get("status", pd.Series(index=df.index, dtype=object)) == "fail"]
        if not failed_df.empty:
            fail_cols = [column for column in ["result_id", "failure_reason", "input_path"] if column in failed_df.columns]
            sections.extend(["", "## Failures", "", _markdown_table(failed_df[fail_cols])])

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


def _candidate_card_html(candidate: pd.Series, phenotype_rows: pd.DataFrame) -> str:
    candidate_id = candidate.get("candidate_id", candidate.get("ct_strict", "candidate"))
    rank = candidate.get("candidate_rank", "")
    context = " / ".join(
        str(candidate.get(column))
        for column in ["dataset_id", "donor_id", "sample_id", "tissue_type", "cell_class"]
        if column in candidate.index and pd.notna(candidate.get(column))
    )
    phenotype_html = _candidate_phenotype_html(phenotype_rows)
    fields = [
        ("Rank", rank),
        ("Context", context),
        ("Cells", candidate.get("n_cells", "")),
        ("TRAV-TRBV", candidate.get("ct_vgene", candidate.get("relaxed_group", ""))),
        ("Collapse Risk", candidate.get("risk_label", "")),
        ("Dominant Fraction", _fmt_float(candidate.get("dominant_fraction"))),
        ("Rank Score", _fmt_float(candidate.get("rank_score"))),
        ("Evidence", candidate.get("evidence_level", "confirmed_clone")),
    ]
    field_html = "\n".join(
        f"""<div><div class="label">{escape(str(label))}</div><div class="value">{escape(str(value))}</div></div>"""
        for label, value in fields
    )
    allowed = escape(str(candidate.get("allowed_claim", "Prioritized strict paired-CDR3 TCR clone candidate.")))
    blocked = escape(str(candidate.get("not_allowed_claim", "Antigen-specific or tumor-reactive clone without orthogonal validation.")))
    return f"""
<section class="card">
  <h2>#{escape(str(rank))} {escape(str(candidate_id))}</h2>
  <div class="grid">{field_html}</div>
  <h3>Phenotype Evidence</h3>
  {phenotype_html}
  <h3>Claim Boundary</h3>
  <p><strong>Allowed:</strong> {allowed}</p>
  <p><strong>Not allowed:</strong> {blocked}</p>
</section>"""


def _phenotype_rows_for_candidate(phenotype_df: pd.DataFrame, candidate_id: object) -> pd.DataFrame:
    if phenotype_df.empty or "candidate_id" not in phenotype_df.columns:
        return pd.DataFrame()
    return phenotype_df[phenotype_df["candidate_id"] == candidate_id].copy()


def _candidate_phenotype_html(rows: pd.DataFrame) -> str:
    if rows.empty:
        return "<p>No candidate phenotype table was provided.</p>"
    status_values = set(rows.get("phenotype_evidence_status", pd.Series(dtype=object)).dropna().astype(str))
    if status_values == {"no_scores_available"}:
        return "<p>No phenotype score was available for this candidate. TCR evidence can prioritize the clone, but phenotype association is not evaluable.</p>"
    columns = [
        column
        for column in ["phenotype_state", "score_column", "direction", "candidate_mean", "reference_mean", "delta_mean"]
        if column in rows.columns
    ]
    return _html_table(rows[columns].head(8))


def _html_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "<p>No rows.</p>"
    header = "".join(f"<th>{escape(str(column))}</th>" for column in df.columns)
    rows = []
    for _, row in df.iterrows():
        cells = "".join(f"<td>{escape(_fmt_cell(value))}</td>" for value in row)
        rows.append(f"<tr>{cells}</tr>")
    return f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(rows)}</tbody></table>"


def _fmt_cell(value: object) -> str:
    if pd.isna(value):
        return ""
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value)


def _fmt_float(value: object) -> str:
    if _is_missing(value):
        return ""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if pd.isna(number):
        return ""
    return f"{number:.3f}"


def _is_missing(value: object) -> bool:
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False
