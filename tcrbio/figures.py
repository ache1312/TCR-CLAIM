from __future__ import annotations

from pathlib import Path

import pandas as pd


COLORS = {
    "blue": "#2563eb",
    "green": "#16a34a",
    "amber": "#d97706",
    "red": "#dc2626",
    "gray": "#6b7280",
    "light": "#f3f4f6",
    "ink": "#111827",
}


def create_benchmark_figures(batch_root: str | Path, out_dir: str | Path | None = None) -> dict[str, Path]:
    """Create lightweight SVG figures from a TCR-CLAIM batch output root."""

    root = Path(batch_root)
    out = root / "figures" if out_dir is None else Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    summary = _read_optional(root / "batch_run_summary.csv")
    paths = {
        "richness_ratio_by_dataset": out / "fig_richness_ratio_by_dataset.svg",
        "apparent_sharing_by_dataset": out / "fig_apparent_sharing_by_dataset.svg",
        "candidate_count_by_dataset": out / "fig_candidate_count_by_dataset.svg",
        "collapse_risk_distribution": out / "fig_collapse_risk_distribution.svg",
        "tcr_claim_flow": out / "fig_tcr_claim_flow.svg",
    }
    _bar_svg(
        summary,
        label_col="result_id",
        value_col="mean_richness_ratio",
        output=paths["richness_ratio_by_dataset"],
        title="Relaxed/Strict Richness Ratio",
        y_label="Mean richness ratio",
        color=COLORS["blue"],
    )
    _bar_svg(
        summary,
        label_col="result_id",
        value_col="apparent_relaxed_only_total",
        output=paths["apparent_sharing_by_dataset"],
        title="Apparent Relaxed-Only Sharing",
        y_label="Calls",
        color=COLORS["amber"],
    )
    _bar_svg(
        summary,
        label_col="result_id",
        value_col="n_candidates",
        output=paths["candidate_count_by_dataset"],
        title="Prioritized Candidates",
        y_label="Candidates",
        color=COLORS["green"],
    )
    _stacked_risk_svg(summary, paths["collapse_risk_distribution"])
    _flow_svg(paths["tcr_claim_flow"])
    return paths


def _bar_svg(
    df: pd.DataFrame,
    *,
    label_col: str,
    value_col: str,
    output: Path,
    title: str,
    y_label: str,
    color: str,
) -> None:
    if df.empty or label_col not in df.columns or value_col not in df.columns:
        output.write_text(_empty_svg(title), encoding="utf-8")
        return
    plot = df[df.get("status", "pass") == "pass"].copy() if "status" in df.columns else df.copy()
    plot[value_col] = pd.to_numeric(plot[value_col], errors="coerce")
    plot = plot.dropna(subset=[value_col]).sort_values(value_col, ascending=False)
    width = max(760, 120 + len(plot) * 64)
    height = 460
    margin_left = 70
    margin_bottom = 120
    margin_top = 56
    chart_w = width - margin_left - 32
    chart_h = height - margin_top - margin_bottom
    max_value = float(plot[value_col].max()) if not plot.empty else 0.0
    max_value = max(max_value, 1e-9)
    bar_w = chart_w / max(len(plot), 1) * 0.7
    gap = chart_w / max(len(plot), 1) * 0.3

    parts = [_svg_header(width, height), _title(title, width), _axis(margin_left, margin_top, chart_w, chart_h, y_label)]
    for i, (_, row) in enumerate(plot.iterrows()):
        value = float(row[value_col])
        x = margin_left + i * (bar_w + gap) + gap / 2
        h = chart_h * value / max_value
        y = margin_top + chart_h - h
        label = _escape(str(row[label_col]))
        parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{h:.1f}" fill="{color}"/>')
        parts.append(f'<text x="{x + bar_w / 2:.1f}" y="{y - 6:.1f}" text-anchor="middle" font-size="11" fill="{COLORS["ink"]}">{value:.2f}</text>')
        parts.append(
            f'<text x="{x + bar_w / 2:.1f}" y="{height - margin_bottom + 22}" text-anchor="end" font-size="11" transform="rotate(-45 {x + bar_w / 2:.1f},{height - margin_bottom + 22})">{label}</text>'
        )
    parts.append("</svg>")
    output.write_text("\n".join(parts), encoding="utf-8")


def _stacked_risk_svg(summary: pd.DataFrame, output: Path) -> None:
    risk_cols = [column for column in ["risk_groups_low", "risk_groups_medium", "risk_groups_high"] if column in summary.columns]
    if summary.empty or not risk_cols:
        output.write_text(_empty_svg("Collapse Risk Distribution"), encoding="utf-8")
        return
    plot = summary[summary.get("status", "pass") == "pass"].copy() if "status" in summary.columns else summary.copy()
    width = max(760, 120 + len(plot) * 72)
    height = 460
    margin_left = 70
    margin_bottom = 120
    margin_top = 56
    chart_w = width - margin_left - 32
    chart_h = height - margin_top - margin_bottom
    max_total = max(float(plot[risk_cols].sum(axis=1).max()), 1.0)
    bar_w = chart_w / max(len(plot), 1) * 0.7
    gap = chart_w / max(len(plot), 1) * 0.3
    colors = {"risk_groups_low": COLORS["green"], "risk_groups_medium": COLORS["amber"], "risk_groups_high": COLORS["red"]}
    labels = {"risk_groups_low": "low", "risk_groups_medium": "medium", "risk_groups_high": "high"}

    parts = [_svg_header(width, height), _title("Collapse Risk Distribution", width), _axis(margin_left, margin_top, chart_w, chart_h, "Groups")]
    for i, (_, row) in enumerate(plot.iterrows()):
        x = margin_left + i * (bar_w + gap) + gap / 2
        y_cursor = margin_top + chart_h
        for col in risk_cols:
            value = float(row.get(col, 0) or 0)
            h = chart_h * value / max_total
            y_cursor -= h
            parts.append(f'<rect x="{x:.1f}" y="{y_cursor:.1f}" width="{bar_w:.1f}" height="{h:.1f}" fill="{colors[col]}"/>')
        label = _escape(str(row["result_id"]))
        parts.append(
            f'<text x="{x + bar_w / 2:.1f}" y="{height - margin_bottom + 22}" text-anchor="end" font-size="11" transform="rotate(-45 {x + bar_w / 2:.1f},{height - margin_bottom + 22})">{label}</text>'
        )
    legend_x = width - 220
    for i, col in enumerate(risk_cols):
        y = 24 + i * 18
        parts.append(f'<rect x="{legend_x}" y="{y}" width="12" height="12" fill="{colors[col]}"/>')
        parts.append(f'<text x="{legend_x + 18}" y="{y + 10}" font-size="12">{labels[col]}</text>')
    parts.append("</svg>")
    output.write_text("\n".join(parts), encoding="utf-8")


def _flow_svg(output: Path) -> None:
    width = 980
    height = 240
    steps = [
        ("scRNA/scTCR", "paired CDR3 + metadata"),
        ("Strict Clone", "confirmed operational clone"),
        ("TRAV-TRBV", "candidate relaxed group"),
        ("Risk + Sharing", "collapse and apparent sharing"),
        ("Candidate Card", "allowed and blocked claims"),
        ("Xenium Roadmap", "spatial validation plan"),
    ]
    parts = [_svg_header(width, height), _title("TCR-CLAIM Workflow", width)]
    box_w = 138
    box_h = 78
    y = 92
    for i, (title, subtitle) in enumerate(steps):
        x = 34 + i * 158
        parts.append(f'<rect x="{x}" y="{y}" width="{box_w}" height="{box_h}" rx="6" fill="{COLORS["light"]}" stroke="#d1d5db"/>')
        parts.append(f'<text x="{x + box_w / 2}" y="{y + 32}" text-anchor="middle" font-size="14" font-weight="700">{_escape(title)}</text>')
        parts.append(f'<text x="{x + box_w / 2}" y="{y + 52}" text-anchor="middle" font-size="11" fill="{COLORS["gray"]}">{_escape(subtitle)}</text>')
        if i < len(steps) - 1:
            x1 = x + box_w + 8
            x2 = x + 158 - 10
            parts.append(f'<line x1="{x1}" y1="{y + box_h / 2}" x2="{x2}" y2="{y + box_h / 2}" stroke="{COLORS["gray"]}" stroke-width="2"/>')
            parts.append(f'<polygon points="{x2},{y + box_h / 2} {x2 - 7},{y + box_h / 2 - 5} {x2 - 7},{y + box_h / 2 + 5}" fill="{COLORS["gray"]}"/>')
    parts.append("</svg>")
    output.write_text("\n".join(parts), encoding="utf-8")


def _svg_header(width: int, height: int) -> str:
    return f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">'


def _title(text: str, width: int) -> str:
    return f'<text x="{width / 2}" y="32" text-anchor="middle" font-size="20" font-weight="700" fill="{COLORS["ink"]}">{_escape(text)}</text>'


def _axis(x: int, y: int, w: int, h: int, y_label: str) -> str:
    return "\n".join(
        [
            f'<line x1="{x}" y1="{y + h}" x2="{x + w}" y2="{y + h}" stroke="{COLORS["gray"]}"/>',
            f'<line x1="{x}" y1="{y}" x2="{x}" y2="{y + h}" stroke="{COLORS["gray"]}"/>',
            f'<text x="18" y="{y + h / 2}" text-anchor="middle" font-size="12" transform="rotate(-90 18,{y + h / 2})">{_escape(y_label)}</text>',
        ]
    )


def _empty_svg(title: str) -> str:
    return "\n".join(
        [
            _svg_header(760, 260),
            _title(title, 760),
            '<text x="380" y="140" text-anchor="middle" font-size="14" fill="#6b7280">No plottable rows.</text>',
            "</svg>",
        ]
    )


def _read_optional(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _escape(value: str) -> str:
    return (
        value.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )
