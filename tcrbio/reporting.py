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
