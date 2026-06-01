from __future__ import annotations

from pathlib import Path

import pandas as pd


DEFAULT_SIGNATURES = {
    "t_cell_core": ["CD3D", "CD3E", "TRAC", "TRBC1", "TRBC2", "CD4", "CD8A", "CD8B"],
    "cytotoxic": ["GZMB", "GZMA", "PRF1", "NKG7", "GNLY", "IFNG"],
    "exhausted": ["PDCD1", "TOX", "LAG3", "HAVCR2", "TIGIT", "CTLA4", "CXCL13"],
    "progenitor_memory": ["TCF7", "IL7R", "CCR7", "SELL", "LEF1"],
    "treg": ["FOXP3", "IL2RA", "IKZF2", "TNFRSF18"],
    "antigen_presentation_ifn": ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "CXCL9", "CXCL10", "STAT1"],
}


def design_xenium_panel_from_candidates(
    candidates: pd.DataFrame | None = None,
    *,
    signatures: dict[str, list[str]] | None = None,
    include_top_cdr3: bool = False,
    max_cdr3_targets: int = 20,
) -> pd.DataFrame:
    """Return a marker/CDR3 target table for a Xenium validation roadmap."""

    signatures = DEFAULT_SIGNATURES if signatures is None else signatures
    rows = []
    for source, genes in signatures.items():
        for gene in genes:
            rows.append({"target": gene, "target_type": "gene", "source": source})

    if include_top_cdr3 and candidates is not None:
        df = pd.DataFrame(candidates)
        cdr3_col = "cdr3_sequence" if "cdr3_sequence" in df.columns else None
        if cdr3_col is not None:
            for _, row in df.dropna(subset=[cdr3_col]).head(max_cdr3_targets).iterrows():
                rows.append(
                    {
                        "target": row[cdr3_col],
                        "target_type": "cdr3_custom_candidate",
                        "source": row.get("entity_id", "candidate"),
                    }
                )

    return pd.DataFrame(rows).drop_duplicates().reset_index(drop=True)


def export_cdr3_fasta_for_xenium(
    candidates: pd.DataFrame,
    output: str | Path,
    *,
    sequence_col: str = "cdr3_sequence",
    id_col: str = "entity_id",
    max_targets: int = 20,
) -> Path:
    """Export candidate CDR3 sequences as FASTA for advanced custom design review."""

    df = pd.DataFrame(candidates)
    if sequence_col not in df.columns:
        raise ValueError(f"Missing sequence_col: {sequence_col}")
    output = Path(output)
    lines = []
    for i, (_, row) in enumerate(df.dropna(subset=[sequence_col]).head(max_targets).iterrows(), start=1):
        identifier = str(row.get(id_col, f"candidate_{i}")).replace(" ", "_")
        sequence = str(row[sequence_col]).replace("\n", "").strip()
        lines.extend([f">{identifier}", sequence])
    output.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    return output
