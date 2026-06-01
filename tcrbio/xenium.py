from __future__ import annotations

from pathlib import Path

import pandas as pd


DEFAULT_SIGNATURES = {
    "t_cell_core": ["CD3D", "CD3E", "TRAC", "TRBC1", "TRBC2", "CD4", "CD8A", "CD8B"],
    "cytotoxic": ["GZMB", "GZMA", "PRF1", "NKG7", "GNLY", "IFNG"],
    "exhausted": ["PDCD1", "TOX", "LAG3", "HAVCR2", "TIGIT", "CTLA4", "CXCL13"],
    "proliferative": ["MKI67", "TOP2A", "STMN1"],
    "progenitor_memory": ["TCF7", "IL7R", "CCR7", "SELL", "LEF1"],
    "treg": ["FOXP3", "IL2RA", "IKZF2", "TNFRSF18"],
    "antigen_presentation_ifn": ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "CXCL9", "CXCL10", "STAT1"],
    "tumor_epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "myeloid": ["LYZ", "C1QA", "C1QB", "CD68", "CD163", "FCGR3A"],
    "stromal": ["COL1A1", "COL1A2", "ACTA2", "TAGLN"],
    "endothelial": ["PECAM1", "VWF"],
}

STATE_TO_SIGNATURE = {
    "cytotoxic": "cytotoxic",
    "exhausted": "exhausted",
    "proliferative": "proliferative",
    "progenitor_memory": "progenitor_memory",
    "treg": "treg",
    "ifn_response": "antigen_presentation_ifn",
}


def design_xenium_panel_from_candidates(
    candidates: pd.DataFrame | None = None,
    *,
    candidate_phenotypes: pd.DataFrame | None = None,
    signatures: dict[str, list[str]] | None = None,
    include_top_cdr3: bool = False,
    max_cdr3_targets: int = 20,
    top_n_candidates: int = 20,
) -> pd.DataFrame:
    """Return a marker/CDR3 target table for a Xenium validation roadmap."""

    signatures = DEFAULT_SIGNATURES if signatures is None else signatures
    rows = []
    for source, genes in signatures.items():
        for gene in genes:
            priority = "core" if source == "t_cell_core" else "roadmap"
            rows.append(
                {
                    "target": gene,
                    "target_type": "gene",
                    "source": source,
                    "priority": priority,
                    "rationale": "Default TCR-CLAIM spatial validation marker.",
                }
            )

    phenotype_df = pd.DataFrame() if candidate_phenotypes is None else pd.DataFrame(candidate_phenotypes)
    if not phenotype_df.empty and {"phenotype_state", "direction"}.issubset(phenotype_df.columns):
        supported = phenotype_df[
            (phenotype_df["direction"] == "enriched")
            & phenotype_df["phenotype_state"].notna()
        ].copy()
        if "candidate_rank" in supported.columns:
            supported = supported.sort_values("candidate_rank")
        states = list(dict.fromkeys(supported["phenotype_state"].astype(str)))
        for state in states:
            signature = STATE_TO_SIGNATURE.get(state, state)
            for gene in signatures.get(signature, []):
                rows.append(
                    {
                        "target": gene,
                        "target_type": "gene",
                        "source": f"candidate_state:{state}",
                        "priority": "candidate_supported",
                        "rationale": "Included because prioritized TCR candidates are enriched for this phenotype score.",
                    }
                )

    if include_top_cdr3 and candidates is not None:
        df = pd.DataFrame(candidates)
        if "candidate_rank" in df.columns:
            df = df.sort_values("candidate_rank")
        cdr3_cols = [column for column in ["cdr3_sequence", "trb_cdr3_nt", "tra_cdr3_nt", "trb_cdr3_aa", "tra_cdr3_aa"] if column in df.columns]
        emitted = 0
        for cdr3_col in cdr3_cols:
            for _, row in df.dropna(subset=[cdr3_col]).head(top_n_candidates).iterrows():
                if emitted >= max_cdr3_targets:
                    break
                rows.append(
                    {
                        "target": row[cdr3_col],
                        "target_type": "cdr3_custom_candidate",
                        "source": row.get("candidate_id", row.get("entity_id", "candidate")),
                        "priority": "experimental",
                        "rationale": "Optional advanced custom CDR3 target requiring assay-specific review.",
                    }
                )
                emitted += 1
            if emitted >= max_cdr3_targets:
                break

    return _deduplicate_panel_rows(pd.DataFrame(rows))


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
        fallback_cols = ["trb_cdr3_nt", "tra_cdr3_nt", "trb_cdr3_aa", "tra_cdr3_aa"] if sequence_col == "cdr3_sequence" else []
        sequence_col = next((column for column in fallback_cols if column in df.columns), sequence_col)
    if sequence_col not in df.columns:
        raise ValueError(f"Missing sequence_col: {sequence_col}")
    output = Path(output)
    lines = []
    for i, (_, row) in enumerate(df.dropna(subset=[sequence_col]).head(max_targets).iterrows(), start=1):
        identifier = str(row.get(id_col, row.get("candidate_id", f"candidate_{i}"))).replace(" ", "_")
        sequence = str(row[sequence_col]).replace("\n", "").strip()
        lines.extend([f">{identifier}", sequence])
    output.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    return output


def _deduplicate_panel_rows(panel: pd.DataFrame) -> pd.DataFrame:
    if panel.empty:
        return panel
    priority_order = {"core": 0, "candidate_supported": 1, "roadmap": 2, "experimental": 3}
    panel = panel.copy()
    panel["_priority_order"] = panel["priority"].map(priority_order).fillna(99)
    panel = panel.sort_values(["target_type", "target", "_priority_order"])
    rows = []
    for (target_type, target), group in panel.groupby(["target_type", "target"], dropna=False, sort=False):
        rows.append(
            {
                "target": target,
                "target_type": target_type,
                "source": ";".join(dict.fromkeys(group["source"].astype(str))),
                "priority": group.sort_values("_priority_order")["priority"].iloc[0],
                "rationale": "; ".join(dict.fromkeys(group["rationale"].astype(str))),
            }
        )
    return pd.DataFrame(rows).reset_index(drop=True)
