from __future__ import annotations

import pandas as pd


def primary_tcr_cells(
    tcr: pd.DataFrame,
    *,
    strict_key: str = "ct_strict",
    relaxed_key: str = "ct_vgene",
    cell_class_col: str = "cell_class",
    allowed_cell_classes: tuple[str, ...] = ("CD4", "CD8"),
) -> pd.DataFrame:
    """Return the benchmark primary universe: CD4/CD8 cells with paired TCR keys.

    If `cell_class_col` is absent, the function keeps all cells with both TCR
    keys. This keeps table-only use cases possible while preserving the
    CD4/CD8 filter whenever the metadata are available.
    """

    df = pd.DataFrame(tcr).copy()
    keep = df[strict_key].notna() & df[relaxed_key].notna()
    if cell_class_col in df.columns:
        keep &= df[cell_class_col].isin(allowed_cell_classes)
    return df[keep].copy()
