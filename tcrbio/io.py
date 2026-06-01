from __future__ import annotations

from pathlib import Path

import pandas as pd


def read_table(path: str | Path, **kwargs) -> pd.DataFrame:
    """Read a CSV/TSV/TXT table, including gzip-compressed files."""

    path = Path(path)
    lower = path.name.lower()
    if lower.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")):
        return pd.read_csv(path, sep="\t", **kwargs)
    return pd.read_csv(path, **kwargs)


def read_10x_contigs(path: str | Path, **kwargs) -> pd.DataFrame:
    """Read a 10x-style filtered_contig_annotations table."""

    return read_table(path, **kwargs)
