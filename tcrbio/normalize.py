from __future__ import annotations

from collections.abc import Mapping

import pandas as pd


MISSING_TEXT = {"", "NA", "NaN", "nan", "NULL", "None", "none"}


def _first_present(df: pd.DataFrame, candidates: list[str | None]) -> str | None:
    for candidate in candidates:
        if candidate and candidate in df.columns:
            return candidate
    return None


def _pull(df: pd.DataFrame, candidates: list[str | None], default=None) -> pd.Series:
    column = _first_present(df, candidates)
    if column is None:
        return pd.Series([default] * len(df), index=df.index)
    return df[column]


def _as_text(series: pd.Series) -> pd.Series:
    out = series.astype("string").str.strip()
    out = out.mask(out.isin(MISSING_TEXT))
    return out.astype(object).where(pd.notna(out), None)


def _as_number(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series
    if pd.api.types.is_numeric_dtype(series):
        return series.notna() & (series != 0)
    text = series.astype("string").str.lower().str.strip()
    true_values = {"true", "t", "yes", "y", "1", "productive"}
    false_values = {"false", "f", "no", "n", "0", "unproductive"}
    return text.map(lambda x: True if x in true_values else False if x in false_values else pd.NA)


def _normalize_chain(series: pd.Series) -> pd.Series:
    text = series.astype("string").str.upper().str.strip().str.replace(r"\s+", "", regex=True)
    mapping = {
        "TRA": "TRA",
        "TCRA": "TRA",
        "ALPHA": "TRA",
        "A": "TRA",
        "TRB": "TRB",
        "TCRB": "TRB",
        "BETA": "TRB",
        "B": "TRB",
    }
    return text.map(mapping).astype(object).where(lambda x: pd.notna(x), None)


def _apply_metadata(
    wide: pd.DataFrame,
    sample_metadata: Mapping[str, object] | pd.DataFrame | None,
) -> pd.DataFrame:
    if sample_metadata is None:
        return wide

    if isinstance(sample_metadata, Mapping):
        out = wide.copy()
        for key, value in sample_metadata.items():
            out[key] = value
        return out

    meta = pd.DataFrame(sample_metadata).copy()
    if len(meta) == 1:
        out = wide.copy()
        for column in meta.columns:
            out[column] = meta.iloc[0][column]
        return out

    barcode_column = _first_present(meta, ["cell_barcode", "barcode", "cell", "cell_id"])
    if barcode_column is None:
        raise ValueError("sample_metadata with multiple rows must include a barcode column")
    meta = meta.rename(columns={barcode_column: "cell_barcode"})
    return wide.merge(meta, on="cell_barcode", how="left")


def normalize_contigs(
    contigs: pd.DataFrame,
    *,
    barcode_col: str | None = None,
    sample_metadata: Mapping[str, object] | pd.DataFrame | None = None,
    dataset_id: str | None = None,
    cancer_type: str | None = None,
    donor_id: str | None = None,
    sample_id: str | None = None,
    tissue_type: str | None = None,
) -> pd.DataFrame:
    """Normalize long TCR contigs to one dominant TRA and TRB per cell.

    This function intentionally does not downgrade paired TRA/TRB requirements.
    Cells lacking either chain remain in the table but do not become strict
    paired clonotypes until `define_clones` is called.
    """

    df = pd.DataFrame(contigs).copy()
    std = pd.DataFrame(
        {
            "cell_barcode": _as_text(_pull(df, [barcode_col, "cell_barcode", "barcode", "cell", "cell_id"])),
            "chain": _normalize_chain(_pull(df, ["chain", "locus"])),
            "v_gene": _as_text(_pull(df, ["v_gene", "v_call", "v", "v_gene_call"])),
            "d_gene": _as_text(_pull(df, ["d_gene", "d_call", "d", "d_gene_call"])),
            "j_gene": _as_text(_pull(df, ["j_gene", "j_call", "j", "j_gene_call"])),
            "c_gene": _as_text(_pull(df, ["c_gene", "c_call", "c", "c_gene_call"])),
            "cdr3_nt": _as_text(_pull(df, ["cdr3_nt", "junction", "junction_nt"])),
            "cdr3_aa": _as_text(_pull(df, ["cdr3", "cdr3_aa", "junction_aa"])),
            "productive": _as_bool(_pull(df, ["productive"], default=True)),
            "full_length": _as_bool(_pull(df, ["full_length", "high_confidence"], default=True)),
            "is_cell": _as_bool(_pull(df, ["is_cell"], default=True)),
            "umi_count": _as_number(_pull(df, ["umis", "umi_count", "consensus_count"])),
            "read_count": _as_number(_pull(df, ["reads", "read_count"])),
        }
    )

    filtered = std[
        std["cell_barcode"].notna()
        & std["chain"].isin(["TRA", "TRB"])
        & (std["productive"].isna() | (std["productive"] == True))
        & (std["full_length"].isna() | (std["full_length"] == True))
        & (std["is_cell"].isna() | (std["is_cell"] == True))
    ].copy()

    if filtered.empty:
        return pd.DataFrame(columns=["cell_barcode"])

    chain_counts = (
        filtered.groupby(["cell_barcode", "chain"], dropna=False)
        .size()
        .unstack(fill_value=0)
        .rename(columns={"TRA": "n_tra_productive", "TRB": "n_trb_productive"})
    )
    for column in ["n_tra_productive", "n_trb_productive"]:
        if column not in chain_counts.columns:
            chain_counts[column] = 0

    dominant = filtered.assign(
        _umi_sort=filtered["umi_count"].fillna(-1),
        _read_sort=filtered["read_count"].fillna(-1),
    ).sort_values(
        [
            "cell_barcode",
            "chain",
            "_umi_sort",
            "_read_sort",
            "v_gene",
            "j_gene",
            "c_gene",
            "cdr3_nt",
            "cdr3_aa",
        ],
        ascending=[True, True, False, False, True, True, True, True, True],
        na_position="last",
    )
    dominant = dominant.groupby(["cell_barcode", "chain"], dropna=False).head(1)

    tra = dominant[dominant["chain"] == "TRA"].rename(
        columns={
            "v_gene": "trav",
            "j_gene": "traj",
            "c_gene": "trac",
            "cdr3_nt": "tra_cdr3_nt",
            "cdr3_aa": "tra_cdr3_aa",
            "umi_count": "tra_umis",
            "read_count": "tra_reads",
        }
    )[
        ["cell_barcode", "trav", "traj", "trac", "tra_cdr3_nt", "tra_cdr3_aa", "tra_umis", "tra_reads"]
    ]
    trb = dominant[dominant["chain"] == "TRB"].rename(
        columns={
            "v_gene": "trbv",
            "d_gene": "trbd",
            "j_gene": "trbj",
            "c_gene": "trbc",
            "cdr3_nt": "trb_cdr3_nt",
            "cdr3_aa": "trb_cdr3_aa",
            "umi_count": "trb_umis",
            "read_count": "trb_reads",
        }
    )[
        ["cell_barcode", "trbv", "trbd", "trbj", "trbc", "trb_cdr3_nt", "trb_cdr3_aa", "trb_umis", "trb_reads"]
    ]

    wide = tra.merge(trb, on="cell_barcode", how="outer")
    wide = wide.merge(chain_counts.reset_index(), on="cell_barcode", how="left")
    wide["n_tra_productive"] = wide["n_tra_productive"].fillna(0).astype(int)
    wide["n_trb_productive"] = wide["n_trb_productive"].fillna(0).astype(int)
    wide["has_multi_tra"] = wide["n_tra_productive"] > 1
    wide["has_multi_trb"] = wide["n_trb_productive"] > 1

    constants = {
        "dataset_id": dataset_id,
        "cancer_type": cancer_type,
        "donor_id": donor_id,
        "sample_id": sample_id,
        "tissue_type": tissue_type,
    }
    wide = _apply_metadata(wide, {k: v for k, v in constants.items() if v is not None})
    wide = _apply_metadata(wide, sample_metadata)
    return wide
