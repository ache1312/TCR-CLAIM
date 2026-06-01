import pandas as pd

import tcrbio as tb


def test_normalize_and_define_dominant_paired_chains():
    contigs = pd.DataFrame(
        [
            ["cell1-1", "TRA", "TRAV8-2", None, "TRAJ18", "TRAC", "AAA", "CAVA", True, 10, 100],
            ["cell1-1", "TRA", "TRAV1-1", None, "TRAJ1", "TRAC", "LOW", "CLOW", True, 1, 10],
            ["cell1-1", "TRB", "TRBV13", None, "TRBJ2-1", "TRBC2", "BBB", "CASB", True, 12, 120],
            ["cell2-1", "TRA", "TRAV8-2", None, "TRAJ18", "TRAC", "AAA", "CAVA", True, 5, 50],
            ["cell2-1", "TRB", "TRBV13", None, "TRBJ2-1", "TRBC2", "CCC", "CASC", True, 5, 50],
            ["cell3-1", "TRA", "TRAV8-2", None, "TRAJ18", "TRAC", "AAA", "CAVA", True, 5, 50],
        ],
        columns=[
            "barcode",
            "chain",
            "v_gene",
            "d_gene",
            "j_gene",
            "c_gene",
            "cdr3_nt",
            "cdr3",
            "productive",
            "umis",
            "reads",
        ],
    )

    tcr = tb.normalize_contigs(contigs, dataset_id="test", sample_id="s1")
    tcr = tb.define_clones(tcr)

    assert len(tcr) == 3
    assert tcr.loc[tcr["cell_barcode"] == "cell1-1", "trav"].iloc[0] == "TRAV8-2"
    assert bool(tcr.loc[tcr["cell_barcode"] == "cell1-1", "has_multi_tra"].iloc[0])
    assert not bool(tcr.loc[tcr["cell_barcode"] == "cell3-1", "has_paired_tcr"].iloc[0])
    paired = tcr[tcr["cell_barcode"].isin(["cell1-1", "cell2-1"])]
    assert set(paired["ct_vgene"]) == {"TRAV8-2_TRBV13"}
    assert paired["ct_strict"].nunique() == 2


def test_collapse_risk_labels_mixed_groups():
    cells = pd.DataFrame(
        {
            "dataset_id": "test",
            "donor_id": "d1",
            "sample_id": "s1",
            "tissue_type": "tumor",
            "cell_class": "CD8",
            "ct_strict": ["strict_a", "strict_a", "strict_b", "strict_c"],
            "ct_vgene": ["TRAV8_TRBV13", "TRAV8_TRBV13", "TRAV8_TRBV13", "TRAV1_TRBV1"],
        }
    )
    risk = tb.collapse_risk(cells)
    merged = risk[risk["ct_vgene"] == "TRAV8_TRBV13"].iloc[0]
    singleton = risk[risk["ct_vgene"] == "TRAV1_TRBV1"].iloc[0]

    assert merged["n_strict_clones"] == 2
    assert merged["dominant_fraction"] == 2 / 3
    assert merged["risk_label"] == "high"
    assert singleton["risk_label"] == "low"


def test_tissue_sharing_separates_strict_backed_from_apparent():
    cells = pd.DataFrame(
        {
            "dataset_id": "test",
            "donor_id": "d1",
            "sample_id": ["tumor", "blood", "tumor", "blood"],
            "tissue_type": ["tumor", "peripheral_blood", "tumor", "peripheral_blood"],
            "cell_class": "CD8",
            "ct_strict": ["strict_shared", "strict_shared", "strict_tumor_only", "strict_blood_only"],
            "ct_vgene": ["TRAV8_TRBV13", "TRAV8_TRBV13", "TRAV1_TRBV1", "TRAV1_TRBV1"],
        }
    )

    sharing = tb.tissue_sharing(cells, thresholds=[1])
    row = sharing.iloc[0]

    assert row["strict_shared_clones"] == 1
    assert row["vgene_shared_groups"] == 2
    assert row["vgene_shared_strict_backed"] == 1
    assert row["vgene_shared_apparent_only"] == 1
    assert row["vgene_sharing_precision_like"] == 0.5


def test_claim_checker_blocks_overclaims():
    cells = pd.DataFrame(
        {
            "ct_strict": ["strict_a", "strict_a", "strict_b"],
            "ct_vgene": ["TRAV8_TRBV13", "TRAV8_TRBV13", "TRAV8_TRBV13"],
        }
    )
    risk = tb.collapse_risk(cells, groupby=[])
    claims = tb.claim_checker(tcr=cells, risk=risk)

    strict_claim = claims[claims["entity_type"] == "strict_clone"].iloc[0]
    relaxed_claim = claims[claims["entity_type"] == "relaxed_vgene_group"].iloc[0]

    assert "Operational paired CDR3" in strict_claim["allowed_claim"]
    assert "Antigen-specific" in strict_claim["not_allowed_claim"]
    assert "Ambiguous" in relaxed_claim["allowed_claim"]


def test_strict_clone_table_summarizes_confirmed_clones():
    cells = pd.DataFrame(
        {
            "dataset_id": "test",
            "donor_id": "d1",
            "sample_id": "s1",
            "tissue_type": "tumor",
            "cell_class": "CD8",
            "ct_strict": ["strict_a", "strict_a", "strict_b"],
            "ct_vgene": ["TRAV8_TRBV13", "TRAV8_TRBV13", "TRAV1_TRBV1"],
        }
    )

    strict = tb.strict_clone_table(cells)

    assert set(strict["ct_strict"]) == {"strict_a", "strict_b"}
    assert strict.loc[strict["ct_strict"] == "strict_a", "n_cells"].iloc[0] == 2
    assert set(strict["evidence_level"]) == {"confirmed_clone"}


def test_xenium_panel_design_has_core_markers():
    panel = tb.design_xenium_panel_from_candidates()
    assert {"CD3D", "CD8A", "PDCD1", "TOX", "HLA-A"}.issubset(set(panel["target"]))


def test_primary_metrics_exclude_non_cd4_cd8_cells():
    cells = pd.DataFrame(
        {
            "dataset_id": "test",
            "donor_id": "d1",
            "sample_id": "s1",
            "tissue_type": "tumor",
            "cell_class": ["CD8", "NK"],
            "ct_strict": ["strict_a", "strict_b"],
            "ct_vgene": ["TRAV8_TRBV13", "TRAV1_TRBV1"],
        }
    )

    primary = tb.primary_tcr_cells(cells)
    agreement = tb.clone_count_agreement(cells, thresholds=[1], include_top10=False)

    assert primary.shape[0] == 1
    assert agreement["strict_clone_count"].iloc[0] == 1
