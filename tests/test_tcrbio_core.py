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


def test_xenium_panel_marks_candidate_supported_state_genes():
    phenotypes = pd.DataFrame(
        {
            "candidate_rank": [1],
            "candidate_id": ["strict_big"],
            "phenotype_state": ["cytotoxic"],
            "direction": ["enriched"],
        }
    )

    panel = tb.design_xenium_panel_from_candidates(candidate_phenotypes=phenotypes)
    gzmb = panel[panel["target"] == "GZMB"].iloc[0]

    assert "candidate_state:cytotoxic" in gzmb["source"]
    assert gzmb["priority"] in {"candidate_supported", "roadmap"}


def test_xenium_cdr3_fasta_uses_candidate_table_sequences(tmp_path):
    candidates = pd.DataFrame(
        {
            "candidate_id": ["strict_big"],
            "trb_cdr3_nt": ["TGTGCCAGC"],
        }
    )
    output = tmp_path / "cdr3.fasta"

    tb.export_cdr3_fasta_for_xenium(candidates, output)

    assert output.read_text(encoding="utf-8") == ">strict_big\nTGTGCCAGC\n"


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


def test_qc_summary_reports_primary_eligibility():
    cells = pd.DataFrame(
        {
            "dataset_id": "test",
            "donor_id": "d1",
            "sample_id": "s1",
            "tissue_type": "tumor",
            "cell_class": ["CD8", "CD4", "NK"],
            "ct_strict": ["strict_a", "strict_b", "strict_c"],
            "ct_vgene": ["TRAV8_TRBV13", "TRAV1_TRBV1", "TRAV2_TRBV2"],
            "has_multi_tra": [False, True, False],
            "has_multi_trb": [False, False, True],
        }
    )

    qc = tb.qc_summary(cells)
    row = qc.iloc[0]

    assert row["n_cells"] == 3
    assert row["n_primary_cd4_cd8_paired"] == 2
    assert row["multi_tra_rate"] == 1 / 3


def test_strict_vs_relaxed_diversity_captures_compression():
    cells = pd.DataFrame(
        {
            "dataset_id": "test",
            "cancer_type": "breast",
            "tissue_type": "tumor",
            "cell_class": "CD8",
            "ct_strict": ["strict_a", "strict_b", "strict_c"],
            "ct_vgene": ["TRAV8_TRBV13", "TRAV8_TRBV13", "TRAV1_TRBV1"],
        }
    )

    diversity = tb.strict_vs_relaxed_diversity(cells)
    row = diversity.iloc[0]

    assert row["richness_strict"] == 3
    assert row["richness_relaxed"] == 2
    assert row["richness_relative_difference"] == -1 / 3


def test_strict_vs_relaxed_diversity_handles_empty_primary_universe():
    cells = pd.DataFrame(
        {
            "cell_class": ["NK"],
            "ct_strict": ["strict_a"],
            "ct_vgene": ["TRAV8_TRBV13"],
        }
    )

    diversity = tb.strict_vs_relaxed_diversity(cells)

    assert diversity.empty


def test_prioritize_candidates_prefers_tumor_cd8_large_low_risk_clone():
    cells = pd.DataFrame(
        {
            "dataset_id": "test",
            "cancer_type": "breast",
            "donor_id": "d1",
            "sample_id": "s1",
            "tissue_type": ["tumor", "tumor", "blood"],
            "cell_class": ["CD8", "CD8", "CD4"],
            "ct_strict": ["strict_big", "strict_big", "strict_other"],
            "ct_vgene": ["TRAV8_TRBV13", "TRAV8_TRBV13", "TRAV1_TRBV1"],
        }
    )

    candidates = tb.prioritize_candidates(cells)

    assert candidates.iloc[0]["ct_strict"] == "strict_big"
    assert candidates.iloc[0]["candidate_rank"] == 1


def test_prioritize_candidates_handles_missing_tissue_context():
    cells = pd.DataFrame(
        {
            "ct_strict": ["strict_big", "strict_big", "strict_other"],
            "ct_vgene": ["TRAV8_TRBV13", "TRAV8_TRBV13", "TRAV1_TRBV1"],
        }
    )

    candidates = tb.prioritize_candidates(cells)

    assert candidates.iloc[0]["ct_strict"] == "strict_big"
    assert not bool(candidates.iloc[0]["tumor_context"])


def test_candidate_phenotype_table_reports_contextual_score_shift():
    cells = pd.DataFrame(
        {
            "dataset_id": "test",
            "cancer_type": "breast",
            "donor_id": "d1",
            "sample_id": "s1",
            "tissue_type": "tumor",
            "cell_class": "CD8",
            "ct_strict": ["strict_big", "strict_big", "strict_other", "strict_other"],
            "ct_vgene": ["TRAV8_TRBV13", "TRAV8_TRBV13", "TRAV1_TRBV1", "TRAV1_TRBV1"],
            "cytotoxic_score": [2.0, 2.2, 0.5, 0.7],
        }
    )
    candidates = tb.prioritize_candidates(cells)

    phenotypes = tb.candidate_phenotype_table(cells, candidates, scores=["cytotoxic_score"])
    row = phenotypes[phenotypes["candidate_id"] == "strict_big"].iloc[0]

    assert row["phenotype_state"] == "cytotoxic"
    assert row["direction"] == "enriched"
    assert round(row["delta_mean"], 2) == 1.5


def test_candidate_phenotype_table_is_explicit_when_scores_are_absent():
    cells = pd.DataFrame(
        {
            "ct_strict": ["strict_big", "strict_big"],
            "ct_vgene": ["TRAV8_TRBV13", "TRAV8_TRBV13"],
        }
    )
    candidates = tb.prioritize_candidates(cells)

    phenotypes = tb.candidate_phenotype_table(cells, candidates)

    assert phenotypes.iloc[0]["phenotype_evidence_status"] == "no_scores_available"
    assert phenotypes.iloc[0]["direction"] == "not_evaluable"


def test_markdown_report_summarizes_candidates_and_claim_boundary(tmp_path):
    report_path = tmp_path / "report.md"
    qc = pd.DataFrame({"n_cells": [3], "n_paired_tcr": [2], "n_primary_cd4_cd8_paired": [2]})
    candidates = pd.DataFrame(
        {
            "candidate_rank": [1],
            "candidate_id": ["strict_big_clone_identifier"],
            "n_cells": [10],
            "tissue_type": ["tumor"],
            "cell_class": ["CD8"],
            "risk_label": ["low"],
            "rank_score": [120.0],
        }
    )
    claims = pd.DataFrame(
        {
            "entity_type": ["prioritized_candidate"],
            "evidence_level": ["confirmed_clone"],
        }
    )

    markdown = tb.report_markdown_summary(output=report_path, qc=qc, candidates=candidates, claims=claims)

    assert "TCR-CLAIM Run Summary" in markdown
    assert "Not supported by TCR-CLAIM alone" in markdown
    assert "strict_big_clone_identifier" in markdown
    assert report_path.exists()
