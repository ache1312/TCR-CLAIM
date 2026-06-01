import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import tcrbio as tb


def main() -> None:
    contigs = pd.DataFrame(
        [
            ["cell1", "TRA", "TRAV8-2", None, "TRAJ18", "TRAC", "AAA", "CAVA", True, 10, 100],
            ["cell1", "TRB", "TRBV13", None, "TRBJ2-1", "TRBC2", "BBB", "CASB", True, 12, 120],
            ["cell2", "TRA", "TRAV8-2", None, "TRAJ18", "TRAC", "AAA", "CAVA", True, 8, 80],
            ["cell2", "TRB", "TRBV13", None, "TRBJ2-1", "TRBC2", "CCC", "CASC", True, 8, 80],
            ["cell3", "TRA", "TRAV1-2", None, "TRAJ33", "TRAC", "DDD", "CAVD", True, 8, 80],
            ["cell3", "TRB", "TRBV20-1", None, "TRBJ1-1", "TRBC1", "EEE", "CASE", True, 8, 80],
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

    metadata = pd.DataFrame(
        {
            "cell_barcode": ["cell1", "cell2", "cell3"],
            "dataset_id": ["demo"] * 3,
            "cancer_type": ["breast"] * 3,
            "donor_id": ["d1"] * 3,
            "sample_id": ["tumor"] * 3,
            "tissue_type": ["tumor"] * 3,
            "cell_class": ["CD8"] * 3,
        }
    )

    tcr = tb.normalize_contigs(contigs, sample_metadata=metadata)
    tcr = tb.define_clones(tcr)
    risk = tb.collapse_risk(tcr)
    claims = tb.claim_checker(tcr=tcr, risk=risk)
    panel = tb.design_xenium_panel_from_candidates()

    out_dir = Path("demo_outputs")
    out_dir.mkdir(exist_ok=True)
    tcr.to_csv(out_dir / "cell_tcr_table.csv", index=False)
    risk.to_csv(out_dir / "vgene_group_table.csv", index=False)
    claims.to_csv(out_dir / "claim_table.csv", index=False)
    panel.to_csv(out_dir / "xenium_panel_roadmap.csv", index=False)
    tb.report_clone_cards(claims, output=out_dir / "clone_cards.html")

    print(f"Wrote demo outputs to {out_dir.resolve()}")


if __name__ == "__main__":
    main()
