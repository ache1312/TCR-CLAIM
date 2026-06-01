"""Auditable biological interpretation of single-cell TCR clone candidates."""

from .io import read_10x_contigs, read_table
from .normalize import normalize_contigs
from .definitions import define_clones
from .filters import primary_tcr_cells
from .metrics import clone_count_agreement, collapse_risk, tissue_sharing
from .qc import qc_summary
from .diversity import clonal_diversity, strict_vs_relaxed_diversity
from .prioritize import prioritize_candidates
from .tables import cell_tcr_table, strict_clone_table
from .phenotype import candidate_phenotype_table, clone_phenotype_association, infer_phenotype_score_columns
from .claims import claim_checker
from .reporting import report_batch_summary, report_candidate_cards, report_clone_cards, report_markdown_summary
from .pipeline import discover_result_dirs, run_tcr_claim_batch, run_tcr_claim_pipeline, summarize_tcr_claim_outputs
from .xenium import design_xenium_panel_from_candidates, export_cdr3_fasta_for_xenium

__all__ = [
    "read_table",
    "read_10x_contigs",
    "normalize_contigs",
    "define_clones",
    "primary_tcr_cells",
    "qc_summary",
    "clonal_diversity",
    "strict_vs_relaxed_diversity",
    "prioritize_candidates",
    "collapse_risk",
    "clone_count_agreement",
    "tissue_sharing",
    "cell_tcr_table",
    "strict_clone_table",
    "candidate_phenotype_table",
    "clone_phenotype_association",
    "infer_phenotype_score_columns",
    "claim_checker",
    "discover_result_dirs",
    "run_tcr_claim_batch",
    "run_tcr_claim_pipeline",
    "summarize_tcr_claim_outputs",
    "report_clone_cards",
    "report_candidate_cards",
    "report_markdown_summary",
    "report_batch_summary",
    "design_xenium_panel_from_candidates",
    "export_cdr3_fasta_for_xenium",
]
