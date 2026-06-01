"""Auditable biological interpretation of single-cell TCR clone candidates."""

from .io import read_10x_contigs, read_table
from .normalize import normalize_contigs
from .definitions import define_clones
from .filters import primary_tcr_cells
from .metrics import clone_count_agreement, collapse_risk, tissue_sharing
from .tables import cell_tcr_table, strict_clone_table
from .phenotype import clone_phenotype_association
from .claims import claim_checker
from .reporting import report_clone_cards
from .xenium import design_xenium_panel_from_candidates, export_cdr3_fasta_for_xenium

__all__ = [
    "read_table",
    "read_10x_contigs",
    "normalize_contigs",
    "define_clones",
    "primary_tcr_cells",
    "collapse_risk",
    "clone_count_agreement",
    "tissue_sharing",
    "cell_tcr_table",
    "strict_clone_table",
    "clone_phenotype_association",
    "claim_checker",
    "report_clone_cards",
    "design_xenium_panel_from_candidates",
    "export_cdr3_fasta_for_xenium",
]
