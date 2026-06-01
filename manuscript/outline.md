# Manuscript Outline: TCR-CLAIM

Working title:

**TCR-CLAIM converts single-cell TCR repertoires into biologically auditable cancer immune clone candidates**

## Central Claim

`TRAV-TRBV` is useful for prioritizing expanded tumor CD8 TCR groups, but it
compresses repertoire diversity and produces apparent sharing. TCR-CLAIM turns
this distinction into an explicit evidence framework for confirmed clones,
candidate groups, and allowed biological claims.

## Main Figures

1. Conceptual overview of strict paired CDR3 clones versus relaxed V-gene groups.
2. Cross-dataset clone-count agreement by threshold and cell class.
3. Repertoire diversity compression under `TRAV-TRBV`.
4. Strict-backed versus apparent tissue sharing.
5. Phenotype associations of expanded tumor CD8 strict clones.
6. TCR-CLAIM software workflow and Xenium validation roadmap.

## Key Results To Report

- 10 result sets, 380,987 paired-TCR cells, 174,118 primary CD4/CD8 paired-TCR cells.
- Tumor CD8 top-10 strict clones: median relative difference 0.0% under `TRAV-TRBV`.
- Tumor CD8 all-clone repertoire: median relative difference -28.5%, showing diversity compression.
- Tumor CD8 `TRAV-TRBV` tissue sharing precision-like around 55-65%, depending on threshold.
- Tumor CD4 relaxed grouping is less stable, especially for expanded clone thresholds.
- Expanded tumor CD8 strict clones show higher cytotoxic/exhausted signatures and lower progenitor-like signatures.

## Methods Sections

- Dataset intake and eligibility.
- Dominant TRA/TRB normalization.
- Strict and relaxed TCR definitions.
- Collapse-risk metrics.
- Tissue-sharing metrics.
- Phenotype-clonality association.
- MAIT-like semi-invariant control.
- TCR-CLAIM software and claim policy.
- Xenium panel-design roadmap.

## Discussion Points

- `TRAV-TRBV` is a prioritization layer, not a clonotype replacement.
- Precision-like is more informative than recall for relaxed tissue sharing.
- Claims about antigen specificity require external validation.
- The TCR-CLAIM framework is useful for translational prioritization and spatial validation design.
