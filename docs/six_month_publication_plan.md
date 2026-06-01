# Six-Month Publication Plan

## Goal

Publish **TCR-CLAIM** as a translational cancer immunology method with companion
software in Python and R.

Final month-6 deliverable:

- preprint;
- journal submission;
- `tcrbio` Python package release;
- R bridge for benchmark interoperability;
- public reproducibility bundle;
- clone-card report demonstrating claims that are allowed and not allowed.

## Month 1: Benchmark Lock

- Freeze `results/cross_dataset_with_io_direct_source_20260504` as the primary
  benchmark result set.
- Audit datasets with zero primary CD4/CD8 paired-TCR cells.
- Audit `GSE200996` raw-cell versus TCR-paired-cell denominator.
- Finalize figure list and manuscript outline.
- Define final terminology:
  - strict paired CDR3 clone;
  - relaxed `TRAV-TRBV` group;
  - collapse risk;
  - strict sharing;
  - apparent sharing;
  - allowed claim.

## Month 2: Software MVP

- Implement `tcrbio` core API:
  - normalization;
  - clone definitions;
  - strict clone table;
  - V-gene collapse risk;
  - tissue sharing;
  - phenotype association;
  - claim checking;
  - clone-card reports.
- Implement R bridge for table exchange with the current benchmark.
- Add synthetic tests matching the benchmark logic.
- Add CLI to process existing `cell_metadata_with_tcr.csv` files.

## Month 3: Reproducibility Bundle

- Convert public and private datasets into a documented reproducibility manifest.
- Prepare public release of IO/breast if permissions allow complete release.
- Generate canonical tables:
  - `cell_tcr_table.csv`;
  - `strict_clone_table.csv`;
  - `vgene_group_table.csv`;
  - `sharing_table.csv`;
  - `claim_table.csv`.
- Validate Python outputs against R benchmark outputs.

## Month 4: Xenium Roadmap

- Generate a Xenium panel roadmap from candidate clone/state results.
- Include core T cell, cytotoxicity, exhaustion, progenitor/memory, Treg,
  antigen-presentation, IFN, tumor, stromal, myeloid and endothelial markers.
- Export optional CDR3 FASTA candidates as experimental advanced custom design
  targets.
- Keep direct CDR3 spatial detection as an experimental validation layer, not a
  primary claim.

## Month 5: Manuscript And Release Candidate

- Write full manuscript.
- Produce final figures and supplemental tables.
- Prepare `tcrbio v0.1.0` release candidate.
- Write documentation:
  - quickstart;
  - method specification;
  - claim policy;
  - benchmark reproduction;
  - Xenium roadmap.
- Prepare citation, license and archive metadata.

## Month 6: Preprint And Submission

- Internal review by TCR biology, cancer immunology, software and reproducibility
  reviewers.
- Freeze code and data versions.
- Archive software and data.
- Upload preprint.
- Submit manuscript.

## Publication Positioning

TCR-CLAIM should not be positioned as a new universal clonotype definition.

It should be positioned as:

> a framework for turning single-cell TCR repertoires into biologically
> defensible cancer immune clone candidates by separating confirmed clones,
> relaxed candidate groups, collapse risk, apparent sharing and unsupported
> biological claims.
