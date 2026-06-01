# TCR-CLAIM Method Specification

**TCR-CLAIM** means **TCR Candidate-Level Assessment and Interpretation Method**.

The method is built around one rule:

> A relaxed TCR group can prioritize biological candidates, but a final clone
> claim requires paired strict CDR3 evidence.

## Evidence Levels

| Evidence level | Meaning | Allowed claim |
|---|---|---|
| `confirmed_clone` | Paired TRA/TRB strict CDR3 clonotype | Operational TCR clonotype |
| `dominant_vgene_group` | `TRAV-TRBV` group dominated by one strict clone | Strong candidate group |
| `mixed_vgene_group` | `TRAV-TRBV` group with multiple strict clones | Ambiguous candidate group |
| `strict_sharing` | Paired CDR3 clone appears in both tissues | Defensible clonal sharing |
| `apparent_sharing` | `TRAV-TRBV` is shared without strict CDR3 support | Relaxed sharing candidate only |

## Required Canonical Tables

TCR-CLAIM is table-first. It can operate on AnnData later, but the MVP requires
these tables:

- `cell_tcr_table`: one row per cell with dominant TRA/TRB and clone keys.
- `qc_summary`: paired-TCR eligibility, CD4/CD8 primary universe and multi-chain
  ambiguity per context.
- `strict_clone_table`: one row per strict clone and context.
- `vgene_group_table`: one row per `TRAV-TRBV` group and context.
- `strict_clonal_diversity`: strict paired-CDR3 repertoire diversity per context.
- `relaxed_clonal_diversity`: relaxed `TRAV-TRBV` repertoire diversity per context.
- `strict_vs_relaxed_diversity`: compression created by using relaxed groups.
- `sharing_table`: one row per tissue pair, context, and threshold.
- `candidate_table`: prioritized strict clone candidates for biological review.
- `candidate_phenotype_table`: descriptive phenotype score shifts for each
  prioritized candidate when score columns are available.
- `claim_table`: one row per auditable biological claim.
- `xenium_panel_roadmap`: spatial validation targets informed by default immune
  signatures and candidate-supported states.
- `xenium_cdr3_targets.fasta`: optional CDR3 target sequences when
  `--include-xenium-cdr3` is explicitly requested.
- `candidate_cards.html`: candidate-centered evidence cards for biological
  review.
- `tcr_claim_report.md`: human-readable run summary.

## Primary Benchmark Universe

The default TCR-CLAIM metrics mirror the benchmark primary universe:

- `cell_class` is `CD4` or `CD8` when the column exists;
- `ct_strict` is present;
- `ct_vgene` is present.

Use the CLI flag `--all-cells` only for exploratory QC when non-CD4/CD8 cells
or missing cell-class annotations should be retained.

## Claim Policy

The method explicitly blocks unsupported overclaims.

Allowed:

- "expanded CD8 tumor-associated strict clonotype";
- "relaxed `TRAV-TRBV` candidate group";
- "strict CDR3-backed tissue sharing";
- "spatially testable T-cell state candidate".

Not allowed without additional validation:

- "antigen-specific clone";
- "tumor-reactive clone";
- "functional killing clone";
- "confirmed shared clonotype" for relaxed-only sharing.

## Candidate Ranking

`candidate_table` is a triage table, not a tumor-reactivity model.

The current score prioritizes:

- larger strict paired-CDR3 clones;
- tumor and CD8 contexts when `focus="tumor_cd8"`;
- relaxed groups dominated by the candidate strict clone;
- lower collapse risk.

The score is intentionally heuristic. It should be used to decide which clones
deserve biological review, not to infer antigen specificity or function.

## Candidate Phenotype Table

`candidate_phenotype_table` connects TCR candidates to scRNA-derived phenotype
scores when those scores are present in the input table or passed with
`--phenotype-scores`.

For each prioritized strict clone candidate and score column, it reports:

- candidate cells in the same biological context;
- reference cells from the same context excluding that candidate;
- candidate mean and median;
- reference mean and median;
- `delta_mean`;
- `direction`: `enriched`, `depleted`, `no_clear_shift` or `not_evaluable`.

If no phenotype scores are available, the table still records
`phenotype_evidence_status = no_scores_available`. This is intentional: absence
of phenotype evidence should be explicit in downstream reports.

Phenotype enrichment is descriptive. It does not validate antigen specificity,
tumor reactivity or functional killing.

## Multi-Dataset Runner

`tcr-claim-batch` runs the same per-dataset workflow across many benchmark
result directories. It discovers folders containing `cell_metadata_with_tcr.csv`
under `--results-root`, or accepts explicit folders with `--results`.

Batch outputs:

- `batch_run_summary.csv`: one row per dataset with QC, clone counts, diversity
  compression, sharing, phenotype-evidence status and candidate counts;
- `batch_report.md`: compact Markdown index for manuscript and supplement
  triage;
- `supplement_tables/`: manuscript-ready aggregate tables;
- `figures/`: lightweight SVG figures generated from batch summaries;
- `per_dataset/<result_id>/`: full TCR-CLAIM outputs for each dataset.

By default, the batch runner continues if a dataset fails and records a `fail`
row. Use `--fail-on-error` when the desired behavior is strict CI-style failure.

The CLI also applies a default `--max-input-rows 2000000` guard. Datasets above
that size are recorded as `skip_large_input` rather than silently producing
multi-GB intermediate files. Use `--max-input-rows 0` only for intentional large
dataset runs or after implementing a streaming/pre-filtered workflow.

## Dataset Audit And Large-Input Prefilter

`tcr-claim-audit` audits benchmark result folders in chunks and reports:

- row count;
- presence of `cell_class`, `ct_strict` and `ct_vgene`;
- paired-TCR count;
- primary CD4/CD8 paired-TCR count;
- `cell_class` label distribution;
- dataset status such as `ok`, `no_paired_tcr` or `zero_primary_universe`.

`tcr-claim-prefilter` writes a reduced `cell_metadata_with_tcr.csv` in chunks.
By default it keeps rows with both `ct_strict` and `ct_vgene` and `cell_class`
in `CD4,CD8`. This is the preferred path for very large inputs such as
`gse200996`, where writing an unfiltered full `cell_tcr_table.csv` can create
multi-GB intermediates.

## Supplemental Tables

The supplement generator turns a batch output root into:

- `supp_dataset_qc.csv`: dataset-level QC, clone counts, diversity compression,
  sharing counts and phenotype-evidence status;
- `supp_collapse_risk.csv`: all relaxed `TRAV-TRBV` groups and collapse-risk
  metrics;
- `supp_diversity_compression.csv`: strict versus relaxed diversity metrics per
  context;
- `supp_sharing_apparent_vs_strict.csv`: strict-backed versus apparent-only
  sharing per tissue pair;
- `supp_candidate_index.csv`: prioritized candidate index across datasets;
- `supp_claim_inventory.csv`: claim counts by entity type and evidence level.

These are descriptive benchmark tables. They are not antigen-specificity or
functional-validation results.

## Reproducible Figures

The figure generator writes SVGs directly from batch outputs:

- relaxed/strict richness ratio by dataset;
- apparent relaxed-only sharing by dataset;
- candidate count by dataset;
- collapse-risk distribution;
- TCR-CLAIM workflow schematic.

The figures are intentionally lightweight and dependency-free. They are meant as
first-pass manuscript/supplement figures that can later be restyled for journal
submission.

## Diversity Compression

`strict_vs_relaxed_diversity` quantifies how much repertoire structure is lost
when strict paired-CDR3 clonotypes are represented as relaxed `TRAV-TRBV` groups.

Key fields:

- `richness_ratio`: relaxed richness divided by strict richness;
- `richness_relative_difference`: relative change in clone count;
- `effective_shannon_ratio`: relaxed/strict effective Shannon diversity;
- `gini_difference`: change in clone-size inequality.

Low relaxed/strict richness ratios indicate strong compression and should
increase caution around relaxed-group biological claims.

## Xenium Roadmap

TCR-CLAIM treats Xenium as a validation roadmap:

1. discover strict clones and candidate groups in scRNA/scTCR;
2. nominate state genes and optional CDR3 targets;
3. design a spatial panel;
4. test spatial location, neighborhoods, and tumor proximity;
5. keep antigen specificity separate from spatial support.

When candidate phenotype scores are available, `xenium_panel_roadmap` marks
matching state signatures as `candidate_supported`. Otherwise, it remains a
default immune/spatial roadmap and does not claim candidate-state support.

CDR3 FASTA export is opt-in. It should be treated as an experimental advanced
custom design input, not as evidence that Xenium will detect the clone or that
the clone is tumor-reactive.
