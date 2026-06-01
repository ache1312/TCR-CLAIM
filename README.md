# tcrbio

`tcrbio` implements the **TCR-CLAIM** method: TCR Candidate-Level Assessment and
Interpretation Method. It is a small MVP package for turning single-cell TCR data
into auditable biological clone candidates. It is designed around the conclusion
of the clonotype benchmark:

> `TRAV-TRBV` is useful for prioritizing candidate TCR groups, especially tumor
> CD8 expanded groups, but it does not replace paired CDR3-based strict
> clonotypes.

The package keeps strict clonotypes, relaxed V-gene groups, collapse risk,
sharing calls, diversity compression, candidate prioritization, phenotype
associations, and claim checking separate.

## Quick Start

```python
import tcrbio as tb

contigs = tb.read_10x_contigs("filtered_contig_annotations.csv")

tcr = tb.normalize_contigs(
    contigs,
    dataset_id="dataset",
    donor_id="donor_1",
    sample_id="tumor_1",
    tissue_type="tumor",
)

tcr = tb.define_clones(tcr)
risk = tb.collapse_risk(tcr)
sharing = tb.tissue_sharing(tcr)
claims = tb.claim_checker(tcr=tcr, risk=risk, sharing=sharing)

tb.report_clone_cards(claims, output="clone_cards.html")
```

## Current MVP Modules

- `io`: read 10x/AIRR-like CSV/TSV tables.
- `normalize`: choose dominant productive TRA/TRB per cell and keep ambiguity
  flags.
- `definitions`: build `ct_strict` and `ct_vgene` with explicit evidence levels.
- `qc`: summarize paired-TCR eligibility and multi-chain ambiguity.
- `diversity`: compare strict versus relaxed repertoire diversity.
- `metrics`: compute collapse risk, clone count agreement, and strict/apparent
  tissue sharing.
- `prioritize`: rank strict candidate clones for biological review.
- `claims`: generate allowed and not-allowed biological claim statements.
- `phenotype`: associate strict clone expansion with phenotype scores.
- `xenium`: build a candidate Xenium marker/CDR3 panel specification.
- `reporting`: render HTML clone cards and Markdown run summaries.
- `r/tcrbio_bridge.R`: lightweight R helpers for exchanging canonical tables.

## Method Name

**TCR-CLAIM** separates what a TCR result can support from what it cannot support:

- confirmed paired-CDR3 clonotypes;
- relaxed `TRAV-TRBV` candidate groups;
- low/medium/high V-gene collapse risk;
- strict versus apparent tissue sharing;
- strict versus relaxed diversity compression;
- prioritized clone candidates;
- allowed versus unsupported biological claims.

## Run Tests

```bash
python -m pytest
```

## CLI

Generate TCR-CLAIM tables from an existing benchmark cell metadata table:

```bash
python scripts/run_tcr_claim_tables.py \
  --input /path/to/cell_metadata_with_tcr.csv \
  --out tcr_claim_outputs
```

Validate Python metrics against R benchmark outputs:

```bash
python scripts/validate_against_r_benchmark.py \
  --results-dir /path/to/results/io_dataset \
  --out validation/io_dataset_validation_summary.csv
```

After installing the package, the equivalent entry points are:

```bash
tcr-claim-tables --input cell_metadata_with_tcr.csv --out tcr_claim_outputs
tcr-claim-validate --results-dir results/io_dataset --out validation.csv
tcr-claim-validate-batch --results results/io_dataset,results/gse121637 --out batch_validation.csv
```

The table CLI writes:

- `cell_tcr_table.csv`
- `qc_summary.csv`
- `strict_clone_table.csv`
- `vgene_group_table.csv`
- `clone_count_agreement.csv`
- `strict_clonal_diversity.csv`
- `relaxed_clonal_diversity.csv`
- `strict_vs_relaxed_diversity.csv`
- `sharing_table.csv`
- `candidate_table.csv`
- `claim_table.csv`
- `xenium_panel_roadmap.csv`
- `clone_cards.html`
- `tcr_claim_report.md`
