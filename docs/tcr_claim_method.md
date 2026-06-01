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
- `claim_table`: one row per auditable biological claim.
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
