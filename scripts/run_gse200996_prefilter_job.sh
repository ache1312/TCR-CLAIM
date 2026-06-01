#!/usr/bin/env bash
set -euo pipefail

RUN_ROOT="${1:?Usage: run_gse200996_prefilter_job.sh RUN_ROOT}"
REPO_ROOT="/mnt/data2/paper_clonotipos"
RESULTS_ROOT="/mnt/data1/Andres/clonotype_definition_benchmark/results"
GSE200996_INPUT="${RESULTS_ROOT}/gse200996/cell_metadata_with_tcr.csv"
FILTERED_DIR="${RUN_ROOT}/gse200996_primary"
FILTERED_INPUT="${FILTERED_DIR}/cell_metadata_with_tcr.csv"

cd "${REPO_ROOT}"
mkdir -p "${RUN_ROOT}" "${FILTERED_DIR}"

echo "[$(date -Is)] Starting TCR-CLAIM long job"
echo "Run root: ${RUN_ROOT}"
echo "Input: ${GSE200996_INPUT}"

echo "[$(date -Is)] Auditing all benchmark datasets"
python scripts/audit_tcr_claim_datasets.py \
  --results-root "${RESULTS_ROOT}" \
  --out "${RUN_ROOT}/dataset_audit" \
  --chunksize 500000

echo "[$(date -Is)] Prefiltering gse200996 to paired CD4/CD8 rows"
python scripts/prefilter_tcr_claim_dataset.py \
  --input "${GSE200996_INPUT}" \
  --out "${FILTERED_INPUT}" \
  --chunksize 250000 \
  --cell-classes CD4,CD8

echo "[$(date -Is)] Running TCR-CLAIM on filtered gse200996"
python scripts/run_tcr_claim_tables.py \
  --input "${FILTERED_INPUT}" \
  --out "${FILTERED_DIR}/tcr_claim_outputs"

echo "[$(date -Is)] Finished TCR-CLAIM long job"
