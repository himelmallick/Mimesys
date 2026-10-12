#!/usr/bin/env bash
#SBATCH --ntasks=1

set -euo pipefail

GRID_FILE="${1:?Missing grid file argument}"
BASE_DIR="${2:?Missing base dir argument}"
OUT_DIR="${3:?Missing output dir argument}"
TASK_OFFSET="${4:-0}"
GLOBAL_TASK_ID=$((SLURM_ARRAY_TASK_ID + TASK_OFFSET))

cd "${HOME}/Mimesys_2.0"

Rscript scripts/slurm/run_tcga_validation_task.R \
  --grid="$GRID_FILE" \
  --base_dir="$BASE_DIR" \
  --out_dir="$OUT_DIR" \
  --task_id="$GLOBAL_TASK_ID" \
  --figure_caps="${FIGURE_CAPS:-}"
