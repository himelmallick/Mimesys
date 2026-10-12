#!/usr/bin/env bash
#SBATCH --ntasks=1

set -euo pipefail

GRID_FILE="${1:?Missing grid file argument}"
OUT_DIR="${2:?Missing output dir argument}"
TASK_OFFSET="${3:-0}"
DATA_PATH="${4:-$HOME/Mimesys_2.0/inst/extdata/StelzerDOS.RData}"
GLOBAL_TASK_ID=$((SLURM_ARRAY_TASK_ID + TASK_OFFSET))

cd "${HOME}/Mimesys_2.0"

Rscript scripts/slurm/run_stelzer_multiview_benchmark_task.R \
  --grid="$GRID_FILE" \
  --out_dir="$OUT_DIR" \
  --data_path="$DATA_PATH" \
  --task_id="$GLOBAL_TASK_ID"
