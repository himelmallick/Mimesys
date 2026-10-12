#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${BASE_DIR:-$HOME/Mimesys_validation}"
PROJECT_DIR="${PROJECT_DIR:-$HOME/Mimesys_2.0}"
OUT_DIR="${OUT_DIR:-$BASE_DIR/stelzer_multiview_cooperative_benchmark}"
GRID_FILE="${GRID_FILE:-$OUT_DIR/stelzer_multiview_benchmark_grid.csv}"
DATA_PATH="${DATA_PATH:-$PROJECT_DIR/inst/extdata/StelzerDOS.RData}"
BATCH_SIZE="${BATCH_SIZE:-1000}"
SEQUENTIAL="${SEQUENTIAL:-1}"
WAIT_FOR_BATCH="${WAIT_FOR_BATCH:-0}"
START_ROW="${START_ROW:-1}"
END_ROW="${END_ROW:-}"
WAIT_POLL_SECONDS="${WAIT_POLL_SECONDS:-120}"
SBATCH_EXTRA="${SBATCH_EXTRA:-}"

TIME="${TIME:-06:00:00}"
MEM="${MEM:-20G}"
CPUS="${CPUS:-1}"

mkdir -p "$OUT_DIR/logs"

if [[ ! -f "$GRID_FILE" ]]; then
  echo "Grid file not found: $GRID_FILE"
  exit 1
fi

N_ROWS=$(($(wc -l < "$GRID_FILE") - 1))
if [[ "$N_ROWS" -le 0 ]]; then
  echo "Grid has no rows: $GRID_FILE"
  exit 1
fi

if [[ -z "$END_ROW" ]]; then
  END_ROW="$N_ROWS"
fi

if [[ "$START_ROW" -lt 1 || "$START_ROW" -gt "$N_ROWS" ]]; then
  echo "START_ROW is outside grid range 1-$N_ROWS: $START_ROW"
  exit 1
fi

if [[ "$END_ROW" -lt "$START_ROW" || "$END_ROW" -gt "$N_ROWS" ]]; then
  echo "END_ROW is outside grid range $START_ROW-$N_ROWS: $END_ROW"
  exit 1
fi

EXTRA_ARGS=()
if [[ -n "$SBATCH_EXTRA" ]]; then
  # shellcheck disable=SC2206
  EXTRA_ARGS=($SBATCH_EXTRA)
fi

submit_batch() {
  local batch="$1"
  local local_end="$2"
  local task_offset="$3"
  shift 3
  local dependency_args=("$@")

  sbatch --parsable \
    "${EXTRA_ARGS[@]}" \
    "${dependency_args[@]}" \
    --job-name="stelzer_mv_${batch}" \
    --output="$OUT_DIR/logs/%x_%A_%a.out" \
    --error="$OUT_DIR/logs/%x_%A_%a.err" \
    --array="1-${local_end}" \
    --time="$TIME" \
    --mem="$MEM" \
    --cpus-per-task="$CPUS" \
    "$PROJECT_DIR/scripts/slurm/stelzer_multiview_benchmark_array.sbatch" "$GRID_FILE" "$OUT_DIR" "$task_offset" "$DATA_PATH"
}

wait_for_job_to_leave_queue() {
  local job_id="$1"
  echo "Waiting for Slurm job $job_id to finish before submitting next batch..."
  while squeue -j "$job_id" -h >/dev/null 2>&1 && [[ -n "$(squeue -j "$job_id" -h 2>/dev/null)" ]]; do
    sleep "$WAIT_POLL_SECONDS"
  done
  echo "Job $job_id is no longer in the queue."
}

START="$START_ROW"
BATCH=1
PREV_JOB_ID=""

while [[ "$START" -le "$END_ROW" ]]; do
  END=$((START + BATCH_SIZE - 1))
  if [[ "$END" -gt "$END_ROW" ]]; then
    END="$END_ROW"
  fi
  LOCAL_END=$((END - START + 1))
  TASK_OFFSET=$((START - 1))

  DEPENDENCY_ARGS=()
  if [[ "$WAIT_FOR_BATCH" != "1" && "$SEQUENTIAL" == "1" && -n "$PREV_JOB_ID" ]]; then
    DEPENDENCY_ARGS=(--dependency="afterok:${PREV_JOB_ID}")
  fi

  echo "Submitting batch $BATCH: grid rows $START-$END as local array 1-$LOCAL_END"
  JOB_ID=$(submit_batch "$BATCH" "$LOCAL_END" "$TASK_OFFSET" "${DEPENDENCY_ARGS[@]}")
  echo "Submitted Slurm job id: $JOB_ID"
  PREV_JOB_ID="$JOB_ID"

  if [[ "$WAIT_FOR_BATCH" == "1" ]]; then
    wait_for_job_to_leave_queue "$JOB_ID"
  fi

  START=$((END + 1))
  BATCH=$((BATCH + 1))
done
