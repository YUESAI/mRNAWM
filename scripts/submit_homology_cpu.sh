#!/bin/bash
# Nibi CPU quota-safe launcher: submit ONE stage or one bounded worker batch.
# Run again after the submitted batch completes. Never queues future stages.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
source scripts/config/paths.sh
bash scripts/init_mrnawm_workspace.sh >/dev/null
command -v sbatch >/dev/null || { echo "ERROR: run on a Nibi login node (sbatch missing)" >&2; exit 2; }
command -v squeue >/dev/null || { echo "ERROR: squeue missing" >&2; exit 2; }
command -v python3 >/dev/null || { echo "ERROR: python3 missing" >&2; exit 2; }
command -v mmseqs >/dev/null || { echo "ERROR: mmseqs missing in active environment" >&2; exit 2; }

# Do not stack a new batch on top of any previous graph/audit tasks.
# The earlier all-at-once launcher may have left pending wm-homology jobs.
active="$(squeue -h -u "$USER" -o "%j" | grep -Ec '^(wm-homology|wm-shard-worker)$' || true)"
if (( active > 0 )); then
  echo "BLOCKED: $active existing mRNAWM Slurm job(s) still pending/running."
  echo 'Inspect: squeue -u "$USER" -o "%.18i %.22j %.10T %.30R"'
  echo "Do not resubmit or cancel automatically; inspect those jobs first."
  exit 3
fi

stage="$(python3 scripts/pipeline/next_stage.py --root "$SPLIT_ROOT")"
case "$stage" in
  DONE)
    echo "DONE: final audit passed; outputs: $SPLIT_ROOT"
    exit 0 ;;
  AUDIT_FAILED)
    echo "BLOCKED: final homology audit did not pass; investigate before training." >&2
    exit 4 ;;
esac

mkdir -p "$LOG_ROOT" "$MRNAWM_ROOT/00_source"
log="$MRNAWM_ROOT/00_source/submissions.tsv"
job=scripts/slurm/homology_sharded_cpu.slurm.sh
submitted=""
workers=""
if [ "$stage" = graph_search ] || [ "$stage" = audit_search ]; then
  job=scripts/slurm/homology_worker_cpu.slurm.sh
  # Request at most 4 tasks. If quota is lower, retry with 2, then 1.
  # The worker's modulo stride matches WORKER_COUNT.
  for n in 4 2 1; do
    last=$((n - 1))
    echo "Trying $n CPU workers, each 12 CPUs, under def-hup-ab ..."
    if submitted="$(sbatch --parsable --account=def-hup-ab --time=11:50:00 \
      --array="0-$last%$n" \
      --output="$LOG_ROOT/$stage-%A_%a.log" \
      --error="$LOG_ROOT/$stage-%A_%a.log" \
      --export="ALL,STAGE=$stage,WORKER_COUNT=$n" "$job")"; then
      workers="$n"
      break
    fi
    submitted=""
  done
else
  # One job with four CPU cores; no future dependency chain.
  if ! submitted="$(sbatch --parsable --account=def-hup-ab --time=11:50:00 \
      --cpus-per-task=4 \
      --output="$LOG_ROOT/$stage-%j.log" \
      --error="$LOG_ROOT/$stage-%j.log" \
      --export="ALL,STAGE=$stage" "$job")"; then
    echo "ERROR: Slurm rejected stage=$stage; no subsequent jobs were submitted." >&2
    exit 5
  fi
fi

if [ -z "$submitted" ]; then
  echo "ERROR: def-hup-ab submission limit still reached. Inspect queue before retry." >&2
  exit 5
fi
if [ ! -s "$log" ]; then
  printf 'utc_time\tstage\tjob_id\tworkers\tinput\toutput\tgit_sha\n' > "$log"
fi
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$(date -u +%FT%TZ)" "$stage" "$submitted" "$workers" \
  "$SOURCE_TSV" "$SPLIT_ROOT" "$(git rev-parse HEAD)" >> "$log"

echo "SUBMITTED: stage=$stage job=$submitted workers=$workers"
echo "INPUT: $SOURCE_TSV"
echo "OUTPUT: $SPLIT_ROOT"
echo "LOGS: $LOG_ROOT"
echo "After this stage's jobs exit, run the SAME command to advance/resume:"
echo "  bash scripts/submit_homology_cpu.sh"
squeue -u "$USER"
