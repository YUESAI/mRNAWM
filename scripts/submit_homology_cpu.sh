#!/bin/bash
# One-command CPU-only, fully ordered Nibi submission. Run from the repo root:
# bash scripts/submit_homology_cpu.sh
set -euo pipefail
repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"
source scripts/config/paths.sh
bash scripts/init_mrnawm_workspace.sh

command -v sbatch >/dev/null || { echo "ERROR: sbatch unavailable. Run this on Nibi." >&2; exit 2; }
command -v mmseqs >/dev/null || {
  echo "ERROR: mmseqs unavailable on PATH; load the installed MMseqs2 environment first." >&2
  exit 2
}
command -v python3 >/dev/null || { echo "ERROR: python3 unavailable" >&2; exit 2; }
script=scripts/slurm/submit_homology_split_cpu.slurm.sh
log_args=(--output="$LOG_ROOT/homology-%x-%j.log" --error="$LOG_ROOT/homology-%x-%j.log")

# Source preparation and Python splitting are not 32-core operations.
prepare_id="$(sbatch --parsable "${log_args[@]}" --cpus-per-task=4 --export=ALL,STAGE=prepare "$script")"
search_id="$(sbatch --parsable "${log_args[@]}" --dependency=afterok:"$prepare_id" \
  --export=ALL,STAGE=search "$script")"
split_id="$(sbatch --parsable "${log_args[@]}" --cpus-per-task=4 \
  --dependency=afterok:"$search_id" --export=ALL,STAGE=split "$script")"
audit_id="$(sbatch --parsable "${log_args[@]}" --dependency=afterok:"$split_id" \
  --export=ALL,STAGE=audit "$script")"

printf 'Nibi CPU jobs\nprepare: %s\nsearch: %s\nsplit: %s\naudit: %s\n' \
  "$prepare_id" "$search_id" "$split_id" "$audit_id"
printf 'Logs: %s\nResults: %s\n' "$LOG_ROOT" "$SPLIT_ROOT"
squeue -u "$USER"
