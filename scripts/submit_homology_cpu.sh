#!/bin/bash
# Nibi: resumable 256-shard global homology graph + 1536-shard audit.
# Start: bash scripts/submit_homology_cpu.sh
# Resume: bash scripts/submit_homology_cpu.sh resume-graph
# Resume audit only: bash scripts/submit_homology_cpu.sh resume-audit
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
source scripts/config/paths.sh
bash scripts/init_mrnawm_workspace.sh
command -v sbatch >/dev/null || { echo "ERROR: sbatch unavailable (Nibi only)" >&2; exit 2; }
command -v python3 >/dev/null || { echo "ERROR: python3 unavailable" >&2; exit 2; }
command -v mmseqs >/dev/null || { echo "ERROR: mmseqs unavailable in environment" >&2; exit 2; }
mmseqs version
mode="${1:-all}"
job=scripts/slurm/homology_sharded_cpu.slurm.sh
mkdir -p "$LOG_ROOT" "$SPLIT_ROOT"
submission="$MRNAWM_ROOT/00_source/last_submission.txt"

simple() {
  local stage="$1" dep="$2"
  if [ -n "$dep" ]; then
    sbatch --parsable --account=def-hup-ab --cpus-per-task=4 \
      --dependency="afterok:$dep" \
      --output="$LOG_ROOT/$stage-%j.log" --error="$LOG_ROOT/$stage-%j.log" \
      --export="ALL,STAGE=$stage" "$job"
  else
    sbatch --parsable --account=def-hup-ab --cpus-per-task=4 \
      --output="$LOG_ROOT/$stage-%j.log" --error="$LOG_ROOT/$stage-%j.log" \
      --export="ALL,STAGE=$stage" "$job"
  fi
}
parallel() {
  local stage="$1" dep="$2" range="$3"
  if [ -n "$dep" ]; then
    sbatch --parsable --account=def-hup-ab --cpus-per-task=12 \
      --dependency="afterok:$dep" --array="$range" \
      --output="$LOG_ROOT/$stage-%A_%a.log" --error="$LOG_ROOT/$stage-%A_%a.log" \
      --export="ALL,STAGE=$stage" "$job"
  else
    sbatch --parsable --account=def-hup-ab --cpus-per-task=12 \
      --array="$range" \
      --output="$LOG_ROOT/$stage-%A_%a.log" --error="$LOG_ROOT/$stage-%A_%a.log" \
      --export="ALL,STAGE=$stage" "$job"
  fi
}
echo "NIBI INPUT: $SOURCE_TSV"
echo "NEW OUTPUT: $SPLIT_ROOT"
echo "LOGS: $LOG_ROOT"
echo "MODE: $mode"
case "$mode" in
  all)
    one=$(simple protein "")
    two=$(simple graph_prepare "$one")
    three=$(parallel graph_search "$two" "0-255%6")
    four=$(simple graph_verify "$three")
    five=$(simple split "$four")
    six=$(simple audit_prepare "$five")
    seven=$(parallel audit_search "$six" "0-1535%6")
    eight=$(simple audit_verify "$seven")
    nine=$(simple audit "$eight")
    ;;
  resume-graph)
    test -f "$SPLIT_ROOT/graph_shards_v1/READY.json" || {
      echo "ERROR: graph not prepared; run 'all' first" >&2; exit 2;
    }
    three=$(parallel graph_search "" "0-255%6")
    four=$(simple graph_verify "$three")
    five=$(simple split "$four")
    six=$(simple audit_prepare "$five")
    seven=$(parallel audit_search "$six" "0-1535%6")
    eight=$(simple audit_verify "$seven")
    nine=$(simple audit "$eight")
    ;;
  resume-audit)
    test -f "$SPLIT_ROOT/audit_shards_v1/READY.json" || {
      echo "ERROR: audit not prepared; cannot resume" >&2; exit 2;
    }
    seven=$(parallel audit_search "" "0-1535%6")
    eight=$(simple audit_verify "$seven")
    nine=$(simple audit "$eight")
    ;;
  *)
    echo "Usage: bash scripts/submit_homology_cpu.sh [all|resume-graph|resume-audit]" >&2
    exit 2
    ;;
esac
{
  printf 'mode=%s\n' "$mode"
  printf 'input=%s\noutput=%s\n' "$SOURCE_TSV" "$SPLIT_ROOT"
  printf 'git_commit=%s\n' "$(git rev-parse HEAD)"
  if [ "$mode" = all ]; then printf 'protein=%s\ngraph_prepare=%s\n' "$one" "$two"; fi
  if [ "$mode" != resume-audit ]; then
    printf 'graph_search=%s\ngraph_verify=%s\nsplit=%s\naudit_prepare=%s\n' \
      "$three" "$four" "$five" "$six"
  fi
  printf 'audit_search=%s\naudit_verify=%s\naudit=%s\n' "$seven" "$eight" "$nine"
} | tee "$submission"
squeue -u "$USER"
