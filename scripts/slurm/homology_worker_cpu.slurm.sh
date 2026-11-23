#!/bin/bash
# Nibi: FOUR bounded workers, not a 256/1536-task submission.
#SBATCH --account=def-hup-ab
#SBATCH --job-name=wm-shard-worker
#SBATCH --time=11:50:00
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=6G
#SBATCH --output=slurm-%A_%a.log
#SBATCH --error=slurm-%A_%a.log
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source scripts/config/paths.sh
: "$STAGE"
: "$SLURM_ARRAY_TASK_ID"
if [ "$STAGE" != graph_search ] && [ "$STAGE" != audit_search ]; then
  echo "Invalid worker stage $STAGE" >&2
  exit 2
fi
WORKERS=4
THREADS="$SLURM_CPUS_PER_TASK"
TOTAL=256
KIND=graph
if [ "$STAGE" = audit_search ]; then
  TOTAL=1536
  KIND=audit
fi
if (( SLURM_ARRAY_TASK_ID < 0 || SLURM_ARRAY_TASK_ID >= WORKERS )); then
  echo "Worker index exceeds $WORKERS" >&2
  exit 2
fi
PYTHON="$(command -v python3)"
MMSEQS="$(command -v mmseqs)"
export OMP_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS=1
mkdir -p "$TMP_ROOT/homology_v2"
sha256sum -c "$SOURCE_MANIFEST"
start_epoch="$(date +%s)"
completed_this_worker=0
echo "Worker=$SLURM_ARRAY_TASK_ID/$WORKERS kind=$KIND CPU=$THREADS"
echo "Each query shard searches the FULL target protein database"
for ((task=SLURM_ARRAY_TASK_ID; task<TOTAL; task+=WORKERS)); do
  elapsed=$(( $(date +%s) - start_epoch ))
  remaining=$(( 41100 - elapsed )) # stop 25 minutes before hard 11h50 limit
  if (( remaining <= 0 )); then
    echo "BUDGET_REACHED: stop taking new shards; rerun the same submit command"
    break
  fi
  # Individual long alignments are stopped before Slurm's hard limit;
  # no done marker will be written for the incomplete shard.
  if timeout --signal=TERM --kill-after=60s "$remaining"s \
    "$PYTHON" scripts/data_pipeline/sharded_mmseqs.py search \
      --kind "$KIND" --root "$SPLIT_ROOT" \
      --tmp-root "$TMP_ROOT/homology_v2" \
      --shards 256 --task "$task" --mmseqs "$MMSEQS" --threads "$THREADS"; then
    completed_this_worker=$((completed_this_worker + 1))
  else
    status=$?
    echo "SHARD_FAILED kind=$KIND shard=$task exit=$status" >&2
    exit "$status"
  fi
done
echo "BATCH_FINISHED kind=$KIND worker=$SLURM_ARRAY_TASK_ID visited=$completed_this_worker"
