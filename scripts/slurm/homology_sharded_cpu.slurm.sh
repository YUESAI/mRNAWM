#!/bin/bash
# Nibi CPU-only resumable sharded search; submit via scripts/submit_homology_cpu.sh
#SBATCH --account=def-hup-ab
#SBATCH --job-name=wm-homology
#SBATCH --time=11:50:00
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=6G
#SBATCH --output=slurm-%A_%a.log
#SBATCH --error=slurm-%A_%a.log
set -euo pipefail

cd "$SLURM_SUBMIT_DIR"
source scripts/config/paths.sh
: "$STAGE"
PYTHON="$(command -v python3)"
MMSEQS="$(command -v mmseqs)"
THREADS="$SLURM_CPUS_PER_TASK"
SHARDS=256
ROOT="$SPLIT_ROOT"
TMP="$TMP_ROOT/homology_v2"
mkdir -p "$ROOT" "$TMP" "$LOG_ROOT" "$MRNAWM_ROOT/00_source/stage_records"
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS=1
echo "STAGE=$STAGE JOB=$SLURM_JOB_ID CPU=$THREADS DATE=$(date -u +%FT%TZ)"
echo "INPUT=$SOURCE_TSV ROOT=$ROOT TMP=$TMP"
sha256sum -c "$SOURCE_MANIFEST"

records="$MRNAWM_ROOT/00_source/stage_records"
case "$STAGE" in
  protein)
    "$PYTHON" scripts/data_pipeline/homology_split.py prepare \
      --input "$SOURCE_TSV" --output-dir "$ROOT"
    "$PYTHON" scripts/pipeline/record_stage.py --stage protein \
      --record "$records/protein_v2.json" \
      --input "normalized=$SOURCE_TSV" \
      --output "proteins=$ROOT/unique_proteins.fasta" \
      --output "prepare_manifest=$ROOT/prepare.json"
    ;;
  graph_prepare)
    "$PYTHON" scripts/data_pipeline/sharded_mmseqs.py prepare --kind graph \
      --root "$ROOT" --shards "$SHARDS" --mmseqs "$MMSEQS" --threads "$THREADS"
    "$PYTHON" scripts/pipeline/record_stage.py --stage graph_prepare \
      --record "$records/graph_prepare_v2.json" \
      --input "proteins=$ROOT/unique_proteins.fasta" \
      --output "graph_ready=$ROOT/graph_shards_v1/READY.json"
    ;;
  graph_search)
    "$PYTHON" scripts/data_pipeline/sharded_mmseqs.py search --kind graph \
      --root "$ROOT" --tmp-root "$TMP" --shards "$SHARDS" \
      --task "$SLURM_ARRAY_TASK_ID" --mmseqs "$MMSEQS" --threads "$THREADS"
    ;;
  graph_verify)
    "$PYTHON" scripts/data_pipeline/sharded_mmseqs.py verify --kind graph \
      --root "$ROOT" --shards "$SHARDS"
    "$PYTHON" scripts/pipeline/record_stage.py --stage graph_verify \
      --record "$records/graph_verify_v2.json" \
      --input "graph_ready=$ROOT/graph_shards_v1/READY.json" \
      --output "complete=$ROOT/graph_shards_v1/COMPLETE.json" \
      --output "hits_list=$ROOT/graph_shards_v1/hits.list"
    ;;
  split)
    "$PYTHON" scripts/data_pipeline/homology_split.py split \
      --input "$SOURCE_TSV" --prepare-manifest "$ROOT/prepare.json" \
      --hits-list "$ROOT/graph_shards_v1/hits.list" \
      --output-dir "$ROOT" --seed 3407
    "$PYTHON" scripts/pipeline/record_stage.py --stage split \
      --record "$records/split_v2.json" \
      --input "normalized=$SOURCE_TSV" \
      --input "graph_complete=$ROOT/graph_shards_v1/COMPLETE.json" \
      --output "train=$ROOT/train.tsv" --output "val=$ROOT/val.tsv" \
      --output "test=$ROOT/test.tsv" --output "manifest=$ROOT/split_manifest.json" \
      --output "cluster_table=$ROOT/clusters.tsv"
    ;;
  audit_prepare)
    "$PYTHON" scripts/data_pipeline/sharded_mmseqs.py prepare --kind audit \
      --root "$ROOT" --shards "$SHARDS" --mmseqs "$MMSEQS" --threads "$THREADS"
    "$PYTHON" scripts/pipeline/record_stage.py --stage audit_prepare \
      --record "$records/audit_prepare_v2.json" \
      --input "split_manifest=$ROOT/split_manifest.json" \
      --output "audit_ready=$ROOT/audit_shards_v1/READY.json"
    ;;
  audit_search)
    : "$AUDIT_PAIR_INDEX"
    if [ "$AUDIT_PAIR_INDEX" -lt 0 ] || [ "$AUDIT_PAIR_INDEX" -gt 5 ]; then
      echo "Invalid audit pair index: $AUDIT_PAIR_INDEX" >&2
      exit 2
    fi
    task=$(( AUDIT_PAIR_INDEX * SHARDS + SLURM_ARRAY_TASK_ID ))
    "$PYTHON" scripts/data_pipeline/sharded_mmseqs.py search --kind audit \
      --root "$ROOT" --tmp-root "$TMP" --shards "$SHARDS" \
      --task "$task" --mmseqs "$MMSEQS" --threads "$THREADS"
    ;;
  audit_verify)
    "$PYTHON" scripts/data_pipeline/sharded_mmseqs.py verify --kind audit \
      --root "$ROOT" --shards "$SHARDS"
    "$PYTHON" scripts/pipeline/record_stage.py --stage audit_verify \
      --record "$records/audit_verify_v2.json" \
      --input "audit_ready=$ROOT/audit_shards_v1/READY.json" \
      --output "complete=$ROOT/audit_shards_v1/COMPLETE.json" \
      --output "hits_list=$ROOT/audit_shards_v1/hits.list"
    ;;
  audit)
    "$PYTHON" scripts/data_pipeline/homology_split.py audit \
      --train "$ROOT/train.tsv" --val "$ROOT/val.tsv" --test "$ROOT/test.tsv" \
      --hits-list "$ROOT/audit_shards_v1/hits.list" \
      --report "$ROOT/leakage_audit.json"
    "$PYTHON" scripts/pipeline/record_stage.py --stage audit \
      --record "$records/audit_v2.json" \
      --input "train=$ROOT/train.tsv" --input "val=$ROOT/val.tsv" \
      --input "test=$ROOT/test.tsv" \
      --input "audit_search_complete=$ROOT/audit_shards_v1/COMPLETE.json" \
      --output "leakage_audit=$ROOT/leakage_audit.json"
    ;;
  *)
    echo "Unsupported stage $STAGE" >&2
    exit 2
    ;;
esac
echo "SUCCESS $STAGE $(date -u +%FT%TZ)"
