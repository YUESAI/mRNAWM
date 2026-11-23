#!/bin/bash
# mRNAWM homology split: Nibi CPU stages, following mRNA_pr's CPU Slurm layout.
# Nibi CPU account confirmed by sbatch: def-hup-ab. No GPU requested.
#SBATCH --account=def-hup-ab
#SBATCH --job-name=mrnawm-split
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=6G
#SBATCH --time=12:00:00
#SBATCH --output=slurm-mrnawm-split-%j.log
#SBATCH --error=slurm-mrnawm-split-%j.log

set -euo pipefail
: "${STAGE:?set STAGE to prepare/search/split/audit}"
cd "$SLURM_SUBMIT_DIR"
source scripts/config/paths.sh
test -L "$SOURCE_TSV" || { echo "Run bash scripts/init_mrnawm_workspace.sh first" >&2; exit 2; }
sha256sum -c "$SOURCE_MANIFEST" || { echo "Source changed. Aborting." >&2; exit 2; }
INPUT="$SOURCE_TSV"
WORK="$SPLIT_ROOT"
THREADS="$SLURM_CPUS_PER_TASK"
PYTHON_BIN="$(command -v python3)"
export OMP_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS=1
mkdir -p "$WORK" "$LOG_ROOT" "$TMP_ROOT"

echo "stage=$STAGE job=$SLURM_JOB_ID host=$(hostname) cpus=$THREADS"
echo "repo=$PWD input=$INPUT output=$WORK"
case "$STAGE" in
  prepare)
    test -s "$INPUT"
    "$PYTHON_BIN" scripts/data_pipeline/homology_split.py prepare \
      --input "$INPUT" --output-dir "$WORK"
    ;;
  search)
    command -v mmseqs
    test -s "$WORK/unique_proteins.fasta"
    mmseqs version > "$WORK/mmseqs_version.txt"
    mkdir -p "$TMP_ROOT/mmseqs_all_vs_all"
    mmseqs easy-search \
      "$WORK/unique_proteins.fasta" "$WORK/unique_proteins.fasta" \
      "$WORK/all_vs_all.tsv" "$TMP_ROOT/mmseqs_all_vs_all" \
      --format-output "query,target,pident,alnlen,qlen,tlen,qcov,tcov" \
      --min-seq-id 0.25 -c 0.0 -s 7.5 --max-seqs 100000 \
      --threads "$THREADS"
    test -s "$WORK/all_vs_all.tsv"
    ;;
  split)
    test -s "$WORK/prepare.json"
    test -s "$WORK/all_vs_all.tsv"
    "$PYTHON_BIN" scripts/data_pipeline/homology_split.py split \
      --input "$INPUT" \
      --prepare-manifest "$WORK/prepare.json" \
      --hits "$WORK/all_vs_all.tsv" \
      --output-dir "$WORK" --seed 3407 \
      --train-fraction 0.8 --val-fraction 0.1
    ;;
  audit)
    command -v mmseqs
    test -s "$WORK/split_manifest.json"
    mkdir -p "$TMP_ROOT/mmseqs_cross_split"
    : > "$WORK/cross_split_hits.tsv"
    for query in train val test; do
      for target in train val test; do
        if [ "$query" = "$target" ]; then
          continue
        fi
        test -s "$WORK/$query.proteins.fasta"
        test -s "$WORK/$target.proteins.fasta"
        mmseqs easy-search \
          "$WORK/$query.proteins.fasta" "$WORK/$target.proteins.fasta" \
          "$WORK/audit-$query-$target.tsv" \
          "$TMP_ROOT/mmseqs_cross_split/$query-$target" \
          --format-output "query,target,pident,alnlen,qlen,tlen,qcov,tcov" \
          --min-seq-id 0.25 -c 0.0 -s 7.5 --max-seqs 100000 \
          --threads "$THREADS"
        cat "$WORK/audit-$query-$target.tsv" >> "$WORK/cross_split_hits.tsv"
      done
    done
    if [ ! -s "$WORK/cross_split_hits.tsv" ]; then
      printf 'query\ttarget\tpident\talnlen\tqlen\ttlen\tqcov\ttcov\n' \
        > "$WORK/cross_split_hits.tsv"
    fi
    "$PYTHON_BIN" scripts/data_pipeline/homology_split.py audit \
      --train "$WORK/train.tsv" \
      --val "$WORK/val.tsv" \
      --test "$WORK/test.tsv" \
      --hits "$WORK/cross_split_hits.tsv" \
      --report "$WORK/leakage_audit.json"
    ;;
  *)
    echo "Unknown STAGE: $STAGE (expected prepare/search/split/audit)" >&2
    exit 2
    ;;
esac
echo "Completed $STAGE at $(date)"
