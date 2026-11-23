#!/bin/bash
# mRNAWM homology split: Nibi CPU stages, following mRNA_pr's CPU Slurm layout.
# The old def-hup-ab_cpu/base settings are Dusky-specific; account is supplied
# with sbatch --account=... for the user's actual Nibi allocation.
#SBATCH --job-name=mrnawm-split
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=11:50:00
#SBATCH --output=slurm-mrnawm-split-%j.out
#SBATCH --error=slurm-mrnawm-split-%j.err

set -euo pipefail
: "$STAGE"
: "$DATA_ROOT"
: "$SLURM_SUBMIT_DIR"
cd "$SLURM_SUBMIT_DIR"

INPUT="$DATA_ROOT/normalized/all_refseq_mammalian_viral_unique.tsv"
WORK="$DATA_ROOT/splits/homology_v1_seed3407"
THREADS="$SLURM_CPUS_PER_TASK"
PYTHON_BIN="$(command -v python3)"
export OMP_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS=1
mkdir -p "$WORK"

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
    mkdir -p "$WORK/mmseqs_tmp"
    mmseqs easy-search \
      "$WORK/unique_proteins.fasta" "$WORK/unique_proteins.fasta" \
      "$WORK/all_vs_all.tsv" "$WORK/mmseqs_tmp" \
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
    mkdir -p "$WORK/audit_tmp"
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
          "$WORK/audit_tmp/$query-$target" \
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
