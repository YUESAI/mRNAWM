#!/usr/bin/env bash
#SBATCH --job-name=mrnawm-homology
#SBATCH --account=def-hup-ab_cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=slurm-mrnawm-homology-%j.log
#SBATCH --error=slurm-mrnawm-homology-%j.log
# CPU only. Set STAGE via sbatch --export.
set -euo pipefail
: "$STAGE"
: "$DATA_ROOT"
INPUT="$DATA_ROOT/normalized/all_refseq_mammalian_viral_unique.tsv"
WORK="$DATA_ROOT/splits/homology_v1_seed3407"
THREADS=16
mkdir -p "$WORK"
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS="$THREADS"
PYTHON=python3
echo "STAGE=$STAGE HOST=$(hostname) WORK=$WORK"
case "$STAGE" in
  prepare)
    test -s "$INPUT"
    "$PYTHON" scripts/data_pipeline/homology_split.py prepare \
      --input "$INPUT" --output-dir "$WORK"
    ;;
  search)
    command -v mmseqs
    test -s "$WORK/unique_proteins.fasta"
    mmseqs version > "$WORK/mmseqs_version.txt"
    mkdir -p "$WORK/tmp"
    mmseqs easy-search \
      "$WORK/unique_proteins.fasta" "$WORK/unique_proteins.fasta" \
      "$WORK/all_vs_all.tsv" "$WORK/tmp" \
      --format-output "query,target,pident,alnlen,qlen,tlen,qcov,tcov" \
      --min-seq-id 0.25 -c 0.0 -s 7.5 --max-seqs 100000 \
      --threads "$THREADS"
    test -s "$WORK/all_vs_all.tsv"
    ;;
  split)
    test -s "$WORK/prepare.json"
    test -s "$WORK/all_vs_all.tsv"
    "$PYTHON" scripts/data_pipeline/homology_split.py split \
      --input "$INPUT" --prepare-manifest "$WORK/prepare.json" \
      --hits "$WORK/all_vs_all.tsv" --output-dir "$WORK" \
      --seed 3407 --train-fraction 0.8 --val-fraction 0.1
    ;;
  audit)
    command -v mmseqs
    test -s "$WORK/split_manifest.json"
    mkdir -p "$WORK/audit_tmp"
    : > "$WORK/cross_split_hits.tsv"
    for query in train val test; do
      for target in train val test; do
        [[ "$query" == "$target" ]] && continue
        test -s "$WORK/$query.proteins.fasta"
        test -s "$WORK/$target.proteins.fasta"
        mmseqs easy-search \
          "$WORK/$query.proteins.fasta" "$WORK/$target.proteins.fasta" \
          "$WORK/audit_${query}_${target}.tsv" "$WORK/audit_tmp/${query}_${target}" \
          --format-output "query,target,pident,alnlen,qlen,tlen,qcov,tcov" \
          --min-seq-id 0.25 -c 0.0 -s 7.5 --max-seqs 100000 \
          --threads "$THREADS"
        cat "$WORK/audit_${query}_${target}.tsv" >> "$WORK/cross_split_hits.tsv"
      done
    done
    if [[ ! -s "$WORK/cross_split_hits.tsv" ]]; then
      printf 'query\ttarget\tpident\talnlen\tqlen\ttlen\tqcov\ttcov\n' > "$WORK/cross_split_hits.tsv"
    fi
    "$PYTHON" scripts/data_pipeline/homology_split.py audit \
      --train "$WORK/train.tsv" --val "$WORK/val.tsv" --test "$WORK/test.tsv" \
      --hits "$WORK/cross_split_hits.tsv" --report "$WORK/leakage_audit.json"
    ;;
  *)
    echo "Invalid STAGE=$STAGE" >&2
    exit 2
    ;;
esac
echo "Finished $STAGE"
