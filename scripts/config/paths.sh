#!/bin/bash
# Single canonical filesystem layout for all Nibi mRNAWM runs.
# Source inputs in the old data root are READ ONLY; no new result is written there.
export LEGACY_DATA_ROOT=/scratch/ylong66/data/mRNA_pr_wm
export MRNAWM_ROOT=/scratch/ylong66/data/mRNA_pr_wm/mRNAWM
export SOURCE_ORIGINAL="$LEGACY_DATA_ROOT/normalized/all_refseq_mammalian_viral_unique.tsv"
export SOURCE_TSV="$MRNAWM_ROOT/00_source/normalized.tsv"
export SOURCE_MANIFEST="$MRNAWM_ROOT/00_source/source.sha256"
export SPLIT_ROOT="$MRNAWM_ROOT/01_split/homology_sharded_v2_seed3407"
export MODEL_DATA_ROOT="$MRNAWM_ROOT/02_model_data"
export FEATURE_CACHE_ROOT="$MRNAWM_ROOT/03_feature_cache"
export CHECKPOINT_ROOT="$MRNAWM_ROOT/04_checkpoints"
export BASELINES_ROOT="$MRNAWM_ROOT/05_baselines"
export EVALUATIONS_ROOT="$MRNAWM_ROOT/06_evaluations"
export LOG_ROOT="$MRNAWM_ROOT/logs"
export TMP_ROOT="$MRNAWM_ROOT/08_tmp"
export EXTERNAL_ROOT="$MRNAWM_ROOT/09_external_models"
# Old runner interfaces may use DATA_ROOT; always point it to NEW experiment outputs.
export DATA_ROOT="$MRNAWM_ROOT"
export SEED=3407
