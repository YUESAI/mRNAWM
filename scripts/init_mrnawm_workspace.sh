#!/bin/bash
# Idempotent Nibi workspace setup. Does NOT download, copy or modify source data.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config/paths.sh"

mkdir -p "$MRNAWM_ROOT/00_source" "$SPLIT_ROOT" "$MODEL_DATA_ROOT" \
  "$FEATURE_CACHE_ROOT" "$CHECKPOINT_ROOT" "$BASELINES_ROOT" \
  "$EVALUATIONS_ROOT" "$LOG_ROOT" "$TMP_ROOT" "$EXTERNAL_ROOT"

if [[ ! -s "$SOURCE_ORIGINAL" ]]; then
  echo "ERROR: required existing normalized RefSeq input is missing:" >&2
  echo "  $SOURCE_ORIGINAL" >&2
  echo "Available files under legacy normalized/:" >&2
  ls -lah "$LEGACY_DATA_ROOT/normalized" >&2 || true
  echo "Nothing was changed in the original data directory." >&2
  exit 2
fi

if [[ -e "$SOURCE_TSV" || -L "$SOURCE_TSV" ]]; then
  if [[ ! -L "$SOURCE_TSV" ]] || \
     [[ "$(readlink -f "$SOURCE_TSV")" != "$(readlink -f "$SOURCE_ORIGINAL")" ]]; then
    echo "ERROR: $SOURCE_TSV already exists but is not the expected symlink" >&2
    exit 3
  fi
else
  ln -s "$SOURCE_ORIGINAL" "$SOURCE_TSV"
fi

current_hash="$(sha256sum "$SOURCE_ORIGINAL" | cut -d' ' -f1)"
if [[ -f "$SOURCE_MANIFEST" ]]; then
  previous_hash="$(cut -d' ' -f1 < "$SOURCE_MANIFEST")"
  if [[ "$previous_hash" != "$current_hash" ]]; then
    echo "ERROR: original source changed since workspace was initialized." >&2
    echo "Old hash: $previous_hash; current hash: $current_hash" >&2
    echo "Use a NEW versioned experiment directory; do not silently overwrite." >&2
    exit 4
  fi
else
  printf '%s  %s\n' "$current_hash" "$SOURCE_ORIGINAL" > "$SOURCE_MANIFEST"
fi

echo "Workspace: $MRNAWM_ROOT"
echo "Read-only source: $SOURCE_ORIGINAL"
echo "Source SHA-256: $current_hash"
echo "Splits: $SPLIT_ROOT"
echo "Checkpoints: $CHECKPOINT_ROOT"
echo "Baselines: $BASELINES_ROOT"
echo "Evaluations: $EVALUATIONS_ROOT"
