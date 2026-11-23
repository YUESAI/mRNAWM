# Nibi CPU Slurm — current

**Single source of truth: [root README](../README.md).** All path constants live in [`scripts/config/paths.sh`](../scripts/config/paths.sh). Your repository checkout is already the working directory; no placeholder `cd` or `DATA_ROOT` assignment is necessary.

Run from the Nibi Git repository root:

```bash
git pull origin main
bash scripts/init_mrnawm_workspace.sh
bash scripts/submit_homology_cpu.sh
```

The first script creates `/scratch/ylong66/data/mRNA_pr_wm/mRNAWM/` and its staged directories, using the previous normalized TSV **read-only**. The submit script queues four ordered CPU tasks (prepare, 32-core MMseqs2 search, split, 32-core independent audit), using Nibi account `rrg-hup`. It never requests a GPU. Logs and outputs are isolated under the new mRNAWM subfolder.

Verify outputs:

```bash
source scripts/config/paths.sh
cat "$SPLIT_ROOT/split_manifest.json"
cat "$SPLIT_ROOT/leakage_audit.json"
ls "$LOG_ROOT"
```

Full-data execution on Nibi has **not** been performed from this conversation. No stage is marked complete solely because code is committed.
