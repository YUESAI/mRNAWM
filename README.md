# mRNAWM — full-length mRNA world model

**Active repository:** `YUESAI/mRNAWM`, **main only**. Run everything on **Nibi** (CPU Slurm account `def-hup-ab`), from the directory where you have already cloned/opened this Git repo. `mRNA_pr` is **read-only**. **Seed = 3407** across split, model, baselines and evaluation.

## Canonical Nibi data paths (never guess a path)

**Existing legacy input (read-only):** `/scratch/ylong66/data/mRNA_pr_wm/`.

**All new experiment artifacts:** `/scratch/ylong66/data/mRNA_pr_wm/mRNAWM/`.

| Stage | Exact output location | Status |
|---|---|---|
| 00 — Input provenance | `mRNAWM/00_source/` (symlink and SHA-256 pointing to legacy `normalized/`) | Initialization script ready; not run |
| 01 — Homology split | `mRNAWM/01_split/homology_v1_seed3407/` | Code ready; Nibi run pending |
| 02 — Model-ready data | `mRNAWM/02_model_data/` | Migration pending |
| 03 — CodonTransformer features | `mRNAWM/03_feature_cache/` | Migration pending |
| 04 — mRNAWM checkpoints | `mRNAWM/04_checkpoints/` | Training code migration pending |
| 05 — Baseline predictions | `mRNAWM/05_baselines/` | Migration pending |
| 06 — Model and baseline evaluations | `mRNAWM/06_evaluations/` | Migration pending |
| 07 — Slurm/task logs | `mRNAWM/07_logs/` | Script ready |
| 08 — Temporary files | `mRNAWM/08_tmp/` | Script ready |
| 09 — External baseline weights | `mRNAWM/09_external_models/` | Migration pending |

The absolute path of every stage is defined **once** in [`scripts/config/paths.sh`](scripts/config/paths.sh). New scripts must load it. Never write newly generated splits, feature caches, baselines, or checkpoints into the original legacy directories.

## Run from the Git repo directory

### Stage 00 — initialize folders; reuse original normalized data

```bash
git pull origin main
bash scripts/init_mrnawm_workspace.sh
```

This creates the stage folders and a **symlink** from `00_source/normalized.tsv` to the original:

```text
/scratch/ylong66/data/mRNA_pr_wm/normalized/all_refseq_mammalian_viral_unique.tsv
```

It writes `00_source/source.sha256` to lock provenance and refuses to proceed if the source file is missing or later changes. **It does not redownload, recalculate, delete, move or edit the old data.** If the old normalized file's name differs, first inspect the original `normalized/` directory and deliberately update the single configuration file; do not silently substitute another dataset.

### Stage 01 — submit CPU protein clustering, split and leakage audit

After Nibi's Python and MMseqs2 are available on PATH:

```bash
bash scripts/submit_homology_cpu.sh
```

The command submits four CPU Slurm jobs under **`def-hup-ab`** in order (`afterok`): `prepare → search → split → audit`. No GPU job is submitted. Prepare/split request 4 CPUs; MMseqs2 search/audit request 32 CPUs, with 6 GB per CPU and up to 12 hours per job. Outputs and logs are always in the canonical stage folders above.

**Exact per-stage input/output manifests** (absolute paths, sizes, SHA-256 hashes, Git commit and Slurm Job ID) are saved only after successful execution in `/scratch/ylong66/data/mRNA_pr_wm/mRNAWM/00_source/stage_records/{prepare,search,split,audit}.json`. Every downstream stage reads the canonical paths defined in `scripts/config/paths.sh` and validates prerequisites. The output contract is:

| Completed stage | Exact outputs consumed by the next stage |
|---|---|
| Prepare | `/scratch/ylong66/data/mRNA_pr_wm/mRNAWM/01_split/homology_v1_seed3407/unique_proteins.fasta`; `prepare.json` |
| Search | `/scratch/ylong66/data/mRNA_pr_wm/mRNAWM/01_split/homology_v1_seed3407/all_vs_all.tsv`; `mmseqs_version.txt` |
| Split | `/scratch/ylong66/data/mRNA_pr_wm/mRNAWM/01_split/homology_v1_seed3407/{train,val,test}.tsv`; `protein_assignments.tsv`; `clusters.tsv`; `split_manifest.json` |
| Audit | `/scratch/ylong66/data/mRNA_pr_wm/mRNAWM/01_split/homology_v1_seed3407/leakage_audit.json`; `cross_split_hits.tsv` |

Verify after completion:

```bash
source scripts/config/paths.sh
cat "$SPLIT_ROOT/split_manifest.json"
cat "$SPLIT_ROOT/leakage_audit.json"
ls -lh "$SPLIT_ROOT"/{train,val,test}.tsv
ls -lh "$LOG_ROOT"
cat "$MRNAWM_ROOT/00_source/stage_records/split.json"
cat "$MRNAWM_ROOT/00_source/stage_records/audit.json"
```

Split protocol: 80/10/10, seed 3407, MMseqs2 alignment plus connected components; 30% identity + 80% reciprocal coverage, with an additional high-identity fragment guard. Cross-split alignment audit is mandatory. For methodological details and caveats, see [`docs/HOMOLOGY_SPLIT.md`](docs/HOMOLOGY_SPLIT.md). A clean heuristic search does not guarantee the absence of all remote homology.

### Stages 02–06 — model training, baselines and evaluation

**Not migrated or run yet.** The approved run sequence is: validate frozen split → prepare model records → build and check CodonTransformer feature coverage → train the *unchanged* mRNAWM → create matched test manifest → rerun *unchanged* baselines → evaluate on identical records/metrics → export results.

When these stages are migrated, **every runnable command goes in this README**, with its hard-coded canonical output path via `paths.sh`, exact Slurm script, resource request, checkpoint name and expected result filenames. Do not call old Dusky runners or mix old checkpoints with this experiment. Report all versions, settings and hashes before inspecting test outcomes.

## Standard Nibi Slurm headers

**CPU (default for new tasks, confirmed on Nibi):**
```bash
#!/bin/bash
#SBATCH --account=def-hup-ab
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=6G
#SBATCH --time=12:00:00
#SBATCH --output=slurm-%j.log
#SBATCH --error=slurm-%j.log
```

**GPU (reference only, user-provided header; `rrg-hup` membership currently denied on Nibi):**
```bash
#!/bin/bash
#SBATCH --account=rrg-hup
#SBATCH --gpus=h100_2g.20gb:1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --output=slurm-%j.log
#SBATCH --error=slurm-%j.log
```

The current authorized CPU account is `def-hup-ab`; the Nibi scheduler rejected `rrg-hup` for this user. The GPU header is retained solely as a reference, not as a confirmed usable allocation. Adjust CPU counts to a task's real parallelism. GPU example requests an H100 MIG slice, **not** the entire H100. Never carry `base/gpu/l40s` partitions or old Dusky account settings forward.

## Reproducibility rules

- **Never modify `mRNA_pr`; develop directly on `mRNAWM/main`.**
- Keep source download/normalization and model/baseline algorithms unchanged; only redesign the homology split and its validation.
- Freeze the normalized input checksum, MMseqs2 version/search settings, component IDs, split manifest, independent audit, Git SHA, training configuration, feature coverage, checkpoint IDs, and matched test manifest.
- **New stage = stage-specific folder + Slurm job + exact command documented here + output verification.** No abandoned paths or deprecated workflow steps.
- Unit tests: `python3 -m unittest discover -s tests -v`.

**Execution status:** Nibi workspace initialization succeeded (user log; source SHA-256 `cbedb5282fbc9a72e60c1b0469ffe5324f795a7908d7ceecc0c22c02e53fd2fe`). Initial CPU submission was blocked by missing MMseqs2; subsequent submission was rejected because `rrg-hup` was not authorized. CPU account corrected to `def-hup-ab` (maximum 12-hour request). No real-data split, training or evaluation has completed.
