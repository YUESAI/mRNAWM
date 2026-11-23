# mRNAWM — full-length mRNA design

**One active repository:** `YUESAI/mRNAWM`, directly on **main**. `YUESAI/mRNA_pr` is **read-only** and supplies reference model/baseline implementations. **All compute runs on Nibi**; never use Dusky-specific paths or resources. Keep data download/parsing and all existing model/baseline algorithms unchanged; replace the old exact-protein train/test grouping with a properly audited homology split. **Seed: 3407.**

## Current pipeline

| Stage | Work and outputs | Code | Status |
|---|---|---|---|
| 1. Download / normalize | RefSeq GBFF → `normalized/all_refseq_mammalian_viral_unique.tsv` | `scripts/data_pipeline/` | Legacy download/parser copied unchanged; Nibi run pending |
| 2. Protein homology search | unique FASTA + MMseqs2 all-vs-all alignments | `homology_split.py prepare/search` | Implemented, production run pending |
| 3. Cluster split + audit | train/val/test TSVs (80/10/10), component assignment, independent cross-split audit | `homology_split.py split/audit` | Implemented, production run pending |
| 4. Model data / feature cache | new split → model records, frozen CodonTransformer features | migrate proven code from `mRNA_pr` | Pending |
| 5. Train mRNAWM | new train/val → checkpoint | migrate unchanged model/training code | Pending |
| 6. Evaluate & run baselines | same frozen test IDs, model/baseline predictions and matched metrics | migrate unchanged baseline/evaluation code | Pending |

Do **not** start the production training or report benchmark results until the homology audit and feature-cache coverage pass. Old test metrics and split-specific cache directories are not the new benchmark.

## Nibi — one-time setup

On a **Nibi login node**, clone this repo and choose your actual data directory. `$SCRATCH` is used if the environment defines it; otherwise set `DATA_ROOT` explicitly. Verify `python3`, Biopython and `mmseqs` are available inside the submitted jobs.

```bash
cd /path/to/mRNAWM
export DATA_ROOT="${SCRATCH:-$HOME/scratch}/mRNAWM"
mkdir -p "$DATA_ROOT"
command -v mmseqs
python3 -c "import Bio; print(Bio.__version__)"
sacctmgr show associations user="$USER" format=Account,Partition -P
```

Set `NIBI_ACCOUNT` to an **account you can actually use**. The account identifiers and `base/gpu/l40s` resources found in the old `mRNA_pr` scripts are Dusky-specific examples, **not confirmed Nibi account/partition/GPU values**.

### 1. Download and normalize — same data pipeline

The code in `scripts/data_pipeline/` is copied unchanged from the working repository; only the execution location changes. Follow [active RefSeq processing commands](scripts/data_pipeline/DOWNLOADS.md). Required unsplit result:

```text
$DATA_ROOT/normalized/all_refseq_mammalian_viral_unique.tsv
```

### 2–3. Submit homology split & leakage audit — CPU by default

From the repository root, **only after the normalized input exists**, run these Nibi commands. Four jobs are linked with `afterok`, so a failed stage blocks its dependents. Default request: **8 CPU cores, 64 GiB RAM, 11h50m per stage**; this is a first pass, not a proven full-dataset MMseqs2 resource estimate.

```bash
export NIBI_ACCOUNT=YOUR_NIBI_ALLOCATION
SCRIPT=scripts/slurm/submit_homology_split_cpu.slurm.sh
J1=$(sbatch --parsable --account="$NIBI_ACCOUNT" --export=ALL,STAGE=prepare "$SCRIPT")
J2=$(sbatch --parsable --account="$NIBI_ACCOUNT" --dependency=afterok:$J1 --export=ALL,STAGE=search "$SCRIPT")
J3=$(sbatch --parsable --account="$NIBI_ACCOUNT" --dependency=afterok:$J2 --export=ALL,STAGE=split "$SCRIPT")
J4=$(sbatch --parsable --account="$NIBI_ACCOUNT" --dependency=afterok:$J3 --export=ALL,STAGE=audit "$SCRIPT")
echo "prepare=$J1 search=$J2 split=$J3 audit=$J4"
squeue -u "$USER"
```

Expected outputs at `$DATA_ROOT/splits/homology_v1_seed3407/`: `train.tsv`, `val.tsv`, `test.tsv`, `clusters.tsv`, `protein_assignments.tsv`, `split_manifest.json`, `leakage_audit.json`.

**Homology rules:** ≥30% aligned identity + ≥80% coverage of **both** proteins; or ≥90% identity, ≥80 amino acids aligned, and ≥80% coverage of **at least one** protein. All qualifying edges are joined into transitive connected components before splitting. The independent audit runs a separate six-direction train/val/test similarity search. These are proposed, testable thresholds; a negative heuristic search is not a proof that distant/short-domain homology is absent. Details: [homology split protocol](docs/HOMOLOGY_SPLIT.md).

## Slurm headers: follow old project layout, not old machine settings

The **structure** and CPU requests follow the working `mRNA_pr` scripts (`run_codonwm_v3_cache_cpu.sh`, `run_mrnawm_train.sh`, and its TAA splitting job). Those scripts still contain `def-hup-ab_cpu/base`, `def-hup-ab_gpu/gpu`, `--gres=gpu:l40s:1`, and old workspace/conda paths: **do not copy these to Nibi**.

**CPU (default):**

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=11:50:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
# submit with sbatch --account="$NIBI_ACCOUNT" script.sh
```

**GPU (reference template only; do not submit unless explicitly asked):**

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=11:50:00
#SBATCH --gpus-per-node=h100:1
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
# submit with sbatch --account="$NIBI_ACCOUNT" script.sh
```

The H100 syntax is a **proposed Nibi GPU template**, not a value verified from `mRNA_pr` or your Nibi allocation. Confirm permitted GPU resources and QoS on Nibi before requesting one.

## Non-negotiable experiment rules

- No edits to `mRNA_pr`; no branch creation in `mRNAWM` without request.
- Seed **3407** explicitly for split/training/sampling/evaluation; the old evaluation default was not 3407.
- Preserve all original transcript/ORF/UTR records within a protein homology component. Freeze input and output checksums and fixed held-out record IDs before training.
- Audit exact and homology overlaps, MMseqs2 coverage/search-cap limits, cluster sizes and dataset proportions. Record software versions, split manifest, feature coverage, checkpoint and Git SHA.
- Give a **CPU Slurm script plus exact sbatch submission commands** whenever writing code to run on Nibi; give GPU commands only on request.
- Unit tests: `python3 -m unittest discover -s tests -v`.

**Current execution status:** synthetic split tests passed in GitHub Actions; no production Nibi data split, model training, baseline generation or evaluation has run through this chat.
