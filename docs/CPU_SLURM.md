# CPU Slurm submission (no GPU)

This project works directly on `main`; no project-specific working branch is needed.

From the **mRNAWM repository root** on the HPC login node, after the *unchanged* RefSeq download/parse/normalize/merge pipeline has created `$DATA_ROOT/normalized/all_refseq_mammalian_viral_unique.tsv`:

```bash
cd /home/ylong66/workspace/mRNAWM
export DATA_ROOT=/home/ylong66/scratch/data/mRNAWM
mkdir -p "$DATA_ROOT/splits/homology_v1_seed3407"
module load mmseqs2  # ONLY if your HPC provides this module; otherwise activate an env with mmseqs
# Verify the tool is available on compute nodes before submitting:
command -v mmseqs

SCRIPT=scripts/slurm/submit_homology_split_cpu.slurm.sh
J1=$(sbatch --parsable --export=ALL,STAGE=prepare "$SCRIPT")
J2=$(sbatch --parsable --dependency=afterok:$J1 --export=ALL,STAGE=search "$SCRIPT")
J3=$(sbatch --parsable --dependency=afterok:$J2 --export=ALL,STAGE=split "$SCRIPT")
J4=$(sbatch --parsable --dependency=afterok:$J3 --export=ALL,STAGE=audit "$SCRIPT")
echo "CPU jobs: prepare=$J1 search=$J2 split=$J3 audit=$J4"
squeue -u "$USER"
```

Every job requests **16 CPU cores, 64 GiB RAM, and 24 hours** with the old project's CPU account `def-hup-ab_cpu`; change the SBATCH account/partition/time for your own HPC if necessary. No `--gres=gpu`, CUDA, or GPU queue is requested. `module load mmseqs2` is illustrative, and depends on installed module naming; an already active MMseqs2 environment also works. Under SLURM the submission environment is normally inherited with `--export=ALL`.

Stage `audit` runs six independent directed cross-split searches. A successful search with zero hits is allowed, but missing files or nonzero exit status fail the job. If the split stage reports unexpectedly large connected components or a >5-percentage-point imbalance, investigate without weakening thresholds to force the split.

Outputs: `$DATA_ROOT/splits/homology_v1_seed3407/{train,val,test}.tsv`, `clusters.tsv`, `protein_assignments.tsv`, `split_manifest.json`, and `leakage_audit.json`. Reaching a passing audit here is a necessary (not sufficient) prerequisite for the downstream training and evaluation pipeline.

**Important:** no real-data Slurm job was submitted by ChatGPT. These are user-run submission commands; full-scale resource estimates should be checked on the HPC. If the MMseqs2 search is too large, increase time/memory or shard the computation while preserving a globally complete similarity graph.
