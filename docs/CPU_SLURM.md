# Nibi CPU submission: homology split

The authoritative, current pipeline is the [root README](../README.md). All jobs run on Nibi, only on `mRNAWM/main`; no Dusky paths or partitions.

**Resource policy:** the user-confirmed account is `rrg-hup`. CPU-only jobs never request `--gpus`. Alignment search and independent audit are parallelized using **32 CPUs × 6 GiB/CPU = 192 GiB/node** and a 12-hour limit. Prepare and split are mostly single-threaded: override to **4 CPUs × 6 GiB** at submission. Requests can be adjusted after real Nibi profiling.

```bash
cd /path/to/mRNAWM
export DATA_ROOT="${SCRATCH:-$HOME/scratch}/mRNAWM"
test -s "$DATA_ROOT/normalized/all_refseq_mammalian_viral_unique.tsv"
command -v mmseqs
SCRIPT=scripts/slurm/submit_homology_split_cpu.slurm.sh

J1=$(sbatch --parsable --cpus-per-task=4 --export=ALL,STAGE=prepare "$SCRIPT")
J2=$(sbatch --parsable --dependency=afterok:$J1 --export=ALL,STAGE=search "$SCRIPT")
J3=$(sbatch --parsable --cpus-per-task=4 --dependency=afterok:$J2 --export=ALL,STAGE=split "$SCRIPT")
J4=$(sbatch --parsable --dependency=afterok:$J3 --export=ALL,STAGE=audit "$SCRIPT")
echo "prepare=$J1 search=$J2 split=$J3 audit=$J4"
squeue -u "$USER"
```

Each stage starts only if the previous stage succeeds. The 32-core MMseqs2 search may exceed 12 hours for a large all-vs-all RefSeq protein set; monitor logs and Slurm accounting, change time/memory/resources if necessary, and keep the same split protocol.

The output directory is `$DATA_ROOT/splits/homology_v1_seed3407/`. The essential validation artifacts are `split_manifest.json` and `leakage_audit.json`. An independent heuristic search cannot prove the absence of every remote homolog or protein-domain overlap.

The original `mRNA_pr` scripts used `def-hup-ab_cpu/base` or `def-hup-ab_gpu/gpu/l40s`, which must not be copied to Nibi. The GPU header supplied by the user is recorded in the root README strictly as a reference, not as a request to submit a GPU job.
