# mRNAWM

Reproducible full-length mRNA world-model project.

**Current scope:** The legacy working model and baseline implementations live in the read-only reference repository `YUESAI/mRNA_pr`. This repository is being rebuilt carefully starting from the dataset split. No changes to existing model algorithms or data download sources are intended.

## Start here

Read [homology-aware protein split protocol](docs/HOMOLOGY_SPLIT.md), including exact MMseqs2 commands, seed **3407**, input and output schema, independent leakage audit, and re-run requirements.

For CPU-only Slurm submissions, see [ready-to-run sbatch commands](docs/CPU_SLURM.md). All development goes directly onto **main**; do not create additional branches without explicit instruction.

The legacy approach used SHA-256 of the exact protein amino-acid sequence as its "cluster"; the new approach groups proteins using similarity-based connected components, with an additional partial-sequence guard.

~~~bash
python -m unittest discover -s tests -v
~~~

**Execution status:** code and synthetic tests have been committed directly to main; full-scale data clustering, independent validation, training and matched baselines have NOT yet been run. Nothing in the documentation should be interpreted as a completed experiment.

## Migration rules

1. Only edit this repository; keep `mRNA_pr` unchanged.
2. Preserve the legacy download/parse pipeline and original model and baseline algorithms. Audit and copy them as separate reviewed migration changes.
3. Version and checksum normalized data, split TSVs, feature caches, checkpoints, baseline outputs and evaluation manifests. Do not report prior test scores as new results.
4. Use `seed=3407` at every split, training, sampling and evaluation entrypoint, explicitly overriding other defaults.
