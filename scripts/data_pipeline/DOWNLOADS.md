# Existing RefSeq data — read-only provenance

**Normal mRNAWM reruns do not download or re-parse RefSeq.** The original data is already present under:

```text
/scratch/ylong66/data/mRNA_pr_wm/raw/
/scratch/ylong66/data/mRNA_pr_wm/normalized/
```

The expected merged normalized TSV is:

```text
/scratch/ylong66/data/mRNA_pr_wm/normalized/all_refseq_mammalian_viral_unique.tsv
```

That filename comes from the old project's normalized merge command. The workspace initializer checks that it **actually exists**; it will refuse to proceed rather than silently use another file. If the file differs, inspect the Nibi directory first and update the single source path in `scripts/config/paths.sh`. Do not rewrite or reprocess old experimental data.

Run from the cloned mRNAWM Git repository:

```bash
bash scripts/init_mrnawm_workspace.sh
```

This stores a symlink and SHA-256 source manifest under:

```text
/scratch/ylong66/data/mRNA_pr_wm/mRNAWM/00_source/
```

All newly calculated products start in `mRNAWM/01_split/`, **not** in the old `splits/` folder. Exact submission steps and later stage paths are maintained in the [root README](../../README.md).

## Rebuild policy

The original downloader and parser Python programs remain in `scripts/data_pipeline/` unchanged for audit and provenance. Re-downloading or re-normalizing is **not part of the current active workflow** and must not happen without explicit instruction. The existing original raw/normalized files are read-only inputs.
