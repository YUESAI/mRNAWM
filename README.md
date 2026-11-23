# mRNAWM — reproducible full-length mRNA design

**Only active code: `YUESAI/mRNAWM/main`.** `mRNA_pr` (including its `v1` branch) is read-only reference. All computing happens on **Nibi**. **Seed = 3407.** Training and baseline algorithms remain unchanged; only dataset homology splitting and consequent reruns are being rebuilt.

## Canonical storage and pipeline

Source RefSeq data (already downloaded and normalized; **do not overwrite**):

```text
/scratch/ylong66/data/mRNA_pr_wm/normalized/all_refseq_mammalian_viral_unique.tsv
```

**All new experiments** are isolated at `/scratch/ylong66/data/mRNA_pr_wm/mRNAWM/`.

| Step | Exact input | Exact output | Status |
|---|---|---|---|
| 00. Input freeze | Previous `normalized/all_refseq_mammalian_viral_unique.tsv` | `mRNAWM/00_source/normalized.tsv` (symlink), `source.sha256` | Completed on Nibi (user log) |
| 01A. Unique proteins | `00_source/normalized.tsv` | `01_split/homology_sharded_v2_seed3407/unique_proteins.fasta`, `prepare.json` | Code ready |
| 01B. Query shards + global DB | `unique_proteins.fasta` | `01_split/homology_sharded_v2_seed3407/graph_shards_v1/READY.json`, query FASTAs and full target DB | Code ready |
| 01C. Homology search | 256 query shards × complete target DB | `graph_shards_v1/hits/full/part_0000.tsv` … `part_0255.tsv`, `.done.json` per shard | Code ready |
| 01D. Verify search | All 256 complete hashed shard outputs | `graph_shards_v1/COMPLETE.json`, `hits.list` | Code ready |
| 01E. Split | `00_source/normalized.tsv` + verified `hits.list` | `01_split/homology_sharded_v2_seed3407/{train,val,test}.tsv`, `split_manifest.json`, `clusters.tsv`, `protein_assignments.tsv` | Code ready |
| 01F. Audit queries/DB | `{train,val,test}.proteins.fasta` | `audit_shards_v1/READY.json`, complete per-split target DBs | Code ready |
| 01G. Audit searches | Six directed cross-split searches, each with 256 query shards | `audit_shards_v1/hits/<split_pair>/part_*.tsv`, completion markers | Code ready |
| 01H. Final audit | Verified 1536 audit shard outputs + split TSVs | `audit_shards_v1/COMPLETE.json`, `hits.list`, `leakage_audit.json` | Code ready |
| 02. Model records | Final audited new split | `mRNAWM/02_model_data/` | Not migrated |
| 03. Frozen protein features | Stage 02 proteins | `mRNAWM/03_feature_cache/` | Not migrated |
| 04. mRNAWM training | New train/val + feature cache | `mRNAWM/04_checkpoints/` | Not migrated |
| 05. Baselines | New fixed test manifest | `mRNAWM/05_baselines/` | Not migrated |
| 06. Evaluation | Same test IDs for model and baselines | `mRNAWM/06_evaluations/` | Not migrated |

Every completed non-array step writes input/output paths, SHA-256, Git commit and Slurm Job ID to `mRNAWM/00_source/stage_records/*_v2.json`. **Only a passing final audit qualifies the split for model training.**

## Nibi — copy/paste commands from your cloned repo

The old one-job all-vs-all MMseqs2 search **timed out** (Job `28555031` on 2026-11-23). It is **not part of the active workflow**. This new run uses a fresh `homology_sharded_v2_seed3407` output directory and does not mix partially produced old results.

Prepare the environment with `mmseqs` and `python3` on PATH. Then, while inside the Git repo:

```bash
git pull --ff-only origin main
bash scripts/init_mrnawm_workspace.sh
python3 -m unittest discover -s tests -v
bash scripts/submit_homology_cpu.sh
```

The final command submits the ordered stages with `afterok` dependencies: protein FASTA → global DB/shards → **256 search array** → strict verify → split → cross-split DB/shards → **six 256-task audit arrays** → strict verify → final audit.

**CPU allocation:** `--account=def-hup-ab`, **11:50:00** per task (Nibi 12-hour maximum). Parallel search/audit: 12 CPU cores × 6 GiB/core per task, up to 6 concurrent tasks (72 CPUs and 432 GiB of aggregate requested memory). The audit is submitted as six 256-task arrays (one active task per audit pair) to stay below scheduler array-size limits. Single-process prepare/verify/split tasks request 4 CPU cores. Adjust concurrency only after observing Nibi scheduling and memory usage. No GPU is requested.

**Retry search after a shard times out:** completed shard checksums are verified and successful shards are skipped; no successful shards are deleted.

```bash
bash scripts/submit_homology_cpu.sh resume-graph
```

**Retry audit after a shard times out:**

```bash
bash scripts/submit_homology_cpu.sh resume-audit
```

Do not submit either retry mode while the same array is still running. If a *single* shard repeatedly times out, reduce its query count (new versioned settings) rather than claiming its partial result is complete. If the split itself fails due to giant transitive components or imbalance, inspect the clusters instead of relaxing the thresholds to force a split.

## Logs, job status, and output inspection

`07_logs` was an arbitrary directory-order number, **not pipeline stage 07**. New successful and failed jobs log to the clearer `mRNAWM/logs/` directory. Older failed logs in `07_logs/` are historical evidence and are not removed automatically.

```bash
source scripts/config/paths.sh
cat "$MRNAWM_ROOT/00_source/last_submission.txt"
ls -lhtr "$LOG_ROOT" | tail -20
sacct -u "$USER" -S 2026-11-23 \
  --format=JobID,JobName%22,State%20,ExitCode,Elapsed,Timelimit,MaxRSS
cat "$SPLIT_ROOT/graph_shards_v1/COMPLETE.json"
cat "$SPLIT_ROOT/audit_shards_v1/COMPLETE.json"
cat "$SPLIT_ROOT/leakage_audit.json"
```

The last three files only exist after their stages finish. `COMPLETED` plus a recorded pass is success; `TIMEOUT`, `FAILED`, `OUT_OF_MEMORY` or missing success markers is not. The Slurm submission script saves the exact job IDs under `00_source/last_submission.txt`.

## Methodological integrity

- **Search graph:** all query shards align against the **complete global protein database** (never independently cluster the shards). Qualifying alignment edges are combined into global connected components.
- **Scientific filtering:** primary identity ≥30% plus 80% coverage of both sequences; fragment guard identity ≥90%, aligned length ≥80 aa, ≥80% coverage of one sequence. Unlike the `mRNA_pr/v1` reciprocal-only 30% protocol, these secondary edges are deliberate and may change component sizes.
- **Audit:** independent all-direction cross-split similarity searches, with the same scientific thresholds. Check `leakage_audit.json` `"passed": true`. MMseqs2 is heuristic and `--max-seqs` finite; successful audits do not prove the absence of all remote or domain-level homology.
- Save frozen input checksums, versions, full set of completed query shard manifests, and all split outputs. The final test set must be used only after the split protocol is fixed.
- GPU jobs are **not requested** unless you explicitly ask. The user-provided GPU header is reference-only (`--gpus=h100_2g.20gb:1`) and cannot be used with the unauthorized `rrg-hup` account without allocation permission.

**Execution status:** code is committed, but the new sharded experiment has **not yet been submitted or completed on Nibi**. The earlier job was TIMEOUT; no old result is treated as a successful cluster.
