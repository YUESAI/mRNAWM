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

## Nibi — exact commands (quota-safe, restartable)

The old one-job all-vs-all search was **TIMEOUT** (Job `28555031`). The first sharded launcher then produced `AssocMaxSubmitJobLimit` because Slurm counts pending array tasks toward association submission quotas: `%6` limits concurrency, **not the number submitted**. **Both old launchers are retired.**

From the current `mRNAWM` Git directory on Nibi:

```bash
git pull --ff-only origin main
python3 -m unittest discover -s tests -v

# Check prior tasks before submitting new work.
squeue -u "$USER" -o "%.18i %.22j %.10T %.30R"

# Submit ONLY the next unfinished step or a bounded search batch.
bash scripts/submit_homology_cpu.sh
```

**Run the same final command again after its submitted job(s) finish.** This is intentionally step-by-step, not a batch of thousands of future jobs. It checks the actual output markers and automatically chooses `protein → graph_prepare → graph_search → graph_verify → split → audit_prepare → audit_search → audit_verify → audit`. It automatically skips successful shard results. No manual selection of stage or checkpoint paths is necessary.

- A *search stage* queues **up to four CPU workers**, not 256 or 1536 Slurm tasks; if Nibi rejects four, it tries two and then one. Each worker processes many shards sequentially and stops taking new shards before the 12-hour limit.
- Other stages submit **one job only**, under the confirmed account `def-hup-ab`.
- All jobs request `11:50:00` or less. Search workers request 12 CPUs / 72 GiB each; up to 48 CPUs active when four workers can be scheduled. A single-process stage requests 4 CPUs.
- A complete batch may still leave unfinished shards. Just run `bash scripts/submit_homology_cpu.sh` again after the batch exits; do not submit while the previous batch is running.
- If a single shard repeatedly exceeds its permitted time, **stop** and increase query sharding in a separate versioned run; no incomplete shard is marked complete.
- The launcher refuses to queue work while any job named `wm-homology` (old launcher) or `wm-shard-worker` (new launcher) is still queued/running. **Do not automatically cancel unrelated jobs.**
- If account quota is still full even for a one-worker array, it aborts with no future tasks queued. Check `squeue`, wait for capacity, and rerun the same command.

**Input:** `/scratch/ylong66/data/mRNA_pr_wm/mRNAWM/00_source/normalized.tsv` (read-only symlink).

**Outputs:** `/scratch/ylong66/data/mRNA_pr_wm/mRNAWM/01_split/homology_sharded_v2_seed3407/`.

**Logs:** `/scratch/ylong66/data/mRNA_pr_wm/mRNAWM/logs/`.

**Submission history (stage, Job ID, Git SHA, input/output paths):** `/scratch/ylong66/data/mRNA_pr_wm/mRNAWM/00_source/submissions.tsv`.

Note that the retired launcher's `00_source/last_submission.txt` might contain incomplete, failed queue submissions; it is **not** a completion report. Do not rely on it for the current run.

## Logs, job status, and output inspection

`07_logs` was an arbitrary directory-order number, **not pipeline stage 07**. New successful and failed jobs log to the clearer `mRNAWM/logs/` directory. Older failed logs in `07_logs/` are historical evidence and are not removed automatically.

```bash
source scripts/config/paths.sh
cat "$MRNAWM_ROOT/00_source/submissions.tsv"
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

**Execution status:** directory initialization succeeded; the old one-job search timed out and the first all-at-once array launcher hit Nibi submission limits. The new quota-safe launcher is committed but has **not** been run on Nibi through this chat. No final split or independent leakage audit is confirmed.
