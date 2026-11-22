# Protein homology split protocol — seed 3407

## Scope
Only YUESAI/mRNAWM is modified. YUESAI/mRNA_pr is reference-only. Keep the existing RefSeq download, parser, normalized schema, mRNAWM model architecture, hyperparameters, baseline algorithms and evaluation metrics unchanged. Data splits, frozen feature-cache coverage, training runs, baseline outputs and all evaluations must be regenerated on the new split.

## Root cause
The old scripts/data_pipeline/make_exact_splits.py assigns sha256(protein_sequence)[:16], i.e. exact-protein grouping, not homology clustering. Even one amino-acid substitution changes the group.

## Proposed, versioned primary protocol
- Seed 3407; 80/10/10 train/val/test. Balance by number of eligible complete-mRNA records, not number of clusters.
- Exact amino-acid sequences first normalized consistently with protein_key in the model: uppercase, '_' converted to stop, trailing stops removed; full SHA-256 protein IDs.
- Sensitive all-against-all MMseqs2 protein alignment. For each hit, make an undirected edge when (identity >=30% and min(query coverage,target coverage) >=80%), OR (identity >=90%, aligned length >=80 amino acids and max coverage >=80%). The second condition guards against high-identity shorter isoforms / partial containment.
- Compute transitive connected components of all passing hits BEFORE any random split. All records for a protein travel together. Cluster ID is lexicographically lowest full SHA-256 inside its component.
- Preserve original normalized TSV schema plus protein_cluster_id and split. Do not discard synonymous CDS alternatives or UTR records.
- Allocation is stable largest-first with seeded tie-breaking and record-count targets. Fail on large imbalance rather than silently changing sequence thresholds.
- These are *proposed* thresholds to be audited and justified against experimental deployment. A domain pair with <80% mutual coverage can escape; distant homology may escape heuristic searches.

## HPC commands
Run after the UNCHANGED download/parse/normalize/merge pipeline completes (input should be unsplit normalized TSV):

~~~bash
export DATA_ROOT=/path/to/mRNAWM_data
export INPUT="$DATA_ROOT/normalized/all_refseq_mammalian_viral_unique.tsv"
export WORK="$DATA_ROOT/splits/homology_v1_seed3407"
python scripts/data_pipeline/homology_split.py prepare --input "$INPUT" --output-dir "$WORK"
mmseqs version > "$WORK/mmseqs_version.txt"
mkdir -p "$WORK/tmp"
mmseqs easy-search "$WORK/unique_proteins.fasta" "$WORK/unique_proteins.fasta" "$WORK/all_vs_all.tsv" "$WORK/tmp" --format-output "query,target,pident,alnlen,qlen,tlen,qcov,tcov" --min-seq-id 0.25 -c 0.0 -s 7.5 --max-seqs 100000 --threads 16
python scripts/data_pipeline/homology_split.py split --input "$INPUT" --prepare-manifest "$WORK/prepare.json" --hits "$WORK/all_vs_all.tsv" --output-dir "$WORK" --seed 3407
~~~

The MMseqs2 pident column is 0–100%; qcov/tcov are 0–1. The program checks this. Keep MMseqs2 version, commands, alignment file and checksums. Never use mmseqs easy-cluster representative membership directly without a cross-cluster leakage check.

## Mandatory independent cross-split check
After splitting, independently run MMseqs2 high-sensitivity searches in BOTH directions for each train/val, train/test and val/test pairing, and concatenate the eight-column TSV results to cross_split_hits.tsv. Keep all underlying search commands, logs, hit counts, exit codes, and verify that they completed without hitting max-seqs truncation.

~~~bash
python scripts/data_pipeline/homology_split.py audit --train "$WORK/train.tsv" --val "$WORK/val.tsv" --test "$WORK/test.tsv" --hits "$WORK/cross_split_hits.tsv" --report "$WORK/leakage_audit.json"
~~~

If there are legitimately zero hits, cross_split_hits.tsv must contain precisely this tab-separated header:
query target pident alnlen qlen tlen qcov tcov

The audit also checks exact protein, group and record-ID intersection. The audit *does not* prove the absence of remote homologs or unseen local sequence similarities: MMseqs2 is a heuristic and a capped search is non-exhaustive. Record maximum train/test identity over all detected candidates, nearest-neighbor identity distributions, passing-hit counts, search-cap saturation and domain-related hard cases. Consider secondary DIAMOND or BLAST checks. Do not claim 'zero homology leakage' without qualification.

## Downstream reproducibility gates
1. Freeze input SHA-256, protein count, normalized pair count, eligibility count, component count and distribution, oversize components, seed, MMseqs2 binary version and options, all hit checksums, group assignments, balance and independent audit.
2. Regenerate test manifest and all ORF JSONL / full-mRNA training inputs. Recreate or reuse immutable protein-key CodonTransformer features only after auditing new-split coverage; old split-specific cache labels are stale.
3. Retain the original mRNAWM model implementation and training hyperparameters; pass --train-input/--val-input on new splits and --seed 3407 explicitly.
4. The old evaluate_mrnawm.py defaults to --seed 9917, so override it with --seed 3407 for the new evaluation runs.
5. Run all original baselines on the SAME frozen test record IDs, with the same budgets, sequence constraints, model versions and seed. Recalculate all metrics and compare on matched records only.
6. Freeze the protocol before inspecting final test metrics. Report dataset and split versions in the paper. Store per-run config and Git commit hash.

## Related work
- Steinegger & Söding, MMseqs2, Nature Biotechnology 2017: https://doi.org/10.1038/nbt.3988
- AlQuraishi, ProteinNet, BMC Bioinformatics 2019: https://doi.org/10.1186/s12859-019-2932-0
- ProteinNet splitting methodology: https://github.com/aqlaboratory/proteinnet/blob/master/docs/splitting_methodology.md
ProteinNet uses a more elaborate sequence/MSA-aware procedure; it is an inspiration, not a claim that the thresholds above are taken directly from that paper.

## Execution status
This commit does NOT mean the production dataset was split. The private full RefSeq data and MMseqs2 HPC runtime were unavailable in this session. No production clusters, coverage report, model run or baselines have yet been generated.
