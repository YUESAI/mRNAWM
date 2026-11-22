#!/usr/bin/env python3
"""Reproducible homology-component split of normalized protein/mRNA TSV records.

Input is the unchanged normalized TSV produced by the legacy RefSeq pipeline.
MMseqs2 supplies alignment hits; this script NEVER infers similarity from hashes.
Only stdlib dependencies are used here. Seed defaults to 3407.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

SPLITS = ("train", "val", "test")
HIT_COLUMNS = ("query", "target", "pident", "alnlen", "qlen", "tlen", "qcov", "tcov")
AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYXBUZOJ")


def csv_limits() -> None:
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit //= 10


def protein_normalize(value: str) -> str:
    # Match the legacy CodonWM feature-store protein_key normalization.
    seq = (value or "").strip().upper().replace("_", "*").rstrip("*")
    if not seq or any(aa not in AA_ALPHABET for aa in seq):
        raise ValueError(f"invalid or missing protein sequence: {seq[:60]!r}")
    return seq


def protein_id(sequence: str) -> str:
    return hashlib.sha256(protein_normalize(sequence).encode("ascii")).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path, required: set[str]):
    csv_limits()
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        if not fields or not required.issubset(fields):
            raise ValueError(f"{path}: missing columns {sorted(required - set(fields or []))}")
        for line_number, row in enumerate(reader, 2):
            try:
                yield fields, row
            except Exception as exc:
                raise ValueError(f"{path}:{line_number}: {exc}") from exc


def eligible(row: dict, protein: str) -> bool:
    """Approximate the unchanged full-mRNA trainer's default eligibility."""
    return bool(row.get("utr5_sequence") and row.get("utr3_sequence")
                and 1 <= len(protein) <= 1022 and row.get("orf_sequence"))


def prepare(args) -> None:
    source, output = Path(args.input), Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    fasta = output / "unique_proteins.fasta"
    seen = set()
    rows = full_length_rows = 0
    with fasta.open("w", encoding="ascii") as dst:
        for _, row in read_tsv(source, {"protein_sequence", "orf_sequence", "utr5_sequence", "utr3_sequence"}):
            seq = protein_normalize(row["protein_sequence"])
            pid = protein_id(seq)
            rows += 1
            full_length_rows += int(eligible(row, seq))
            if pid not in seen:
                seen.add(pid)
                dst.write(f">{pid}\n{seq}\n")
    if not rows:
        raise ValueError("empty normalized input TSV")
    manifest = {"input_sha256": sha256_file(source), "input_rows": rows,
                "eligible_rows": full_length_rows, "unique_proteins": len(seen),
                "fasta_sha256": sha256_file(fasta)}
    (output / "prepare.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


class UnionFind:
    def __init__(self, ids):
        self.parent = {p: p for p in ids}
        self.size = {p: 1 for p in ids}

    def find(self, p):
        while self.parent[p] != p:
            self.parent[p] = self.parent[self.parent[p]]
            p = self.parent[p]
        return p

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a == b:
            return False
        if self.size[a] < self.size[b] or (self.size[a] == self.size[b] and a > b):
            a, b = b, a
        self.parent[b] = a
        self.size[a] += self.size[b]
        del self.size[b]
        return True


def alignment_hits(paths, ids, main_id, main_cov, fragment_id, fragment_cov, fragment_aa):
    """Yield qualifying pairs, checking MMseqs2 output identifiers and ranges."""
    for path_string in paths:
        path = Path(path_string)
        if not path.exists() or path.stat().st_size == 0:
            raise ValueError(f"missing/empty MMseqs2 hits: {path}")
        with path.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, 1):
                cols = line.rstrip("\n").split("\t")
                if cols == list(HIT_COLUMNS):  # optional header
                    continue
                if len(cols) != 8:
                    raise ValueError(f"{path}:{line_number}: expected {HIT_COLUMNS}")
                q, t = cols[:2]
                if q not in ids or t not in ids:
                    raise ValueError(f"{path}:{line_number}: protein not in input: {q}, {t}")
                try:
                    ident, aln_len, qlen, tlen, qcov, tcov = map(float, cols[2:])
                except ValueError as exc:
                    raise ValueError(f"{path}:{line_number}: invalid alignment values") from exc
                if not (0 <= ident <= 100 and 0 <= qcov <= 1 and 0 <= tcov <= 1
                        and aln_len > 0 and qlen > 0 and tlen > 0):
                    raise ValueError(f"{path}:{line_number}: out-of-range alignment values")
                if q == t:
                    continue
                identity = ident / 100.0  # MMseqs2 pident is a percentage, unlike qcov/tcov.
                full = identity >= main_id and min(qcov, tcov) >= main_cov
                fragment = (identity >= fragment_id and aln_len >= fragment_aa
                            and max(qcov, tcov) >= fragment_cov)
                if full or fragment:
                    yield q, t, "full" if full else "fragment"


def collect_input(path):
    counts, full_counts = Counter(), Counter()
    fields = None
    rows = 0
    for fields, row in read_tsv(path, {"protein_sequence", "orf_sequence", "utr5_sequence", "utr3_sequence"}):
        seq = protein_normalize(row["protein_sequence"])
        pid = protein_id(seq)
        counts[pid] += 1
        full_counts[pid] += int(eligible(row, seq))
        rows += 1
    if not rows:
        raise ValueError("empty normalized input")
    if "split" in fields or "protein_cluster_id" in fields:
        raise ValueError("Expected UNSPLIT normalized input, not an existing split")
    return fields, counts, full_counts


def choose_splits(cluster_weights, train, val, seed):
    """Largest-first group allocation with seeded, stable tie breaking."""
    target = {"train": train, "val": val, "test": 1 - train - val}
    assigned = Counter()
    rng = random.Random(seed)
    components = sorted(cluster_weights)
    rng.shuffle(components)
    components.sort(key=lambda cluster: -cluster_weights[cluster])  # stable ties
    total = sum(cluster_weights.values())
    choices = {}
    for cluster in components:
        split = max(SPLITS, key=lambda name: target[name] * total - assigned[name])
        choices[cluster] = split
        assigned[split] += cluster_weights[cluster]
    return choices


def split_dataset(args) -> None:
    source, output = Path(args.input), Path(args.output_dir)
    if args.seed != 3407:
        raise ValueError("The registered split seed is fixed to 3407.")
    if not (0 < args.train_fraction < 1 and 0 < args.val_fraction < 1
            and args.train_fraction + args.val_fraction < 1):
        raise ValueError("invalid split fractions")
    fields, counts, full_counts = collect_input(source)
    prepared = json.loads(Path(args.prepare_manifest).read_text())
    if (sha256_file(source) != prepared["input_sha256"]
            or sum(counts.values()) != prepared["input_rows"]
            or len(counts) != prepared["unique_proteins"]):
        raise ValueError("input differs from prepare.json; regenerate FASTA and hits")
    uf = UnionFind(counts)
    passing = Counter()
    unions = 0
    for q, t, label in alignment_hits(args.hits, counts,
            args.min_identity, args.min_coverage,
            args.fragment_identity, args.fragment_coverage, args.fragment_min_aa):
        passing[label] += 1
        unions += int(uf.union(q, t))
    if len(counts) >= 20 and unions == 0:
        raise ValueError("No similar proteins grouped: check MMseqs2 output, thresholds, and IDs")

    # Connected components, NOT centroid membership. The representative key is
    # the lexically smallest full SHA-256, independent of MMseqs row order.
    root_by_pid = {p: uf.find(p) for p in counts}
    representative = {}
    for pid, root in root_by_pid.items():
        if root not in representative or pid < representative[root]:
            representative[root] = pid
    component_by_pid = {pid: representative[root] for pid, root in root_by_pid.items()}
    pairs = Counter()
    possible = Counter()
    proteins = Counter()
    for pid, cluster in component_by_pid.items():
        pairs[cluster] += counts[pid]
        possible[cluster] += full_counts[pid]
        proteins[cluster] += 1
    if len(possible) < 3 or sum(possible.values()) == 0:
        raise ValueError("Need >=3 homology components and eligible full-mRNA records")
    weights = possible if args.balance_on == "eligible" else pairs
    choices = choose_splits(weights, args.train_fraction, args.val_fraction, args.seed)
    actual = Counter()
    for component, label in choices.items():
        actual[label] += weights[component]
    total_weight = sum(weights.values())
    expected = {"train": args.train_fraction, "val": args.val_fraction,
                "test": 1 - args.train_fraction - args.val_fraction}
    deviation = {name: abs(actual[name] / total_weight - expected[name]) for name in SPLITS}
    if max(deviation.values()) > args.max_fraction_deviation:
        raise ValueError(f"Unbalanced components: deviations={deviation}; "
                         "inspect oversized homology components before splitting")

    output.mkdir(parents=True, exist_ok=True)
    handles = {}
    protein_fasta = {}
    try:
        for name in SPLITS:
            handles[name] = (output / f"{name}.tsv").open("w", newline="", encoding="utf-8")
            protein_fasta[name] = (output / f"{name}.proteins.fasta").open("w", encoding="ascii")
        writers = {name: csv.DictWriter(handles[name], fieldnames=fields + ["protein_cluster_id", "split"],
                                         delimiter="\t") for name in SPLITS}
        for writer in writers.values():
            writer.writeheader()
        seen = set()
        row_count = Counter()
        eligible_count = Counter()
        for _, row in read_tsv(source, {"protein_sequence"}):
            seq = protein_normalize(row["protein_sequence"])
            pid = protein_id(seq)
            cluster = component_by_pid[pid]
            label = choices[cluster]
            writerow = dict(row, protein_cluster_id=cluster, split=label)
            writers[label].writerow(writerow)
            row_count[label] += 1
            eligible_count[label] += int(eligible(row, seq))
            if pid not in seen:
                seen.add(pid)
                protein_fasta[label].write(f">{pid}\n{seq}\n")
    finally:
        for f in list(handles.values()) + list(protein_fasta.values()):
            f.close()

    with (output / "protein_assignments.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["protein_sha256", "protein_cluster_id", "split", "pair_count", "eligible_pair_count"])
        for pid in sorted(counts):
            cluster = component_by_pid[pid]
            writer.writerow([pid, cluster, choices[cluster], counts[pid], full_counts[pid]])
    with (output / "clusters.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["protein_cluster_id", "split", "n_proteins", "n_pairs", "n_eligible_pairs"])
        for cluster in sorted(pairs):
            writer.writerow([cluster, choices[cluster], proteins[cluster], pairs[cluster], possible[cluster]])

    singleton = sum(size == 1 for size in proteins.values())
    summary = {
        "seed": args.seed, "input_sha256": prepared["input_sha256"],
        "prepare_manifest": str(args.prepare_manifest),
        "hits_sha256": {str(path): sha256_file(Path(path)) for path in args.hits},
        "settings": {"min_identity": args.min_identity, "min_coverage": args.min_coverage,
                     "fragment_identity": args.fragment_identity,
                     "fragment_coverage": args.fragment_coverage,
                     "fragment_min_aa": args.fragment_min_aa,
                     "balance_on": args.balance_on, "train_fraction": args.train_fraction,
                     "val_fraction": args.val_fraction},
        "rows": dict(row_count), "eligible_rows": dict(eligible_count),
        "proteins": len(counts), "clusters": len(pairs),
        "singleton_clusters": singleton, "singleton_fraction": singleton / len(pairs),
        "largest_cluster_proteins": max(proteins.values()),
        "largest_cluster_pairs": max(pairs.values()),
        "qualifying_alignment_rows": dict(passing),
        "connected_component_merges": unions,
        "achieved_balance": {s: actual[s] / total_weight for s in SPLITS},
        "absolute_balance_deviation": deviation,
    }
    (output / "split_manifest.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def audit(args) -> None:
    proteins = {}
    clusters = {}
    record_ids = {}
    counts = Counter()
    errors = []
    for split in SPLITS:
        for _, row in read_tsv(Path(getattr(args, split)), {"protein_sequence", "protein_cluster_id", "split", "record_id"}):
            pid = protein_id(row["protein_sequence"])
            cluster = row["protein_cluster_id"]
            record_id = row["record_id"]
            counts[split] += 1
            if row["split"] != split or not cluster:
                errors.append(f"{split}: invalid split/cluster column")
            if pid in proteins and proteins[pid] != split:
                errors.append(f"exact protein leakage: {pid}")
            if cluster in clusters and clusters[cluster] != split:
                errors.append(f"cluster leakage: {cluster}")
            if record_id and record_id in record_ids and record_ids[record_id] != split:
                errors.append(f"record_id leakage: {record_id}")
            proteins[pid] = split
            clusters[cluster] = split
            if record_id:
                record_ids[record_id] = split
    seen = Counter()
    samples = []
    for q, t, kind in alignment_hits(args.hits, proteins,
            args.min_identity, args.min_coverage,
            args.fragment_identity, args.fragment_coverage, args.fragment_min_aa):
        if proteins[q] != proteins[t]:
            seen[kind] += 1
            if len(samples) < 20:
                samples.append({"query": q, "target": t, "kind": kind,
                                "splits": [proteins[q], proteins[t]]})
    if any(seen.values()):
        errors.append("Qualifying cross-split homology alignment detected")
    report = {"seed": 3407, "rows": dict(counts), "unique_proteins": len(proteins),
              "cross_split_alignment_rows": dict(seen), "samples": samples,
              "errors": errors[:100], "passed": not errors,
              "audit_hits_sha256": {str(p): sha256_file(Path(p)) for p in args.hits}}
    out = Path(args.report)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if errors:
        raise SystemExit(1)


def common_hit_args(parser):
    parser.add_argument("--hits", action="append", required=True, help="MMseqs2 TSV; repeat for multiple searches.")
    parser.add_argument("--min-identity", type=float, default=0.30)
    parser.add_argument("--min-coverage", type=float, default=0.80)
    parser.add_argument("--fragment-identity", type=float, default=0.90)
    parser.add_argument("--fragment-coverage", type=float, default=0.80)
    parser.add_argument("--fragment-min-aa", type=int, default=80)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--input", required=True)
    p.add_argument("--output-dir", required=True)
    p.set_defaults(func=prepare)
    p = sub.add_parser("split")
    p.add_argument("--input", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--prepare-manifest", required=True)
    p.add_argument("--train-fraction", type=float, default=0.80)
    p.add_argument("--val-fraction", type=float, default=0.10)
    p.add_argument("--seed", type=int, default=3407)
    p.add_argument("--balance-on", choices=("eligible", "all"), default="eligible")
    p.add_argument("--max-fraction-deviation", type=float, default=0.05)
    common_hit_args(p)
    p.set_defaults(func=split_dataset)
    p = sub.add_parser("audit")
    for name in SPLITS:
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--report", required=True)
    common_hit_args(p)
    p.set_defaults(func=audit)
    args = parser.parse_args()
    if args.command in ("split", "audit"):
        if not (0 <= args.min_identity <= 1 and 0 <= args.min_coverage <= 1
                and 0 <= args.fragment_identity <= 1 and 0 <= args.fragment_coverage <= 1
                and args.fragment_min_aa > 0):
            parser.error("homology thresholds out of range")
    args.func(args)


if __name__ == "__main__":
    main()
