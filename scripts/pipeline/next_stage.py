#!/usr/bin/env python3
"""Inspect successful artifacts and choose ONE next Nibi CPU job.

No Slurm submissions happen here. The task scheduler calls this probe before
queuing a single stage so that large arrays are never left pending.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

GRAPH_SHARDS = 256
AUDIT_SHARDS = 6 * GRAPH_SHARDS

def count(root: Path, kind: str) -> int:
    hits = root / ("graph_shards_v1" if kind == "graph" else "audit_shards_v1") / "hits"
    return sum(1 for _ in hits.rglob("*.done.json")) if hits.is_dir() else 0

def ok_json(path: Path, key: str | None = None) -> bool:
    if not path.is_file():
        return False
    try:
        value = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        return False
    return isinstance(value, dict) and (key is None or value.get(key) == "completed")

def next_step(root: Path) -> tuple[str, dict]:
    meta = {"graph_shards_completed": count(root, "graph"),
            "graph_shards_total": GRAPH_SHARDS,
            "audit_shards_completed": count(root, "audit"),
            "audit_shards_total": AUDIT_SHARDS}
    if not (root / "unique_proteins.fasta").is_file() or not (root / "prepare.json").is_file():
        return "protein", meta
    if not (root / "graph_shards_v1/READY.json").is_file():
        return "graph_prepare", meta
    if not ok_json(root / "graph_shards_v1/COMPLETE.json", "status"):
        if meta["graph_shards_completed"] < GRAPH_SHARDS:
            return "graph_search", meta
        return "graph_verify", meta
    if not all((root / name).is_file() for name in (
            "split_manifest.json", "train.tsv", "val.tsv", "test.tsv",
            "clusters.tsv", "protein_assignments.tsv",
            "train.proteins.fasta", "val.proteins.fasta", "test.proteins.fasta")):
        return "split", meta
    if not (root / "audit_shards_v1/READY.json").is_file():
        return "audit_prepare", meta
    if not ok_json(root / "audit_shards_v1/COMPLETE.json", "status"):
        if meta["audit_shards_completed"] < AUDIT_SHARDS:
            return "audit_search", meta
        return "audit_verify", meta
    report = root / "leakage_audit.json"
    if report.is_file():
        try:
            audit = json.loads(report.read_text())
        except (OSError, json.JSONDecodeError):
            return "AUDIT_FAILED", meta
        return ("DONE" if audit.get("passed") is True else "AUDIT_FAILED"), meta
    return "audit", meta

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--json", action="store_true")
    args = p.parse_args()
    step, stats = next_step(args.root)
    if args.json:
        print(json.dumps({"next_stage": step, **stats}))
    else:
        print(step)

if __name__ == "__main__":
    main()
