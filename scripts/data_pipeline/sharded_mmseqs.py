#!/usr/bin/env python3
"""Checkpointed, global-target MMseqs2 searches split by QUERY only.

Graph: each of N query shards searches the FULL protein target DB.
Audit: each of N query shards per split searches each of the TWO other
entire split target DBs. Shards never define independent homology universes.

Checkpoint success requires a complete output, hashed inputs, and a .done.json.
No incomplete result is ever promoted to a finished shard.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

FORMAT = "query,target,pident,alnlen,qlen,tlen,qcov,tcov"
HEADER = FORMAT.replace(",", "\t") + "\n"
PAIRS = ("train_val", "train_test", "val_train", "val_test", "test_train", "test_val")
VERSION = "query-global-target-v1"
SETTINGS = {
    "min_seq_id": 0.25, "min_coverage": 0.0, "sensitivity": 7.5,
    "evalue": 0.001, "max_seqs": 100000, "format": FORMAT
}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_json(path: Path):
    return json.loads(path.read_text())


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        temp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)


def fasta_iter(fasta: Path):
    name, parts = None, []
    with fasta.open(encoding="ascii") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                if name is None:
                    raise ValueError("invalid FASTA: missing header")
                parts.append(line)
    if name is not None:
        yield name, "".join(parts)


def shard_fasta(source: Path, queries: Path, n: int):
    queries.mkdir(parents=True, exist_ok=True)
    outputs = [queries / f"part_{i:04d}.faa" for i in range(n)]
    partials = [p.with_suffix(".partial") for p in outputs]
    handles = [p.open("w", encoding="ascii") for p in partials]
    seen = set()
    counts = [0] * n
    try:
        for protein_id, seq in fasta_iter(source):
            if (not seq or len(protein_id) != 64 or
                    hashlib.sha256(seq.encode("ascii")).hexdigest() != protein_id or
                    protein_id in seen):
                raise ValueError(f"invalid, duplicated or mismatched protein: {protein_id}")
            seen.add(protein_id)
            i = int(protein_id[:16], 16) % n
            handles[i].write(f">{protein_id}\n{seq}\n")
            counts[i] += 1
    finally:
        for h in handles:
            h.close()
    if not seen:
        raise ValueError(f"empty FASTA {source}")
    for tmp, final in zip(partials, outputs):
        tmp.replace(final)
    return {"proteins": len(seen), "counts": counts, "source_sha256": sha(source)}


def base_path(args):
    return args.root / ("graph_shards_v1" if args.kind == "graph" else "audit_shards_v1")


def prepare(args):
    base = base_path(args)
    prepared = base / "READY.json"
    if args.kind == "graph":
        prep = read_json(args.root / "prepare.json")
        fasta = args.root / "unique_proteins.fasta"
        if prep["fasta_sha256"] != sha(fasta):
            raise RuntimeError("FASTA does not match prepare manifest")
        sources = {"full": {"path": str(fasta.resolve()), "sha256": sha(fasta),
                            "expected_proteins": prep["unique_proteins"]}}
    else:
        manifest = read_json(args.root / "split_manifest.json")
        sources = {"split_manifest_sha256": sha(args.root / "split_manifest.json")}
        for split in ("train", "val", "test"):
            fasta = args.root / f"{split}.proteins.fasta"
            sources[split] = {"path": str(fasta.resolve()), "sha256": sha(fasta)}
        sources["expected_proteins"] = manifest["proteins"]

    ver = subprocess.run([args.mmseqs, "version"], check=True,
                         capture_output=True, text=True).stdout.strip()
    config = {"version": VERSION, "kind": args.kind, "nshards": args.shards,
              "sources": sources, "seed": 3407, "mmseqs_version": ver,
              "search_settings": SETTINGS}
    if prepared.exists():
        prior = read_json(prepared)
        if prior["config"] != config:
            raise RuntimeError("Existing READY.json has different data/parameters. Refusing to mix runs.")
        labels = ("full",) if args.kind == "graph" else ("train", "val", "test")
        for label in labels:
            if not Path(str(base / "db" / label) + ".dbtype").is_file():
                raise RuntimeError("Missing global target database: " + label)
            for i in range(args.shards):
                if not (base / "queries" / label / f"part_{i:04d}.faa").is_file():
                    raise RuntimeError(f"Missing query shard: {label}/{i}")
        print(f"PREPARE reused {prepared}", flush=True)
        return
    if any(base.rglob("*.done.json")):
        raise RuntimeError("Done markers exist but READY.json is missing; manual inspection required.")
    counts = {}
    labels = ("full",) if args.kind == "graph" else ("train", "val", "test")
    for label in labels:
        src = Path(sources[label]["path"])
        counts[label] = shard_fasta(src, base / "queries" / label, args.shards)
        if sources[label].get("expected_proteins") is not None:
            if counts[label]["proteins"] != sources[label]["expected_proteins"]:
                raise ValueError(f"expected protein count mismatch for {label}")
        db = base / "db" / label
        db.parent.mkdir(parents=True, exist_ok=True)
        if not Path(str(db) + ".dbtype").is_file():
            subprocess.run([args.mmseqs, "createdb", str(src), str(db),
                            "--threads", str(args.threads)], check=True)
        if not Path(str(db) + ".dbtype").is_file():
            raise RuntimeError("MMseqs target DB missing: " + str(db))
    if args.kind == "audit":
        if sum(counts[s]["proteins"] for s in labels) != sources["expected_proteins"]:
            raise RuntimeError("Audit query FASTAs do not cover all split proteins.")
    write_json(prepared, {"config": config, "counts": counts})
    print(f"PREPARE completed: {prepared}", flush=True)


def task_paths(args):
    if args.kind == "graph":
        if not 0 <= args.task < args.shards:
            raise ValueError("graph task out of bounds")
        query_label, target_label, shard, group = "full", "full", args.task, "full"
    else:
        if not 0 <= args.task < 6 * args.shards:
            raise ValueError("audit task out of bounds")
        group = PAIRS[args.task // args.shards]
        query_label, target_label = group.split("_")
        shard = args.task % args.shards
    base = base_path(args)
    query = base / "queries" / query_label / f"part_{shard:04d}.faa"
    db = base / "db" / target_label
    hit = base / "hits" / group / f"part_{shard:04d}.tsv"
    done = hit.with_suffix(".done.json")
    return query, db, hit, done, group, shard


def check_done(done: Path, hit: Path, query: Path, ready_hash: str):
    if not done.exists():
        return False
    meta = read_json(done)
    if (not hit.is_file() or meta.get("ready_sha256") != ready_hash
            or meta.get("query_sha256") != sha(query)
            or meta.get("hits_sha256") != sha(hit)):
        raise RuntimeError(f"Inconsistent done marker: {done}")
    return True


def search(args):
    base = base_path(args)
    ready = base / "READY.json"
    conf = read_json(ready)["config"]
    if conf["kind"] != args.kind or conf["nshards"] != args.shards:
        raise RuntimeError("Search and prepare configurations differ")
    query, db, hit, done, group, shard = task_paths(args)
    ready_hash = sha(ready)
    if check_done(done, hit, query, ready_hash):
        print(f"SKIP completed {group}/{shard}", flush=True)
        return
    if not Path(str(db) + ".dbtype").is_file():
        raise FileNotFoundError("Global target DB missing: " + str(db))
    hit.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.tmp_root / args.kind / group / f"part_{shard:04d}"
    tmp.parent.mkdir(parents=True, exist_ok=True)
    if tmp.exists():
        if tmp.is_symlink() or not tmp.is_dir():
            raise RuntimeError("Refusing unsafe temporary path: " + str(tmp))
        shutil.rmtree(tmp)  # only this shard's incomplete intermediate files
    partial = hit.with_suffix(".partial")
    partial.unlink(missing_ok=True)
    if query.stat().st_size == 0:
        partial.write_text(HEADER)
    else:
        cmd = [args.mmseqs, "easy-search", str(query), str(db),
               str(partial), str(tmp), "--format-output", FORMAT,
               "--min-seq-id", str(SETTINGS["min_seq_id"]),
               "-c", str(SETTINGS["min_coverage"]),
               "-s", str(SETTINGS["sensitivity"]),
               "-e", str(SETTINGS["evalue"]),
               "--max-seqs", str(SETTINGS["max_seqs"]),
               "--threads", str(args.threads)]
        print(" ".join(cmd), flush=True)
        subprocess.run(cmd, check=True)
        if not partial.exists():
            raise RuntimeError("MMseqs search output missing: " + str(partial))
        if partial.stat().st_size == 0:
            partial.write_text(HEADER)
    partial.replace(hit)
    write_json(done, {"status": "completed", "query_sha256": sha(query),
                      "ready_sha256": ready_hash, "hits_sha256": sha(hit),
                      "output": str(hit.resolve()), "job_id": os.getenv("SLURM_JOB_ID"),
                      "task": args.task})
    print(f"DONE {group}/{shard}: {hit}", flush=True)


def verify(args):
    base = base_path(args)
    ready = base / "READY.json"
    conf = read_json(ready)["config"]
    if conf["kind"] != args.kind or conf["nshards"] != args.shards:
        raise RuntimeError("Verification config differs from search")
    paths = []
    job_count = args.shards if args.kind == "graph" else 6 * args.shards
    for task in range(job_count):
        args.task = task
        query, db, hit, done, group, shard = task_paths(args)
        if not check_done(done, hit, query, sha(ready)):
            raise RuntimeError(f"INCOMPLETE shard: {group}/{shard} ({done})")
        paths.append(str(hit.resolve()))
    hits_list = base / "hits.list"
    partial = hits_list.with_suffix(".partial")
    partial.write_text("\n".join(paths) + "\n")
    partial.replace(hits_list)
    write_json(base / "COMPLETE.json", {"status": "completed",
              "kind": args.kind, "shards": args.shards, "jobs": job_count,
              "ready_sha256": sha(ready),
              "hits_list_sha256": sha(hits_list),
              "hits_list": str(hits_list.resolve())})
    print(f"VERIFIED {job_count} shards: {hits_list}", flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=("prepare", "search", "verify"))
    p.add_argument("--kind", required=True, choices=("graph", "audit"))
    p.add_argument("--root", required=True, type=Path)
    p.add_argument("--tmp-root", type=Path)
    p.add_argument("--shards", required=True, type=int)
    p.add_argument("--task", type=int)
    p.add_argument("--mmseqs", default="mmseqs")
    p.add_argument("--threads", type=int, default=8)
    args = p.parse_args()
    if not 1 <= args.shards <= 4096:
        p.error("shards must be 1..4096")
    if args.command == "search" and (args.task is None or args.tmp_root is None):
        p.error("search requires --task and --tmp-root")
    {"prepare": prepare, "search": search, "verify": verify}[args.command](args)


if __name__ == "__main__":
    main()
