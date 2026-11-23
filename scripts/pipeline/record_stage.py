#!/usr/bin/env python3
"""Store exact input/output paths + hashes on successful Nibi pipeline steps."""
import argparse
import hashlib
import json
import os
import socket
import subprocess
from datetime import datetime, timezone
from pathlib import Path

def entry(raw):
    label, sep, filename = raw.partition("=")
    if not sep or not label or not filename.startswith("/"):
        raise ValueError("expected LABEL=/absolute/path")
    path = Path(filename)
    if not path.is_file():
        raise FileNotFoundError(path)
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return label, {"path": str(path), "resolved_path": str(path.resolve()),
                   "bytes": path.stat().st_size, "sha256": digest.hexdigest()}

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--stage", required=True)
    p.add_argument("--record", required=True, type=Path)
    p.add_argument("--input", action="append", required=True)
    p.add_argument("--output", action="append", required=True)
    args = p.parse_args()
    inputs = dict(entry(x) for x in args.input)
    outputs = dict(entry(x) for x in args.output)
    if len(inputs) != len(args.input) or len(outputs) != len(args.output):
        raise ValueError("duplicate artifact label")
    try:
        git_sha = subprocess.check_output(["git", "rev-parse", "HEAD"],
                                          text=True, stderr=subprocess.DEVNULL).strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        git_sha = None
    report = {
        "stage": args.stage, "status": "completed", "seed": 3407,
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_sha, "host": socket.gethostname(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "slurm_account": os.environ.get("SLURM_JOB_ACCOUNT"),
        "slurm_cpus_per_task": os.environ.get("SLURM_CPUS_PER_TASK"),
        "inputs": inputs, "outputs": outputs,
    }
    args.record.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.record.with_name(args.record.name + "." + str(os.getpid()) + ".tmp")
    try:
        temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        os.replace(temporary, args.record)
    finally:
        temporary.unlink(missing_ok=True)
    print("Wrote input/output manifest:", args.record)

if __name__ == "__main__":
    main()
