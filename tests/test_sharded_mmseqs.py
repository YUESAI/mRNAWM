"""Regression tests: global target sharding, resumable markers, strict coverage."""
import hashlib
import importlib.util
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "data_pipeline" / "sharded_mmseqs.py"
spec = importlib.util.spec_from_file_location("sharded_mmseqs", SCRIPT)
sh = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sh)


class Args:
    pass


class ShardingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.root.joinpath("prepare.json").write_text("{}")
        self.proteins = []
        fasta = self.root / "unique_proteins.fasta"
        with fasta.open("w") as out:
            for i in range(19):
                seq = "M" + chr(ord("A") + i) * 100
                pid = hashlib.sha256(seq.encode()).hexdigest()
                self.proteins.append((pid, seq))
                out.write(f">{pid}\n{seq}\n")
        (self.root / "prepare.json").write_text(json.dumps({
            "fasta_sha256": sh.sha(fasta), "unique_proteins": len(self.proteins)
        }))
        self.args = Args()
        self.args.kind = "graph"
        self.args.root = self.root
        self.args.shards = 4
        self.args.mmseqs = "mock-mmseqs"
        self.args.threads = 2
        self.args.tmp_root = self.root / "tmp"

    def fake_run(self, cmd, **kwargs):
        if cmd[1] == "version":
            from subprocess import CompletedProcess
            return CompletedProcess(cmd, 0, stdout="version-test\n")
        if cmd[1] == "createdb":
            Path(cmd[3] + ".dbtype").write_text("db")
            from subprocess import CompletedProcess
            return CompletedProcess(cmd, 0)
        if cmd[1] == "easy-search":
            self.assertEqual(cmd[3], str(self.root / "graph_shards_v1/db/full"))
            Path(cmd[4]).write_text(sh.HEADER)
            from subprocess import CompletedProcess
            return CompletedProcess(cmd, 0)
        raise AssertionError(cmd)

    def test_all_shards_search_same_full_target_and_resume(self):
        with patch.object(sh.subprocess, "run", side_effect=self.fake_run):
            sh.prepare(self.args)
            base = self.root / "graph_shards_v1"
            ids = []
            for i in range(4):
                ids += [name for name, _ in sh.fasta_iter(base / f"queries/full/part_{i:04d}.faa")]
            self.assertEqual(set(ids), set(p for p, _ in self.proteins))
            self.assertEqual(len(ids), len(set(ids)))
            for i in range(4):
                self.args.task = i
                sh.search(self.args)
                sh.search(self.args)  # skip completed shard
            sh.verify(self.args)
            self.assertEqual(len((base / "hits.list").read_text().splitlines()), 4)
            self.assertEqual(sh.read_json(base / "COMPLETE.json")["jobs"], 4)
            # A damaged output must never be silently reused.
            hit = base / "hits/full/part_0000.tsv"
            hit.write_text("tampered\n")
            with self.assertRaises(RuntimeError):
                sh.verify(self.args)

    def test_missing_shard_blocks_global_merge(self):
        with patch.object(sh.subprocess, "run", side_effect=self.fake_run):
            sh.prepare(self.args)
            self.args.task = 0
            sh.search(self.args)
            with self.assertRaisesRegex(RuntimeError, "INCOMPLETE shard"):
                sh.verify(self.args)

    def test_audit_tasks_cover_all_six_directions(self):
        self.args.kind = "audit"
        self.args.shards = 3
        observed = []
        for i in range(18):
            self.args.task = i
            q, db, _, _, label, shard = sh.task_paths(self.args)
            self.assertNotEqual(label.split("_")[0], label.split("_")[1])
            self.assertTrue(str(db).endswith("/db/" + label.split("_")[1]))
            observed.append((label, shard))
        self.assertEqual(set(observed), {(p, i) for p in sh.PAIRS for i in range(3)})


if __name__ == "__main__":
    unittest.main()
