"""Synthetic state progression for safe, quota-limited stage submissions."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

script = Path(__file__).resolve().parents[1] / "scripts/pipeline/next_stage.py"
spec = importlib.util.spec_from_file_location("next_stage", script)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class NextStageTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)

    def touch(self, name, value="x"):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(value)

    def stage(self):
        return mod.next_step(self.root)[0]

    def test_requires_all_graph_and_audit_shards(self):
        self.assertEqual(self.stage(), "protein")
        self.touch("unique_proteins.fasta")
        self.touch("prepare.json", "{}")
        self.assertEqual(self.stage(), "graph_prepare")
        self.touch("graph_shards_v1/READY.json", "{}")
        self.assertEqual(self.stage(), "graph_search")
        for i in range(256):
            self.touch(f"graph_shards_v1/hits/full/part_{i:04d}.done.json", "{}")
        self.assertEqual(self.stage(), "graph_verify")
        self.touch("graph_shards_v1/COMPLETE.json", '{"status":"completed"}')
        self.assertEqual(self.stage(), "split")
        for file in (
            "split_manifest.json", "train.tsv", "val.tsv", "test.tsv",
            "clusters.tsv", "protein_assignments.tsv",
            "train.proteins.fasta", "val.proteins.fasta", "test.proteins.fasta"
        ):
            self.touch(file, "{}")
        self.assertEqual(self.stage(), "audit_prepare")
        self.touch("audit_shards_v1/READY.json", "{}")
        self.assertEqual(self.stage(), "audit_search")
        for pair in ("train_val", "train_test", "val_train", "val_test",
                     "test_train", "test_val"):
            for i in range(256):
                self.touch(f"audit_shards_v1/hits/{pair}/part_{i:04d}.done.json", "{}")
        self.assertEqual(self.stage(), "audit_verify")
        self.touch("audit_shards_v1/COMPLETE.json", '{"status":"completed"}')
        self.assertEqual(self.stage(), "audit")
        self.touch("leakage_audit.json", '{"passed":true}')
        self.assertEqual(self.stage(), "DONE")
        self.touch("leakage_audit.json", '{"passed":false}')
        self.assertEqual(self.stage(), "AUDIT_FAILED")


if __name__ == "__main__":
    unittest.main()
