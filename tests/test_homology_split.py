"""Small synthetic homology/leakage regression suite (stdlib only)."""
import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "data_pipeline" / "homology_split.py"
FIELDS = ["record_id", "protein_sequence", "orf_sequence", "utr5_sequence", "utr3_sequence"]
HITS = ["query", "target", "pident", "alnlen", "qlen", "tlen", "qcov", "tcov"]


def run(*args, succeed=True):
    completed = subprocess.run([sys.executable, str(SCRIPT), *map(str, args)],
                               capture_output=True, text=True)
    if succeed and completed.returncode:
        raise AssertionError(completed.stdout + "\n" + completed.stderr)
    return completed


class HomologySplitTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / "merged.tsv"
        self.work = self.root / "work"
        self.work.mkdir()

    def rows(self, proteins):
        with self.source.open("w", newline="") as f:
            writer = csv.DictWriter(f, delimiter="\t", fieldnames=FIELDS)
            writer.writeheader()
            for n, protein in enumerate(proteins):
                writer.writerow(dict(record_id=f"r{n}", protein_sequence=protein,
                                     orf_sequence="ATG" * len(protein),
                                     utr5_sequence="AC", utr3_sequence="GT"))

    def ids(self):
        records = (self.work / "unique_proteins.fasta").read_text().splitlines()
        return [line[1:] for line in records if line.startswith(">")]

    def hits(self, edges):
        path = self.work / "hits.tsv"
        with path.open("w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(HITS)
            for edge in edges:
                writer.writerow(edge)
        return path

    def split(self, hits, max_dev=0.5):
        return run("split", "--input", self.source,
                   "--prepare-manifest", self.work / "prepare.json",
                   "--hits", hits, "--output-dir", self.work,
                   "--max-fraction-deviation", max_dev)

    def test_connected_components_and_repeatability(self):
        self.rows(["M" + "A"*119, "M" + "C"*119, "M" + "D"*119,
                   "M" + "E"*119, "M" + "F"*119, "M" + "G"*119,
                   "M" + "A"*119])
        run("prepare", "--input", self.source, "--output-dir", self.work)
        a,b,c,d,e,f = self.ids()
        # A--B--C must transitively form one cluster.
        hits = self.hits([[a,b,45,120,120,120,1,1],
                          [b,c,45,120,120,120,1,1],
                          [d,e,10,120,120,120,1,1],
                          [e,f,98,20,120,120,.167,.167]])
        self.split(hits)
        assn = {}
        with (self.work / "protein_assignments.tsv").open() as f:
            for row in csv.DictReader(f, delimiter="\t"):
                assn[row["protein_sha256"]] = row
        self.assertEqual(assn[a]["protein_cluster_id"], assn[b]["protein_cluster_id"])
        self.assertEqual(assn[a]["protein_cluster_id"], assn[c]["protein_cluster_id"])
        self.assertEqual(assn[a]["split"], assn[c]["split"])
        self.assertEqual(len(assn), 6)
        self.assertEqual(assn[a]["pair_count"], "2")
        first = (self.work / "protein_assignments.tsv").read_bytes()
        self.split(hits)
        self.assertEqual(first, (self.work / "protein_assignments.tsv").read_bytes())

    def test_near_duplicate_fragment_guard_and_audit_failure(self):
        self.rows(["M" + "A"*119, "M" + "C"*119, "M" + "D"*119,
                   "M" + "E"*119, "M" + "F"*119])
        run("prepare", "--input", self.source, "--output-dir", self.work)
        a,b,c,d,e = self.ids()
        hits = self.hits([[a,b,95,100,120,120,.83,.83],
                          [c,d,95,90,120,120,.75,.75],
                          [d,e,91,90,120,120,.75,.75]])
        self.split(hits)
        assn = {row["protein_sha256"]: row
                for row in csv.DictReader((self.work/"protein_assignments.tsv").open(),
                                          delimiter="\t")}
        self.assertEqual(assn[a]["protein_cluster_id"], assn[b]["protein_cluster_id"])
        self.assertNotEqual(assn[c]["protein_cluster_id"], assn[d]["protein_cluster_id"])
        audit_hits = self.hits([[a,c,95,110,120,120,.92,.92]])
        result = run("audit", "--train", self.work/"train.tsv", "--val", self.work/"val.tsv",
                     "--test", self.work/"test.tsv", "--hits", audit_hits,
                     "--report", self.work/"audit.json", succeed=False)
        # If the extra edge happens to lie within a split, audit is allowed to pass.
        self.assertEqual(result.returncode != 0, assn[a]["split"] != assn[c]["split"])
        self.assertEqual(json.loads((self.work/"audit.json").read_text())["passed"],
                         assn[a]["split"] == assn[c]["split"])

    def test_changed_source_is_rejected(self):
        self.rows(["MAAA", "MCCC", "MDDD"])
        run("prepare", "--input", self.source, "--output-dir", self.work)
        a,b,c = self.ids()
        hits = self.hits([[a,b,60,4,4,4,1,1]])
        with self.source.open("a") as f:
            f.write("extra\tMEEE\tATG\tA\tG\n")
        result = run("split", "--input", self.source,
                     "--prepare-manifest", self.work/"prepare.json",
                     "--hits", hits, "--output-dir", self.work, succeed=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("differs from prepare.json", result.stderr)


if __name__ == "__main__":
    unittest.main()
