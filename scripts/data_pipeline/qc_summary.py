"""Summarize normalized mRNA design TSV files."""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import Counter
from pathlib import Path

from codon import gc_content
from schema import raise_csv_field_limit


def safe_len(value: str | None) -> int:
    return len(value or "")


def summarize_lengths(values: list[int]) -> str:
    if not values:
        return "n=0"
    values_sorted = sorted(values)
    return (
        f"n={len(values)} "
        f"min={values_sorted[0]} "
        f"p50={int(statistics.median(values_sorted))} "
        f"mean={statistics.mean(values_sorted):.1f} "
        f"max={values_sorted[-1]}"
    )


def summarize_file(path: Path) -> dict[str, object]:
    raise_csv_field_limit()
    rows = 0
    protein_lengths = []
    orf_lengths = []
    utr5_lengths = []
    utr3_lengths = []
    source_counts = Counter()
    organism_counts = Counter()
    proteins = set()
    orfs = set()
    protein_orf_pairs = set()
    gc_sum = 0.0
    gc_count = 0

    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows += 1
            protein = row.get("protein_sequence", "")
            orf = row.get("orf_sequence", "")
            utr5 = row.get("utr5_sequence", "")
            utr3 = row.get("utr3_sequence", "")
            protein_lengths.append(safe_len(protein))
            orf_lengths.append(safe_len(orf))
            if utr5:
                utr5_lengths.append(safe_len(utr5))
            if utr3:
                utr3_lengths.append(safe_len(utr3))
            source_counts.update([row.get("source_db", "")])
            if row.get("organism"):
                organism_counts.update([row.get("organism", "")])
            proteins.add(protein)
            orfs.add(orf)
            protein_orf_pairs.add((protein, orf))
            if orf:
                gc_sum += gc_content(orf)
                gc_count += 1

    return {
        "path": str(path),
        "rows": rows,
        "protein_unique": len(proteins),
        "orf_unique": len(orfs),
        "protein_orf_unique": len(protein_orf_pairs),
        "with_utr5": len(utr5_lengths),
        "with_utr3": len(utr3_lengths),
        "protein_lengths": summarize_lengths(protein_lengths),
        "orf_lengths": summarize_lengths(orf_lengths),
        "utr5_lengths": summarize_lengths(utr5_lengths),
        "utr3_lengths": summarize_lengths(utr3_lengths),
        "orf_gc_mean": gc_sum / gc_count if gc_count else 0.0,
        "source_counts": source_counts,
        "top_organisms": organism_counts.most_common(10),
    }


def print_summary(summary: dict[str, object]) -> None:
    print(f"\n== {summary['path']} ==")
    print(f"rows: {summary['rows']}")
    print(f"unique proteins: {summary['protein_unique']}")
    print(f"unique ORFs: {summary['orf_unique']}")
    print(f"unique protein+ORF pairs: {summary['protein_orf_unique']}")
    print(f"with 5UTR: {summary['with_utr5']}")
    print(f"with 3UTR: {summary['with_utr3']}")
    print(f"protein lengths: {summary['protein_lengths']}")
    print(f"ORF lengths: {summary['orf_lengths']}")
    print(f"5UTR lengths: {summary['utr5_lengths']}")
    print(f"3UTR lengths: {summary['utr3_lengths']}")
    print(f"mean ORF GC: {summary['orf_gc_mean']:.4f}")
    print(f"sources: {dict(summary['source_counts'])}")
    print("top organisms:")
    for organism, count in summary["top_organisms"]:
        print(f"  {count}\t{organism}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", help="Normalized TSV files or directories.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    for input_path in args.inputs:
        path = Path(input_path)
        if not path.exists():
            print(f"\n== {path} ==")
            print("warning: file does not exist; skipped")
            continue
        if path.is_dir():
            files = sorted(path.rglob("*.tsv"))
            if not files:
                print(f"\n== {path} ==")
                print("warning: directory contains no TSV files; skipped")
                continue
            for file_path in files:
                print_summary(summarize_file(file_path))
            continue
        print_summary(summarize_file(path))


if __name__ == "__main__":
    main()
