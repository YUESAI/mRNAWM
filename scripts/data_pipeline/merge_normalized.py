"""Merge normalized TSV files from files or directories with exact deduplication."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from schema import MrnaDesignRecord, NORMALIZED_COLUMNS, iter_records


def stable_hash(values: list[str]) -> str:
    import hashlib

    digest = hashlib.sha256()
    for value in values:
        digest.update(value.encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest()


def dedup_key(record: MrnaDesignRecord, mode: str) -> str:
    if mode == "protein":
        return stable_hash([record.protein_sequence])
    if mode == "orf":
        return stable_hash([record.orf_sequence])
    if mode == "protein_orf":
        return stable_hash([record.protein_sequence, record.orf_sequence])
    if mode == "accession":
        return stable_hash([record.source_db, record.transcript_id, record.protein_id])
    raise ValueError(f"Unknown dedup mode: {mode}")


def discover_tsv_inputs(paths: list[str]) -> list[Path]:
    inputs: list[Path] = []
    for raw_path in paths:
        path = Path(raw_path)
        if not path.exists():
            print(f"warning: input path does not exist; skipped: {path}")
            continue
        if path.is_dir():
            inputs.extend(sorted(path.rglob("*.tsv")))
        else:
            inputs.append(path)
    return inputs


def merge_normalized(args: argparse.Namespace) -> int:
    inputs = discover_tsv_inputs(args.inputs)
    if not inputs:
        raise SystemExit("No TSV inputs found.")

    seen_keys: set[str] = set()
    total = 0
    duplicates = 0
    empty_files = 0
    written = 0

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=NORMALIZED_COLUMNS, delimiter="\t")
        writer.writeheader()

        for input_path in inputs:
            file_records = 0
            for record in iter_records(input_path):
                file_records += 1
                total += 1
                key = dedup_key(record, args.mode)
                if key in seen_keys:
                    duplicates += 1
                    continue
                seen_keys.add(key)
                writer.writerow(record.to_row())
                written += 1
            if file_records == 0:
                empty_files += 1

    print(f"Read {total} records from {len(inputs)} files.")
    print(f"Skipped {duplicates} duplicate records by {args.mode}.")
    print(f"Skipped {empty_files} empty TSV files.")
    print(f"Wrote {written} unique records to {args.output}.")
    return written


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="Normalized TSV files or directories containing TSV files.",
    )
    parser.add_argument("--output", required=True, help="Output normalized TSV.")
    parser.add_argument(
        "--mode",
        choices=["protein", "orf", "protein_orf", "accession"],
        default="protein_orf",
        help="Exact deduplication key.",
    )
    return parser


def main() -> None:
    merge_normalized(build_parser().parse_args())


if __name__ == "__main__":
    main()
