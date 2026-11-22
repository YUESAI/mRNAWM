"""Normalized record schema for mRNA design datasets."""

from __future__ import annotations

import csv
import json
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable


NORMALIZED_COLUMNS = [
    "record_id",
    "source_db",
    "source_accession",
    "transcript_id",
    "protein_id",
    "gene_name",
    "taxonomy_id",
    "organism",
    "protein_sequence",
    "orf_sequence",
    "utr5_sequence",
    "utr3_sequence",
    "is_antigen",
    "disease",
    "pathogen",
    "metadata_json",
]

RECORD_COLUMNS = [column for column in NORMALIZED_COLUMNS if column != "metadata_json"]


def raise_csv_field_limit() -> None:
    """Allow very long biological sequence fields in TSV/CSV files."""

    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit = int(limit / 10)


@dataclass
class MrnaDesignRecord:
    """A database-independent protein-ORF-UTR training record."""

    record_id: str
    source_db: str
    source_accession: str = ""
    transcript_id: str = ""
    protein_id: str = ""
    gene_name: str = ""
    taxonomy_id: str = ""
    organism: str = ""
    protein_sequence: str = ""
    orf_sequence: str = ""
    utr5_sequence: str = ""
    utr3_sequence: str = ""
    is_antigen: str = ""
    disease: str = ""
    pathogen: str = ""
    metadata: dict[str, str | int | float | bool | None] = field(default_factory=dict)

    def to_row(self) -> dict[str, str]:
        raw = asdict(self)
        metadata = raw.pop("metadata")
        row = {column: "" for column in NORMALIZED_COLUMNS}
        for key, value in raw.items():
            row[key] = "" if value is None else str(value)
        row["metadata_json"] = json.dumps(metadata, sort_keys=True, ensure_ascii=True)
        return row


def write_records(path: str | Path, records: Iterable[MrnaDesignRecord]) -> int:
    """Write normalized records as TSV and return the number of written rows."""

    raise_csv_field_limit()
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    count = 0
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=NORMALIZED_COLUMNS, delimiter="\t")
        writer.writeheader()
        for record in records:
            writer.writerow(record.to_row())
            count += 1
    return count


def read_records(path: str | Path) -> list[MrnaDesignRecord]:
    """Read normalized TSV records."""

    return list(iter_records(path))


def iter_records(path: str | Path):
    """Iterate normalized TSV records without loading the full file."""

    raise_csv_field_limit()
    with Path(path).open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            metadata_json = row.pop("metadata_json", "{}") or "{}"
            metadata = json.loads(metadata_json)
            record_kwargs = {column: row.get(column, "") for column in RECORD_COLUMNS}
            yield MrnaDesignRecord(metadata=metadata, **record_kwargs)
