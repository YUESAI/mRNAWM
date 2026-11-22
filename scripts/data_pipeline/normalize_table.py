"""Normalize existing CSV/TSV protein-mRNA tables into the new schema."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from codon import split_transcript_by_cds, validate_orf
from schema import raise_csv_field_limit
from schema import MrnaDesignRecord, write_records


DEFAULT_COLUMN_ALIASES = {
    "protein_sequence": ["Protein Sequence", "Protein_Sequence", "protein_sequence"],
    "mrna_sequence": ["mRNA Sequence", "mRNA_Sequence", "cDNA_Sequence", "transcript_sequence"],
    "orf_sequence": ["ORF Sequence", "ORF_Sequence", "CDS Sequence", "CDS_Sequence", "orf_sequence"],
    "utr5_sequence": ["5' UTR", "UTR_5", "UTR5", "utr5_sequence"],
    "utr3_sequence": ["3' UTR", "UTR_3", "UTR3", "utr3_sequence"],
    "cds_start": ["CDS Start", "CDS_Start", "cds_start"],
    "cds_end": ["CDS End", "CDS_End", "cds_end"],
    "codon_start": ["Codon Start", "Codon_Start", "codon_start"],
    "source_accession": ["Seq ID", "mRNA ID", "Transcript ID", "transcript_id"],
    "transcript_id": ["Transcript ID", "Transcript_ID", "mRNA ID", "transcript_id"],
    "protein_id": ["Protein ID", "Protein_ID", "protein_id"],
    "gene_name": ["Gene Name", "Gene_Name", "gene_name"],
    "organism": ["Organism", "organism"],
    "taxonomy_id": ["Taxonomy ID", "Taxonomy_ID", "taxid", "taxonomy_id"],
    "is_antigen": ["Antigen Protein", "Is_Antigen", "is_antigen"],
    "disease": ["Disease", "disease"],
    "pathogen": ["Pathogen", "pathogen"],
}


def pick(row: dict[str, str], field: str) -> str:
    for column in DEFAULT_COLUMN_ALIASES[field]:
        if column in row and row[column] not in (None, ""):
            return str(row[column])
    return ""


def sniff_dialect(path: Path) -> csv.Dialect:
    with path.open("r", newline="", encoding="utf-8") as handle:
        sample = handle.read(4096)
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t")
    except csv.Error:
        return csv.excel


def parse_int(value: str, default: int = 0) -> int:
    value = str(value).strip().replace("<", "").replace(">", "")
    if value == "":
        return default
    return int(float(value))


def normalize_row(
    row: dict[str, str],
    row_index: int,
    source_db: str,
    require_stop: bool,
) -> MrnaDesignRecord | None:
    protein_sequence = pick(row, "protein_sequence").strip().upper().replace("*", "")
    mrna_sequence = pick(row, "mrna_sequence")
    orf_sequence = pick(row, "orf_sequence")
    utr5_sequence = pick(row, "utr5_sequence")
    utr3_sequence = pick(row, "utr3_sequence")

    if not protein_sequence:
        return None

    if not orf_sequence and mrna_sequence and pick(row, "cds_start") and pick(row, "cds_end"):
        cds_start = parse_int(pick(row, "cds_start"))
        cds_end = parse_int(pick(row, "cds_end"))
        codon_start = parse_int(pick(row, "codon_start"), default=1)
        utr5_sequence, orf_sequence, utr3_sequence = split_transcript_by_cds(
            mrna_sequence,
            cds_start + codon_start - 1,
            cds_end + codon_start - 1,
        )

    validation = validate_orf(protein_sequence, orf_sequence, require_stop=require_stop)
    if not validation.is_valid:
        return None

    source_accession = pick(row, "source_accession")
    transcript_id = pick(row, "transcript_id") or source_accession
    protein_id = pick(row, "protein_id")
    record_id = f"{source_db}:{transcript_id or source_accession or row_index}"

    return MrnaDesignRecord(
        record_id=record_id,
        source_db=source_db,
        source_accession=source_accession,
        transcript_id=transcript_id,
        protein_id=protein_id,
        gene_name=pick(row, "gene_name"),
        taxonomy_id=pick(row, "taxonomy_id"),
        organism=pick(row, "organism"),
        protein_sequence=protein_sequence,
        orf_sequence=orf_sequence,
        utr5_sequence=utr5_sequence,
        utr3_sequence=utr3_sequence,
        is_antigen=pick(row, "is_antigen"),
        disease=pick(row, "disease"),
        pathogen=pick(row, "pathogen"),
        metadata={
            "source_row": row_index,
            "validation_reason": validation.reason,
        },
    )


def normalize_table(args: argparse.Namespace) -> int:
    raise_csv_field_limit()
    input_path = Path(args.input)
    dialect = sniff_dialect(input_path)
    records: list[MrnaDesignRecord] = []
    skipped = 0

    with input_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, dialect=dialect)
        for row_index, row in enumerate(reader, start=1):
            record = normalize_row(row, row_index, args.source_db, args.require_stop)
            if record is None:
                skipped += 1
                continue
            records.append(record)

    written = write_records(args.output, records)
    print(f"Wrote {written} records to {args.output}; skipped {skipped} rows.")
    return written


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input CSV/TSV table.")
    parser.add_argument("--output", required=True, help="Output normalized TSV.")
    parser.add_argument("--source-db", required=True, help="Source database name.")
    parser.add_argument(
        "--no-require-stop",
        dest="require_stop",
        action="store_false",
        help="Allow ORF sequences without a terminal stop codon.",
    )
    parser.set_defaults(require_stop=True)
    return parser


def main() -> None:
    normalize_table(build_parser().parse_args())


if __name__ == "__main__":
    main()
