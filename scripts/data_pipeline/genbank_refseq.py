"""Parse local GenBank/RefSeq transcript files into normalized records."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import Iterator

from codon import split_transcript_by_cds, validate_orf
from schema import MrnaDesignRecord, write_records


def open_maybe_gzip(path: str | Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def is_transcript_record(record_id: str, molecule_type: str) -> bool:
    """Heuristic for records whose sequence is already an RNA/transcript."""

    record_id = record_id.upper()
    molecule_type = molecule_type.lower()
    return (
        "mrna" in molecule_type
        or "rna" in molecule_type
        or record_id.startswith(("NM_", "XM_", "NR_", "XR_"))
    )


def iter_genbank_records(path: str | Path) -> Iterator[MrnaDesignRecord]:
    """Yield records from a GenBank file.

    Biopython is imported lazily so the rest of the pipeline can run without it.
    """

    try:
        from Bio import SeqIO
    except ImportError as exc:
        raise SystemExit("Biopython is required: pip install biopython") from exc

    try:
        from Bio.Seq import UndefinedSequenceError
    except ImportError:
        UndefinedSequenceError = ValueError

    skipped_undefined = 0
    skipped_invalid_orf = 0

    with open_maybe_gzip(path) as handle:
        for record in SeqIO.parse(handle, "genbank"):
            try:
                transcript_sequence = str(record.seq)
            except UndefinedSequenceError:
                skipped_undefined += 1
                continue
            organism = record.annotations.get("organism", "")
            molecule_type = str(record.annotations.get("molecule_type", ""))
            is_transcript = is_transcript_record(record.id, molecule_type)
            taxonomy_id = ""
            for feature in record.features:
                for value in feature.qualifiers.get("db_xref", []):
                    if value.startswith("taxon:"):
                        taxonomy_id = value.split(":", 1)[1]
                        break
                if taxonomy_id:
                    break

            for feature_index, feature in enumerate(record.features):
                if feature.type != "CDS":
                    continue
                protein_sequence = feature.qualifiers.get("translation", [""])[0]
                if not protein_sequence:
                    continue

                cds_start = int(feature.location.start)
                cds_end = int(feature.location.end)
                orf = str(feature.extract(record.seq))
                utr5 = ""
                utr3 = ""
                if is_transcript:
                    utr5, _, utr3 = split_transcript_by_cds(
                        transcript_sequence,
                        cds_start,
                        cds_end,
                    )
                validation = validate_orf(protein_sequence, orf, require_stop=False)
                if not validation.is_valid:
                    validation = validate_orf(protein_sequence, orf, require_stop=True)
                if not validation.is_valid:
                    skipped_invalid_orf += 1
                    continue

                protein_id = feature.qualifiers.get("protein_id", [""])[0]
                transcript_id = record.id
                gene_name = feature.qualifiers.get("gene", [""])[0]
                record_id = f"{transcript_id}:cds:{feature_index}"
                db_xrefs = feature.qualifiers.get("db_xref", [])

                yield MrnaDesignRecord(
                    record_id=record_id,
                    source_db="",
                    source_accession=record.id,
                    transcript_id=transcript_id,
                    protein_id=protein_id,
                    gene_name=gene_name,
                    taxonomy_id=taxonomy_id,
                    organism=organism,
                    protein_sequence=protein_sequence.strip().upper().replace("*", ""),
                    orf_sequence=orf,
                    utr5_sequence=utr5,
                    utr3_sequence=utr3,
                    metadata={
                        "definition": record.description,
                        "molecule_type": molecule_type,
                        "is_transcript_record": is_transcript,
                        "cds_start": cds_start,
                        "cds_end": cds_end,
                        "db_xref": ";".join(db_xrefs),
                        "validation_reason": validation.reason,
                    },
                )

    if skipped_undefined:
        print(f"Skipped {skipped_undefined} records with undefined sequence content.")
    if skipped_invalid_orf:
        print(f"Skipped {skipped_invalid_orf} CDS features with invalid ORF/protein translation.")


def parse_genbank(args: argparse.Namespace) -> int:
    def records_with_source():
        count = 0
        for record in iter_genbank_records(args.input):
            record.source_db = args.source_db
            yield record
            count += 1
            if args.limit and count >= args.limit:
                break

    written = write_records(args.output, records_with_source())
    print(f"Wrote {written} records to {args.output}.")
    return written


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Local .gbff/.gbff.gz file.")
    parser.add_argument("--output", required=True, help="Output normalized TSV.")
    parser.add_argument("--source-db", default="refseq", help="Source database label.")
    parser.add_argument("--limit", type=int, default=0, help="Optional max records.")
    return parser


def main() -> None:
    parse_genbank(build_parser().parse_args())


if __name__ == "__main__":
    main()
