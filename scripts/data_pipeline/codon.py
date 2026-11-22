"""Codon utilities for protein-ORF-UTR data processing."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable


CODON_TO_AA = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

AA_TO_CODONS = {}
for codon, aa in CODON_TO_AA.items():
    AA_TO_CODONS.setdefault(aa, []).append(codon)

STOP_CODONS = set(AA_TO_CODONS["*"])
DNA_BASES = set("ACGT")


@dataclass(frozen=True)
class OrfValidation:
    """Validation result for a coding sequence and protein pair."""

    is_valid: bool
    reason: str = ""
    translated: str = ""


def normalize_dna(sequence: str | None) -> str:
    """Return upper-case DNA sequence, converting U to T and removing whitespace."""

    if sequence is None:
        return ""
    return "".join(str(sequence).upper().replace("U", "T").split())


def split_codons(sequence: str) -> list[str]:
    sequence = normalize_dna(sequence)
    return [sequence[i : i + 3] for i in range(0, len(sequence), 3)]


def translate_orf(orf_sequence: str, include_stop: bool = False) -> str:
    """Translate a DNA/RNA ORF using the standard genetic code."""

    codons = split_codons(orf_sequence)
    aas: list[str] = []
    for codon in codons:
        if len(codon) != 3:
            aas.append("X")
            continue
        aa = CODON_TO_AA.get(codon, "X")
        if aa == "*" and not include_stop:
            break
        aas.append(aa)
    return "".join(aas)


def validate_orf(
    protein_sequence: str | None,
    orf_sequence: str | None,
    require_stop: bool = True,
) -> OrfValidation:
    """Validate that an ORF encodes the provided protein sequence."""

    protein = (protein_sequence or "").strip().upper().replace("*", "")
    orf = normalize_dna(orf_sequence)

    if not protein:
        return OrfValidation(False, "missing_protein")
    if not orf:
        return OrfValidation(False, "missing_orf")
    if set(orf) - DNA_BASES:
        return OrfValidation(False, "invalid_nucleotide")
    if len(orf) % 3 != 0:
        return OrfValidation(False, "orf_length_not_multiple_of_3")

    codons = split_codons(orf)
    if require_stop:
        expected_len = (len(protein) + 1) * 3
        if len(orf) != expected_len:
            return OrfValidation(False, "orf_length_does_not_match_protein_plus_stop")
        if codons[-1] not in STOP_CODONS:
            return OrfValidation(False, "missing_terminal_stop")
        translated = "".join(CODON_TO_AA.get(c, "X") for c in codons[:-1])
    else:
        expected_len = len(protein) * 3
        if len(orf) != expected_len:
            return OrfValidation(False, "orf_length_does_not_match_protein")
        translated = "".join(CODON_TO_AA.get(c, "X") for c in codons)

    if "X" in translated:
        return OrfValidation(False, "unknown_codon", translated)
    if translated != protein:
        return OrfValidation(False, "translation_mismatch", translated)

    return OrfValidation(True, "valid", translated)


def split_transcript_by_cds(
    transcript_sequence: str,
    cds_start: int,
    cds_end: int,
    one_based: bool = False,
) -> tuple[str, str, str]:
    """Split transcript into 5'UTR, ORF/CDS, and 3'UTR.

    `cds_start` is inclusive and `cds_end` is exclusive when `one_based=False`.
    With `one_based=True`, coordinates are converted from GenBank-style
    one-based inclusive coordinates.
    """

    transcript = normalize_dna(transcript_sequence)
    start = int(cds_start)
    end = int(cds_end)
    if one_based:
        start -= 1
    start = max(0, start)
    end = min(len(transcript), end)
    if end < start:
        raise ValueError(f"Invalid CDS interval: start={cds_start}, end={cds_end}")
    return transcript[:start], transcript[start:end], transcript[end:]


def legal_codons_for_amino_acid(amino_acid: str) -> list[str]:
    """Return legal synonymous codons for an amino acid."""

    return list(AA_TO_CODONS.get(amino_acid.upper(), []))


def gc_content(sequence: str) -> float:
    seq = normalize_dna(sequence)
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def count_motifs(sequence: str, motifs: Iterable[str]) -> dict[str, int]:
    seq = normalize_dna(sequence)
    return {motif.upper(): seq.count(motif.upper()) for motif in motifs}

