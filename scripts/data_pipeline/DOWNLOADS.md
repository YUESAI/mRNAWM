# Download and Processing Commands

This project should rebuild data from structured sources rather than keyword
search. Use taxonomy IDs, curated accessions, and source metadata whenever
possible.

## Recommended First-Pass Sources

### 1. RefSeq Human Transcript GBFF

Purpose:

- high-quality host transcript records
- protein-CDS pairs
- 5'UTR and 3'UTR extraction from transcript records
- baseline human codon/UTR distribution

Source:

- NCBI RefSeq human annotation FTP
- transcript file: `GRCh38_latest_rna.gbff.gz`

### 2. Broad RefSeq Vertebrate Mammalian RNA GBFF

Purpose:

- broad mammalian protein-CDS-UTR records
- improves general codon/UTR pretraining without relying only on human genes
- avoids overfitting the general model to one or two organisms

Source:

- NCBI RefSeq release FTP
- release-level files: `vertebrate_mammalian.*.rna.gbff.gz`

### 3. RefSeq Viral GenBank Records

Purpose:

- pathogen/viral protein-ORF records
- closer to vaccine antigen inputs

Note:

- viral/genomic records usually do not provide host-style 5'UTR/3'UTR.
- The parser extracts CDS/ORF with `feature.extract()` and leaves UTR blank
  unless the source record is a transcript record.

### 4. UniProt / IEDB Annotation

Purpose:

- antigen, disease, pathogen, epitope metadata
- annotation join onto protein sequences or accessions

These are annotation sources, not primary mRNA/ORF sources.

## Install Tools

```bash
python3 -m pip install biopython
```

Optional but useful later:

```bash
python3 -m pip install pandas
```

## Download Commands

Set the server-side data root and create folders:

```bash
export DATA_ROOT=/home/UWO/ylong66/data/mRNA_pr_wm

mkdir -p ${DATA_ROOT}/raw/refseq_human ${DATA_ROOT}/raw/refseq_viral
mkdir -p ${DATA_ROOT}/normalized ${DATA_ROOT}/splits/exact_protein
```

Download RefSeq human transcript records:

```bash
curl -fL \
  https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.gbff.gz \
  -o ${DATA_ROOT}/raw/refseq_human/GRCh38_latest_rna.gbff.gz
gzip -t ${DATA_ROOT}/raw/refseq_human/GRCh38_latest_rna.gbff.gz
```

Optional quick viral test download:

```bash
curl -fL \
  https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.genomic.gbff.gz \
  -o ${DATA_ROOT}/raw/refseq_viral/viral.1.genomic.gbff.gz
gzip -t ${DATA_ROOT}/raw/refseq_viral/viral.1.genomic.gbff.gz
```

For a larger viral set, list all chunks first:

```bash
curl -L https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/ \
  | grep -o 'viral\.[0-9][0-9]*\.genomic\.gbff\.gz' \
  | sort -u
```

Then download the listed chunks as needed.

## Broader RefSeq Release Downloads

For a larger general pretraining corpus, download release-level RefSeq chunks.
These files cover many species and are better suited to broad representation
learning than only human and mouse.

Mammalian RNA transcript chunks:

```bash
mkdir -p ${DATA_ROOT}/raw/refseq_release/vertebrate_mammalian_rna
python3 scripts/data_pipeline/download_refseq_release.py \
  --base-url https://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_mammalian/ \
  --pattern "*.rna.gbff.gz" \
  --output-dir ${DATA_ROOT}/raw/refseq_release/vertebrate_mammalian_rna
gzip -t ${DATA_ROOT}/raw/refseq_release/vertebrate_mammalian_rna/*.rna.gbff.gz
```

Viral genomic chunks:

```bash
mkdir -p ${DATA_ROOT}/raw/refseq_release/viral_genomic
python3 scripts/data_pipeline/download_refseq_release.py \
  --base-url https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/ \
  --pattern "*.genomic.gbff.gz" \
  --output-dir ${DATA_ROOT}/raw/refseq_release/viral_genomic
gzip -t ${DATA_ROOT}/raw/refseq_release/viral_genomic/*.genomic.gbff.gz
```

Optional non-viral pathogen/generalization chunks:

```bash
mkdir -p ${DATA_ROOT}/raw/refseq_release/fungi_genomic ${DATA_ROOT}/raw/refseq_release/protozoa_genomic
python3 scripts/data_pipeline/download_refseq_release.py \
  --base-url https://ftp.ncbi.nlm.nih.gov/refseq/release/fungi/ \
  --pattern "*.genomic.gbff.gz" \
  --output-dir ${DATA_ROOT}/raw/refseq_release/fungi_genomic
python3 scripts/data_pipeline/download_refseq_release.py \
  --base-url https://ftp.ncbi.nlm.nih.gov/refseq/release/protozoa/ \
  --pattern "*.genomic.gbff.gz" \
  --output-dir ${DATA_ROOT}/raw/refseq_release/protozoa_genomic
gzip -t ${DATA_ROOT}/raw/refseq_release/fungi_genomic/*.genomic.gbff.gz
gzip -t ${DATA_ROOT}/raw/refseq_release/protozoa_genomic/*.genomic.gbff.gz
```

## Processing Commands

Parse RefSeq human transcript GBFF:

```bash
python3 scripts/data_pipeline/genbank_refseq.py \
  --input ${DATA_ROOT}/raw/refseq_human/GRCh38_latest_rna.gbff.gz \
  --output ${DATA_ROOT}/normalized/refseq_human.tsv \
  --source-db refseq_human
```

Parse broader RefSeq release chunks:

```bash
mkdir -p ${DATA_ROOT}/normalized/refseq_release_vertebrate_mammalian
for f in ${DATA_ROOT}/raw/refseq_release/vertebrate_mammalian_rna/*.rna.gbff.gz; do
  base=$(basename "$f" .gz)
  python3 scripts/data_pipeline/genbank_refseq.py \
    --input "$f" \
    --output ${DATA_ROOT}/normalized/refseq_release_vertebrate_mammalian/${base}.tsv \
    --source-db refseq_vertebrate_mammalian
done

mkdir -p ${DATA_ROOT}/normalized/refseq_release_viral
for f in ${DATA_ROOT}/raw/refseq_release/viral_genomic/*.genomic.gbff.gz; do
  base=$(basename "$f" .gz)
  python3 scripts/data_pipeline/genbank_refseq.py \
    --input "$f" \
    --output ${DATA_ROOT}/normalized/refseq_release_viral/${base}.tsv \
    --source-db refseq_viral
done
```

Parse RefSeq viral GBFF:

```bash
python3 scripts/data_pipeline/genbank_refseq.py \
  --input ${DATA_ROOT}/raw/refseq_viral/viral.1.genomic.gbff.gz \
  --output ${DATA_ROOT}/normalized/refseq_viral_1.tsv \
  --source-db refseq_viral
```

Merge and exact-deduplicate. Prefer this directory-aware command because it
skips missing directories and empty TSV files:

```bash
python3 scripts/data_pipeline/merge_normalized.py \
  --inputs ${DATA_ROOT}/normalized/refseq_human.tsv ${DATA_ROOT}/normalized/refseq_release_vertebrate_mammalian ${DATA_ROOT}/normalized/refseq_release_viral \
  --output ${DATA_ROOT}/normalized/all_refseq_unique.tsv \
  --mode protein_orf
```

Create first-pass exact protein-cluster splits:

```bash
python3 scripts/data_pipeline/make_exact_splits.py \
  --input ${DATA_ROOT}/normalized/all_refseq_unique.tsv \
  --output-dir ${DATA_ROOT}/splits/exact_protein \
  --train-fraction 0.8 \
  --val-fraction 0.1 \
  --seed 13
```

## Notes

- Human and broader mammalian transcript files are good for ORF+UTR records.
- Viral genomic records are good for ORF/pathogen records, but UTRs are usually
  absent or not comparable to mammalian mRNA UTRs.
- For antigen-specific fine-tuning, join UniProt/IEDB metadata after this
  normalized base table exists.
- For final publication-quality splits, replace exact protein clusters with
  MMseqs2/CD-HIT similarity clusters.
