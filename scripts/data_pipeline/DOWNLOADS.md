# RefSeq input preparation (Nibi)

This document lists **only the active download → normalize → merge path** for the mRNAWM homology split. The original download and parsing Python implementations have been copied **without changing their behavior** from `mRNA_pr`. Run on Nibi in the project's dedicated data directory; do not use the obsolete `exact_protein` splitter.

Install Biopython in your Nibi Python environment if needed, then set:

```bash
cd /path/to/mRNAWM
export DATA_ROOT="${SCRATCH:-$HOME/scratch}/mRNAWM"
mkdir -p "$DATA_ROOT/raw/refseq_human" \
         "$DATA_ROOT/raw/refseq_release/vertebrate_mammalian_rna" \
         "$DATA_ROOT/raw/refseq_release/viral_genomic" \
         "$DATA_ROOT/normalized/refseq_release_vertebrate_mammalian" \
         "$DATA_ROOT/normalized/refseq_release_viral"
python3 -c "import Bio; print(Bio.__version__)"
```

## Download (same URLs, formats and downloader as old project)

```bash
curl -fL \
  https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.gbff.gz \
  -o "$DATA_ROOT/raw/refseq_human/GRCh38_latest_rna.gbff.gz"

python3 scripts/data_pipeline/download_refseq_release.py \
  --base-url https://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_mammalian/ \
  --pattern "*.rna.gbff.gz" \
  --output-dir "$DATA_ROOT/raw/refseq_release/vertebrate_mammalian_rna"

python3 scripts/data_pipeline/download_refseq_release.py \
  --base-url https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/ \
  --pattern "*.genomic.gbff.gz" \
  --output-dir "$DATA_ROOT/raw/refseq_release/viral_genomic"
```

Record downloaded file names, sizes, checksums, and the NCBI release/download date. Download nodes must have outbound internet access; if Nibi compute nodes do not, download using a permitted network-enabled environment before submitting parse jobs.

## Parse

```bash
python3 scripts/data_pipeline/genbank_refseq.py \
  --input "$DATA_ROOT/raw/refseq_human/GRCh38_latest_rna.gbff.gz" \
  --output "$DATA_ROOT/normalized/refseq_human.tsv" \
  --source-db refseq_human

for f in "$DATA_ROOT"/raw/refseq_release/vertebrate_mammalian_rna/*.rna.gbff.gz; do
  test -f "$f" || { echo "Missing mammalian GBFF files" >&2; exit 1; }
  base=$(basename "$f" .gz)
  python3 scripts/data_pipeline/genbank_refseq.py \
    --input "$f" \
    --output "$DATA_ROOT/normalized/refseq_release_vertebrate_mammalian/$base.tsv" \
    --source-db refseq_vertebrate_mammalian
done

for f in "$DATA_ROOT"/raw/refseq_release/viral_genomic/*.genomic.gbff.gz; do
  test -f "$f" || { echo "Missing viral GBFF files" >&2; exit 1; }
  base=$(basename "$f" .gz)
  python3 scripts/data_pipeline/genbank_refseq.py \
    --input "$f" \
    --output "$DATA_ROOT/normalized/refseq_release_viral/$base.tsv" \
    --source-db refseq_viral
done
```

## Merge and QC

```bash
python3 scripts/data_pipeline/merge_normalized.py \
  --inputs "$DATA_ROOT/normalized/refseq_human.tsv" \
           "$DATA_ROOT/normalized/refseq_release_vertebrate_mammalian" \
           "$DATA_ROOT/normalized/refseq_release_viral" \
  --output "$DATA_ROOT/normalized/all_refseq_mammalian_viral_unique.tsv" \
  --mode protein_orf

python3 scripts/data_pipeline/qc_summary.py \
  "$DATA_ROOT/normalized/all_refseq_mammalian_viral_unique.tsv"

sha256sum "$DATA_ROOT/normalized/all_refseq_mammalian_viral_unique.tsv" \
  > "$DATA_ROOT/normalized/all_refseq_mammalian_viral_unique.tsv.sha256"
```

Then follow the **root README's CPU-only Slurm commands** for homology-aware splitting and leakage audit. Keep all data records, including rows lacking UTRs; the split helper reports full-mRNA-eligible counts separately. Viral genomic sequences may lack UTRs, so those rows are not automatically treated as full-mRNA model-ready examples.

**Execution status:** commands documented; full download/parse/merge not launched here. Do not assume existing Nibi data until verified.
