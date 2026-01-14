# EchoGO NEWS

## EchoGO 0.1.2 (2026-01-14)

### Added
- Native input resolver support for **reference-based RNA-seq experiment outputs** (DESeq2 + GOseq precomputed).
  - Recognizes the standard reference-based layout with:
    - `allcounts_table.txt`
    - `dge_<CONTRAST>.csv`
    - `<CONTRAST>.GOseq.enriched.tsv`
    - `Trinotate_for_EchoGO.tsv` and/or `<reference_label>_eggNOG_for_EchoGO.tsv`

### Changed
- Reference-based mode now supports cases where GOseq `gene_ids` are already **gene symbols/names** (not transcript IDs).
  - If GOseq `gene_ids` do not match annotation `transcript_id`, EchoGO assumes the GOseq IDs are already usable names (intended for reference-based pipelines).

### Documentation
- `echogo_scaffold()` now documents **both** de novo and reference-based layouts in `input/_README.txt`.
- Added/updated vignette: *Reference-based RNA-seq inputs for EchoGO*.

## EchoGO 0.1.1
- Previous release.
