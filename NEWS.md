# EchoGO NEWS

## EchoGO 0.1.3 (2026-07-16)

### Fixed
- Restore portable canonical-name resolution for g:Profiler: genuine SwissProt
  gene symbol, then eggNOG preferred name, then portable native symbol. Raw
  transcript/contig, seed-ortholog, Ensembl protein, and accession fallbacks are
  retained in mapping provenance but excluded from submitted vectors.
- Keep reference-preparation fields semantically separate instead of copying
  eggNOG preferred or seed values into the SwissProt column; taxonomy labels are
  provenance rather than numeric resolver cutoffs.
- Require `GO.db` and `rrvgo` for the complete pipeline and report missing RRvGO directly.
- Install OrgDb packages into the active first library, including project-local `renv` libraries.
- Make the packaged demo deterministic by shipping only the canonical GOseq TSV, passing all demo files explicitly, and cleaning stale results by default.
- Rank automatic input patterns, prefer GOseq TSV files, and reject ambiguous matches.
- Preserve correctly quoted comma-separated `gene_ids` values when reading CSV files and detect spill-column corruption.
- Validate that quickstart produces a non-empty consensus table and, when requested, an HTML report.
- Replace the malformed reference-based Markdown guide with a fully rendered, consistently styled R Markdown vignette.
- Filter generated GOseq tables by adjusted FDR rather than raw p-value.
- Honor existing DE `significant` flags and reuse cached reference maps on staged reruns.
- Keep interactive network widgets portable by staging their companion
  dependency directories, expose the existing PDF networks as report fallbacks,
  and report self-contained widget failures before using a dependency-backed
  HTML fallback.
- Ship every local asset referenced by the packaged frozen-demo HTML report so
  the report renders completely after installation.

### Added
- `run_rrvgo` controls in `run_full_echogo()` and `run_echogo_pipeline()` for deliberately lighter runs.
- `echogo_prepare_reference_inputs()` for a resumable GFF3/GTF + protein FASTA + eggNOG-mapper workflow, with optional species-native OrgDb GO mapping.
- Offline regression coverage for demo inputs, parsing, resolution, dependency errors, installer libraries, frozen outputs, and quickstart wiring.
- Synthetic end-to-end reference preparation coverage from annotation import through GOseq and EchoGO loading.

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
