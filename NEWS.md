# EchoGO NEWS

## EchoGO 0.1.4

### Interpretation model
- Retire `consensus_score`, `consensus_score_all`, cross-context composite
  significance, and score-era interpretation. Exact-term evidence now reports
  **Target only**, **Target + context**, and **Context-derived hypotheses**, with
  target GOseq retained as the experimental anchor.
- Keep adjusted p-values and fold enrichment local to their originating GOseq or
  g:Profiler source; context recurrence is descriptive provenance, not a
  probability, conservation claim, or independent replication.

### Explicit biological configuration
- Require researcher-selected enrichment contexts for own-data workflows.
- Require a queried `target_context` or an explicit no-target declaration;
  target GOseq remains the experimental anchor in either case.
- Add optional context-selection rationale and explicit default-domain
  exploratory opt-in.
- Resolve the DE foreground and tested-feature background by the same procedure
  and submit the same portable vectors across selected annotation contexts.
  Record context-specific recognition and effective query/domain metadata;
  do not construct species-specific ortholog vectors.
- Keep default-domain exploratory evidence outside primary profiles, recurrence,
  semantic products, networks and evaluation products.

### RRvGO semantic summaries
- Prefer `semantic_reference_orgdb` and require an explicit
  `semantic_reference_role` when RRvGO is enabled; retain `orgdb` as a
  warning-bearing compatibility alias.
- Freeze semantic similarity to `Rel` and separate target-supported semantic
  themes from alternative-context hypothesis semantic themes.

### Offline quickstart and outputs
- Rebuild the deterministic demonstration around explicit v0.1.4 configuration
  and cached offline responses; remove the advertised v0.1.3 frozen score-era
  result tree.
- Lead successful runs to a biological HTML report, exact-term evidence,
  representability audit, and separated semantic products where available.
- Make the full offline demonstration easy to discover with
  `echogo_quickstart(run_demo = TRUE, full = TRUE, live_gprofiler = FALSE)`.
- Improve report navigation and Key Findings, expose semantic plots and network
  assets in reports, and preserve complete demo trees after package installation.
- Add a readable startup panel with optional colour and a monochrome fallback.

### Documentation and migration
- Add a biology-first README, reference-rich zebrafish and annotation-limited
  Gammarus examples, and a v0.1.3-to-v0.1.4 migration guide.
- Rename public report wording to matched experimental-background and
  default-domain exploratory terminology.
- Separate historical validation enrichment from current v0.1.4 interpretation
  in the accompanying compendium and document the deterministic,
  recurrence-neutral Figure 2 selection and its reproducibility limits.

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
- Validate that quickstart produces a non-empty score-era consensus table (legacy v0.1.3 output) and, when requested, an HTML report.
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
