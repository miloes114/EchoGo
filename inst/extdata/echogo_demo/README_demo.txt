EchoGO v0.1.3 deterministic demonstration dataset
=================================================

This synthetic scientific fixture contains no research measurements or
machine-specific paths.

Contract checks represented by the fixture:
- 120 unique tested genes in counts and the full DE table
- 24 significant genes using the explicit logical 'significant' column
- 118 annotation-matched genes, two intentionally unmapped tested genes
- one intentional many-to-one canonical-name collision
- GOseq rows with explicit 24/120 denominators; the first term has fold 5
- a strict resolved foreground subset of the resolved background
- deterministic cached g:Profiler responses for offline quickstart testing

Offline quickstart (default):
  echogo_quickstart(run_demo = TRUE)

Optional live integration run:
  echogo_quickstart(run_demo = TRUE, live_gprofiler = TRUE)

Regeneration:
  source('data-raw/build_demo_v013.R')
