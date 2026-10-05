EchoGO v0.1.4 deterministic demonstration dataset
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
- deterministic cached matched-background g:Profiler responses for offline quickstart testing

The built-in example explicitly declares:
- zebrafish (`drerio`) as its target annotation context;
- mouse and human as researcher-selected alternative annotation contexts;
- default-domain exploration disabled for the basic quickstart;
- `org.Dr.eg.db` as a target semantic reference, used only when `full = TRUE`;
- RRvGO semantic method `Rel` and semantic-reference role `target_reference`.

The cached responses are frozen demonstration evidence. They are not claims that
mouse or human were independently tested in this experiment.

Offline quickstart (default):
  echogo_quickstart(run_demo = TRUE)

Optional live integration run:
  echogo_quickstart(run_demo = TRUE, live_gprofiler = TRUE)

Regeneration is intentionally not part of the installed package. The quickstart
copies these inputs and cached responses, then regenerates current scoreless
v0.1.4 outputs locally.
