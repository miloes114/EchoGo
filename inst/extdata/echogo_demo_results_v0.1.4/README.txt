EchoGO v0.1.4 deterministic demonstration
===========================================

This package ships frozen demo inputs and cached matched-background g:Profiler
responses, not a pre-rendered score-era result tree.

Generate the current biological report offline:

  echogo_quickstart(run_demo = TRUE)

For the optional separated RRvGO semantic products and default-domain
exploratory tier, after installing the explicitly declared semantic dependencies:

  echogo_quickstart(run_demo = TRUE, full = TRUE)

The standard demonstration declares zebrafish as the target annotation context,
uses mouse and human as researcher-selected alternative annotation contexts, and
uses org.Dr.eg.db as a target semantic reference only for the optional RRvGO
step. Cached responses are deterministic and no g:Profiler network request is
made unless live_gprofiler = TRUE is explicitly supplied.
