EchoGO demo dataset (distilled from a completed full run)
---------------------------------------------------------
Files:
 - GOseq_enrichment_demo.csv : GOseq subset with gene_ids RESTRICTED to mapped transcripts
 - Trinotate_demo.tsv        : TRINITY -> SYMBOL mapping (synthetic, Metazoa-tagged)
 - DE_results_demo.tsv       : DE subset for mapped transcripts (if DE source found)
 - counts_demo.tsv           : Counts subset for mapped transcripts (if counts source found)

Run the demo:
 demo_dir <- system.file('extdata','echogo_demo', package='EchoGO')
 echogo_quickstart(run_demo = TRUE)
