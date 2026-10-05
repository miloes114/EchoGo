% Generated manually for the v0.1.4 exported interface.
\name{echogo_hypothesis_drilldown}
\alias{echogo_hypothesis_drilldown}
\title{Retrieve a traceable EchoGO term evidence trail}
\usage{
echogo_hypothesis_drilldown(evidence_table, source_provenance, term_id,
  annotation_provenance = NULL)
}
\arguments{
\item{evidence_table}{Canonical exact-term evidence table.}
\item{source_provenance}{Long-form source provenance table.}
\item{term_id}{GO term identifier.}
\item{annotation_provenance}{Optional gene-level annotation provenance table.}
}
\value{A list containing the exact term evidence, matching source rows,
normalized context-selection provenance, contributing genes, and any matching
annotation provenance.}
\description{Returns the inspectable evidence trail needed to follow up a
term, especially an \code{ALTERNATIVE_CONTEXT} hypothesis. It does not make a
mechanistic or combined-significance claim.}
