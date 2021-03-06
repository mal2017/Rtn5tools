% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statsBam.R
\name{reads_aligned}
\alias{reads_aligned}
\title{Count aligned reads in a bam file.}
\usage{
reads_aligned(bamfiles, idx, primary_align_only = F, return_tibble = T,
  millions = F, ...)
}
\value{
A list of integer values corresponding to counts.
}
\description{
Convenience wrapper for a few
}
\details{
When primary_align_only is set to FALSE and no extra arguments
are given this will use an index if one exists.
}
