% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc_alignment.R
\name{qc_alignment}
\alias{qc_alignment}
\title{Generate some alignment statistics on reads}
\usage{
qc_alignment(samtools = "samtools", bam.files = NULL,
  filtered.bam.files = NULL, threads = 1,
  remove.mitochondrial = "MT")
}
\arguments{
\item{samtools}{character string specifying path to samtools. [DEFAULT = "samtools"]}

\item{bam.files}{character vector specifying path to pre-filtered bam files. [DEFAULT = NULL]}

\item{filtered.bam.files}{character vector specifying path to filtered bam files. [DEFAULT = NULL]}

\item{threads}{positive integer specifying the number of cores to use}

\item{remove.mitochondrial}{character string.  If set, this will count reads mapping to the
mitochondrial genome.  The string should match the reference name for the mitochondrial genome
in the alignment file.  Examples include "ChrM", "M" and "MT".}
}
\description{
Generate some alignment statistics on reads
}
