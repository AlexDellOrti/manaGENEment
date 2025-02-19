% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_functions.R
\name{getProtSeq}
\alias{getProtSeq}
\title{Get Protein SEQ Generic Function}
\usage{
getProtSeq(geneobject)
}
\arguments{
\item{geneobject}{An object of class ProteinCodingGene}
}
\value{
An object of class AAString representing the protein sequence
}
\description{
Retrieves the Protein sequence from a Protein Coding Gene object
}
\examples{
library(GenomicRanges)
library(Biostrings)
library(IRanges)

pcg <- ProteinCodingGene(
  geneID = "ENST000001",
  h_symbol = "ACTB",
  geneName = "Beta Actin",
  description = "A major component of the cytoskeleton",
  geneStructure = GenomicRanges::GRanges(
    seqnames = "chr7",
    ranges = IRanges::IRanges(start = 5527175, end = 5535161),
    strand = "+"
  ),
  RnaID = "ENST_RNA_001",
  RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA"),
  proteinID = "ENSP000001",
  proteinSEQ = AAString("MEEEIAALVIDNGSG")
)

getProtSeq(pcg)
}
