% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_constructors.R
\name{miRNAGene}
\alias{miRNAGene}
\title{Create a miRNAGene Object}
\usage{
miRNAGene(
  geneID,
  h_symbol,
  geneName,
  description,
  geneStructure,
  RnaID,
  RnaSEQ,
  mirna_seed_SEQ
)
}
\arguments{
\item{geneID}{Character, Ensembl or NCBI Gene identifier}

\item{h_symbol}{Character, Gene HUGO symbol}

\item{geneName}{Character, Full gene name}

\item{description}{Character, Gene description}

\item{geneStructure}{GRanges, Genomic Ranges object representing the gene's structure}

\item{RnaID}{Character, RNA identifier}

\item{RnaSEQ}{RNAString, RNA sequence}

\item{mirna_seed_SEQ}{RNAString, miRNA seed sequence}
}
\value{
miRNAGene S4 object.
}
\description{
This constructor function creates a new miRNAGene class object with input validation.
If any of the attributes are missing, the default value defined by the class prototype parameter
will be assigned and a warning message will be printed.
}
\examples{
# Load necessary libraries
library(GenomicRanges)
library(Biostrings)
library(IRanges)

# Creating a miRNAGene object (with all attributes)
mir_object <- miRNAGene(
  geneID = "ENST000003",
  h_symbol = "MIR1",
  geneName = "MicroRNA 1",
  description = "A microRNA gene involved in post-transcriptional regulation",
  geneStructure = GenomicRanges::GRanges(
    seqnames = "chr3",
    ranges = IRanges::IRanges(start = 2500000, end = 2500500),
    strand = "+"
  ),
  RnaID = "ENST_RNA_003",
  RnaSEQ = Biostrings::RNAString("AUGCUAGCUA"),
  mirna_seed_SEQ = Biostrings::RNAString("AUCGUAA")
)
}
\seealso{
\code{\linkS4class{Gene}}, \code{\linkS4class{ProteinCodingGene}}, \code{\linkS4class{lncRNAGene}}
}
\author{
Alexander Dell'Orti \cr Politecnico di Milano\cr Maintainer: Alexander
Dell'Orti \cr E-Mail: \href{mailto:alexander.dellorti@mail.polimi.it}{alexander.dellorti@mail.polimi.it}
}
