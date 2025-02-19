% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_constructors.R
\name{ProteinCodingGene}
\alias{ProteinCodingGene}
\title{Create a ProteinCodingGene Object}
\usage{
ProteinCodingGene(
  geneID,
  h_symbol,
  geneName,
  description,
  geneStructure,
  RnaID,
  RnaSEQ,
  proteinID,
  proteinSEQ
)
}
\arguments{
\item{geneID}{Character, Ensembl or NCBI Gene identifier.}

\item{h_symbol}{Character, Gene HUGO symbol.}

\item{geneName}{Character, Full gene name.}

\item{description}{Character, Gene description.}

\item{geneStructure}{GRanges, Genomic Ranges object representing the gene's structure.}

\item{RnaID}{Character, RNA identifier.}

\item{RnaSEQ}{RNAString, RNA sequence.}

\item{proteinID}{Character, Protein identifier.}

\item{proteinSEQ}{AAString, Protein amino acid sequence.}
}
\value{
ProteinCodingGene S4 object
}
\description{
This constructor function creates a new ProteinCodingGene class object with input validation.
If any of the attributes are missing, the default value defined by the class prototype parameter
will be assigned and a warning message will be printed containing the list of missing parameters.
This behavior will allow the user to insert entries with incomplete information, to be filled in
a second moment. At least one argument has to be given.
}
\examples{
library(GenomicRanges)
library(Biostrings)
library(IRanges)

# Creating a ProteinCodingGene object (with no missing attributes)
pcg_object <- ProteinCodingGene(
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
  RnaSEQ = Biostrings::RNAString("AUGGCUAGCUAGCUAGCUA"),
  proteinID = "ENSP000001",
  proteinSEQ = Biostrings::AAString("MEEEIAALVIDNGSG")
)
}
\seealso{
\code{\linkS4class{Gene}}, \code{\linkS4class{lncRNAGene}}, \code{\linkS4class{miRNAGene}}
}
\author{
Alexander Dell'Orti \cr Politecnico di Milano\cr Maintainer: Alexander
Dell'Orti \cr E-Mail: \href{mailto:alexander.dellorti@mail.polimi.it}{alexander.dellorti@mail.polimi.it}
}
