% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/0_gene_classes.R
\docType{class}
\name{Gene-class}
\alias{Gene-class}
\title{Gene Class Initialization (Virtual class object)}
\description{
A virtual S4 class representing a generic Gene. It defines
slots common to various gene types.
}
\section{Slots}{

\describe{
\item{\code{geneID}}{Character, Ensembl or NCBI Gene identifier}

\item{\code{h_symbol}}{Character, Gene HUGO symbol}

\item{\code{geneName}}{Character, Full gene name}

\item{\code{description}}{Character, Gene description}

\item{\code{geneStructure}}{GRanges, Genomic Ranges object representing the gene's structure}

\item{\code{RnaID}}{Character, RNA identifier}

\item{\code{RnaSEQ}}{RNAString, RNA sequence in RNAString format}
}}

\examples{
# Gene class is virtual and objects cannot be created directly from the user.
# Attempting to create a Gene object will result in an error.
}
\seealso{
\code{\linkS4class{ProteinCodingGene}}, \code{\linkS4class{lncRNAGene}}, \code{\linkS4class{miRNAGene}}
}
\author{
Alexander Dell'Orti \cr Politecnico di Milano\cr Maintainer: Alexander
Dell'Orti \cr E-Mail: \href{mailto:alexander.dellorti@mail.polimi.it}{alexander.dellorti@mail.polimi.it}
}
