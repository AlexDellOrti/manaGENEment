#' manaGENEment - Management of Gene Objects Using S4 Classes
#'
#' The \code{manaGENEment} package provides S4 classe objects and corresponding 
#' methods for managing informations various gene types, each represented by a class, 
#' including protein-coding genes, lncRNA genes, and miRNA genes.
#'
#' \tabular{ll}{
#' Package: \tab manaGENEment\cr
#' Type: \tab Package\cr
#' Version: \tab 0.1.0\cr
#' Date: \tab 2025-01-18\cr
#' License: \tab GPL (>=2)\cr
#' }
#'
#' @name manaGENEment-package
#' @aliases manaGENEment-package manaGENEment
#' @docType package
#' @author Alexander Dell'Orti \cr
#' Politecnico di Milano\cr
#' Maintainer: Alexander Dell'Orti\cr
#' E-Mail: <alexander.dellorti@@mail.polimi.it>
#' @examples
#' # Example usage of manaGENEment package functions
#' library(manaGENEment)
#' 
#' # Create a ProteinCodingGene object
#' pcg <- ProteinCodingGene(
#'   geneID = "ENST000001",
#'   h_symbol = "ACTB",
#'   geneName = "Beta Actin",
#'   description = "One of six different actin isoforms and a major component of the cytoskeleton",
#'   geneStructure = GenomicRanges::GRanges(
#'     seqnames = "chr7",
#'     ranges = IRanges::IRanges(start = 5527175, end = 5535161),
#'     strand = "+"
#'   ),
#'   RnaID = "ENST_RNA_001",
#'   RnaSEQ = Biostrings::RNAString("AUGGCUAGCUAGCUAGCUA"),
#'   proteinID = "ENSP000001",
#'   proteinSEQ = Biostrings::AAString("MEEEIAALVIDNGSG")
#' )
#' 
#' # Get Gene ID
#' getGeneID(pcg)
#' 
#' # Create a lncRNAGene object
#' lnc <- lncRNAGene(
#'   geneID = "ENST000002",
#'   h_symbol = "HOTAIR",
#'   geneName = "HOX Transcript Antisense RNA",
#'   description = "Involved in epigenetic regulation",
#'   geneStructure = GenomicRanges::GRanges(
#'     seqnames = "chr12",
#'     ranges = IRanges::IRanges(start = 1150000, end = 1160000),
#'     strand = "-"
#'   ),
#'   RnaID = "ENST_RNA_002",
#'   RnaSEQ = Biostrings::RNAString("AUGGCUAGCUAGCUAGCUA")
#' )
#' 
#' # Compute Length Product for lncRNAGene
#' len_lnc <- lengthProduct(lnc)
#' print(len_lnc)
NULL