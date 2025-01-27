# Protein Coding Gene constructor
#' Create a ProteinCodingGene Object
#'
#' This constructor function creates a new ProteinCodingGene class object with input validation.
#' If any of the attributes are missing, the default value defined by the class prototype parameter
#' will be assigned and a warning message will be printed containing the list of missing parameters.
#' This behavior will allow the user to insert entries with incomplete information, to be filled in 
#' a second moment. At least one argument has to be given.
#' 
#' 
#' @author Alexander Dell'Orti \cr Politecnico di Milano\cr Maintainer: Alexander
#' Dell'Orti \cr E-Mail: <alexander.dellorti@@mail.polimi.it>
#'
#' @param geneID Character, Ensembl or NCBI Gene identifier.
#' @param h_symbol Character, Gene HUGO symbol.
#' @param geneName Character, Full gene name.
#' @param description Character, Gene description.
#' @param geneStructure GRanges, Genomic Ranges object representing the gene's structure.
#' @param RnaID Character, RNA identifier.
#' @param RnaSEQ RNAString, RNA sequence.
#' @param proteinID Character, Protein identifier.
#' @param proteinSEQ AAString, Protein amino acid sequence.
#'
#' @return ProteinCodingGene S4 object
#'
#' @export
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom Biostrings RNAString AAString
#' @importFrom IRanges IRanges
#' @importFrom methods new
#' 
#'
#' @seealso \code{\linkS4class{Gene}}, \code{\linkS4class{lncRNAGene}}, \code{\linkS4class{miRNAGene}}
#'
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#'
#' # Creating a ProteinCodingGene object (with no missing attributes)
#' pcg_object <- ProteinCodingGene(
#'   geneID = "ENST000001",
#'   h_symbol = "ACTB",
#'   geneName = "Beta Actin",
#'   description = "A major component of the cytoskeleton",
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
ProteinCodingGene <- function(geneID, h_symbol, geneName, description, geneStructure, RnaID, RnaSEQ, proteinID, proteinSEQ){

  proto <- methods::getClass("ProteinCodingGene")@prototype
  missing_param= c()
  
  if (missing(geneID)) {
    geneID <- proto@geneID
    missing_param <- c(missing_param , "geneID")
  }
  
  if (missing(h_symbol)) {
    h_symbol <- proto@h_symbol
    missing_param <- c(missing_param , "h_symbol")
  }
  if (missing(geneName)) {
    geneName <- proto@geneName
    missing_param <- c(missing_param , "geneName")
  }
  if (missing(description)) {
    description <- proto@description
    missing_param <- c(missing_param , "description")
  }
  if (missing(geneStructure)) {
    geneStructure <- proto@geneStructure
    missing_param <- c(missing_param , "geneStructure")
  }
  if (missing(RnaID)) {
    RnaID <- proto@RnaID
    missing_param <- c(missing_param , "RnaID")
  }
  if (missing(RnaSEQ)) {
    RnaSEQ <- proto@RnaSEQ
    missing_param <- c(missing_param , "RnaSEQ")
  }
  if (missing(proteinID)) {
    proteinID <- proto@proteinID
    missing_param <- c(missing_param , "proteinID")
  }
  if (missing(proteinSEQ)) {
    proteinSEQ <- proto@proteinSEQ
    missing_param <- c(missing_param , "proteinSEQ")
  }
  
  if (length(missing_param) == length(methods::getClass("ProteinCodingGene")@slots)){
    stop("Insert at least one argument")
  }
    
  if (length(missing_param) >0) {
    warning(paste("Warning! Missing parameter(s): " , paste(missing_param , collapse = ",")))
  }
    
  new("ProteinCodingGene" , 
      geneID = geneID, 
      h_symbol = h_symbol, 
      geneName = geneName, 
      description = description, 
      geneStructure = geneStructure,
      RnaID = RnaID,
      RnaSEQ = RnaSEQ,
      proteinID = proteinID, 
      proteinSEQ = proteinSEQ)
  

}





# long non coding RNA class constructor

#' Create a lncRNAGene Object
#'
#' This constructor function creates a new lncRNAGene class object with input validation.
#' If any of the attributes are missing, the default value defined by the class prototype parameter
#' will be assigned and a warning message will be printed.
#'
#' @author Alexander Dell'Orti \cr Politecnico di Milano\cr Maintainer: Alexander
#' Dell'Orti \cr E-Mail: <alexander.dellorti@@mail.polimi.it>
#'
#' @param geneID Character, Ensembl or NCBI Gene identifier
#' @param h_symbol Character, Gene HUGO symbol
#' @param geneName Character, Full gene name
#' @param description Character, Gene description
#' @param geneStructure GRanges, Genomic Ranges object representing the gene's structure
#' @param RnaID Character, RNA identifier
#' @param RnaSEQ RNAString, RNA sequence
#'
#' @return lncRNAGene S4 object
#'
#' @export
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom Biostrings RNAString
#' @importFrom IRanges IRanges
#' @importFrom methods new
#'
#' @seealso \code{\linkS4class{Gene}}, \code{\linkS4class{ProteinCodingGene}}, \code{\linkS4class{miRNAGene}}
#'
#' @examples
#' # Load necessary libraries
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#'
#' # Creating a lncRNAGene object (with all attributes)
#' lnc_object <- lncRNAGene(
#'   geneID = "ENST000001",
#'   h_symbol = "LNC1",
#'   geneName = "Long Non-Coding RNA 1",
#'   description = "A long non-coding RNA gene involved in regulation",
#'   geneStructure = GenomicRanges::GRanges(
#'     seqnames = "chr2",
#'     ranges = IRanges::IRanges(start = 1500000, end = 1501000),
#'     strand = "-"
#'   ),
#'   RnaID = "ENST_RNA_002",
#'   RnaSEQ = Biostrings::RNAString("AUCGAUCGAUGCA")
#' )
#'
lncRNAGene <- function(geneID, h_symbol, geneName, description, geneStructure, RnaID, RnaSEQ){
  
  proto <- methods::getClass("lncRNAGene")@prototype
  missing_param= c()
  
  if (missing(geneID)) {
    geneID <- proto@geneID
    missing_param <- c(missing_param , "geneID")
  }
  
  if (missing(h_symbol)) {
    h_symbol <- proto@h_symbol
    missing_param <- c(missing_param , "h_symbol")
  }
  if (missing(geneName)) {
    geneName <- proto@geneName
    missing_param <- c(missing_param , "geneName")
  }
  if (missing(description)) {
    description <- proto@description
    missing_param <- c(missing_param , "description")
  }
  if (missing(geneStructure)) {
    geneStructure <- proto@geneStructure
    missing_param <- c(missing_param , "geneStructure")
  }
  if (missing(RnaID)) {
    RnaID <- proto@RnaID
    missing_param <- c(missing_param , "RnaID")
  }
  if (missing(RnaSEQ)) {
    RnaSEQ <- proto@RnaSEQ
    missing_param <- c(missing_param , "RnaSEQ")
  }

  
  if (length(missing_param) == length(methods::getClass("lncRNAGene")@slots)){
    stop("Insert at least one argument")
  }
  
  if (length(missing_param) >0) {
    warning(paste("Warning! Missing parameter(s):" , paste(missing_param , collapse = ",")))
  }
  
  new("lncRNAGene" , 
      geneID = geneID, 
      h_symbol = h_symbol, 
      geneName = geneName, 
      description = description, 
      geneStructure = geneStructure,
      RnaID = RnaID,
      RnaSEQ = RnaSEQ)
}





# microRNA class constructor 
#' Create a miRNAGene Object
#'
#' This constructor function creates a new miRNAGene class object with input validation.
#' If any of the attributes are missing, the default value defined by the class prototype parameter
#' will be assigned and a warning message will be printed.
#' 
#' @author Alexander Dell'Orti \cr Politecnico di Milano\cr Maintainer: Alexander
#' Dell'Orti \cr E-Mail: <alexander.dellorti@@mail.polimi.it>
#'
#' @param geneID Character, Ensembl or NCBI Gene identifier
#' @param h_symbol Character, Gene HUGO symbol
#' @param geneName Character, Full gene name
#' @param description Character, Gene description
#' @param geneStructure GRanges, Genomic Ranges object representing the gene's structure
#' @param RnaID Character, RNA identifier
#' @param RnaSEQ RNAString, RNA sequence
#' @param mirna_seed_SEQ RNAString, miRNA seed sequence
#'
#' @return miRNAGene S4 object.
#'
#' @export
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom Biostrings RNAString
#' @importFrom IRanges IRanges
#' @importFrom methods new
#'
#' @seealso \code{\linkS4class{Gene}}, \code{\linkS4class{ProteinCodingGene}}, \code{\linkS4class{lncRNAGene}}
#'
#' @examples
#' # Load necessary libraries
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#'
#' # Creating a miRNAGene object (with all attributes)
#' mir_object <- miRNAGene(
#'   geneID = "ENST000003",
#'   h_symbol = "MIR1",
#'   geneName = "MicroRNA 1",
#'   description = "A microRNA gene involved in post-transcriptional regulation",
#'   geneStructure = GenomicRanges::GRanges(
#'     seqnames = "chr3",
#'     ranges = IRanges::IRanges(start = 2500000, end = 2500500),
#'     strand = "+"
#'   ),
#'   RnaID = "ENST_RNA_003",
#'   RnaSEQ = Biostrings::RNAString("AUGCUAGCUA"),
#'   mirna_seed_SEQ = Biostrings::RNAString("AUCGUAA")
#' )
miRNAGene <- function(geneID, h_symbol, geneName, description, geneStructure, RnaID, RnaSEQ, mirna_seed_SEQ){
  
  proto <- methods::getClass("miRNAGene")@prototype
  missing_param= c()
  
  if (missing(geneID)) {
    geneID <- proto@geneID
    missing_param <- c(missing_param , "geneID")
  }
  
  if (missing(h_symbol)) {
    h_symbol <- proto@h_symbol
    missing_param <- c(missing_param , "h_symbol")
  }
  if (missing(geneName)) {
    geneName <- proto@geneName
    missing_param <- c(missing_param , "geneName")
  }
  if (missing(description)) {
    description <- proto@description
    missing_param <- c(missing_param , "description")
  }
  if (missing(geneStructure)) {
    geneStructure <- proto@geneStructure
    missing_param <- c(missing_param , "geneStructure")
  }
  if (missing(RnaID)) {
    RnaID <- proto@RnaID
    missing_param <- c(missing_param , "RnaID")
  }
  if (missing(RnaSEQ)) {
    RnaSEQ <- proto@RnaSEQ
    missing_param <- c(missing_param , "RnaSEQ")
  }
  if (missing(mirna_seed_SEQ)) {
    mirna_seed_SEQ <- proto@mirna_seed_SEQ
    missing_param <- c(missing_param , "mirna_seed_SEQ")
  }
  
  if (length(missing_param) == length(methods::getClass("miRNAGene")@slots)){
    stop("Insert at least one argument")
  }
  
  if (length(missing_param) >0) {
    warning(paste("Warning! Missing parameter(s):" , paste(missing_param , collapse = ",")))
  }
  
  new("miRNAGene" , 
      geneID = geneID, 
      h_symbol = h_symbol, 
      geneName = geneName, 
      description = description, 
      geneStructure = geneStructure,
      RnaID = RnaID,
      RnaSEQ = RnaSEQ,
      mirna_seed_SEQ = mirna_seed_SEQ)
}