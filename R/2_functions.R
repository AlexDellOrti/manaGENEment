

# =====================================================================
#   ***                     Generic Functions                     ***
# =====================================================================


#' Get Gene ID Generic Function
#'
#' Retrieves the Gene ID from a Gene object
#'
#' @param geneobject An object of class Gene
#'
#' @return A character string representing the Gene ID
#' 
#' @export
#' 
#' @importFrom methods new 
#' @importFrom IRanges IRanges
#'
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#' 
#' pcg <- ProteinCodingGene(
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
#'   RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA"),
#'   proteinID = "ENSP000001",
#'   proteinSEQ = AAString("MEEEIAALVIDNGSG")
#' )
#' 
#' getGeneID(pcg)
setGeneric("getGeneID", function(geneobject) standardGeneric("getGeneID"))



#' Set Gene ID Generic Function
#'
#' Sets a new Gene ID for a Gene object
#'
#' @param geneobject An object of class Gene
#' @param newid A character string representing the new Gene ID
#'
#' @return The updated Gene object with the new Gene ID
#'
#' @export
#' 
#' @importFrom methods new
#' @importFrom IRanges IRanges
#'
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#' 
#' pcg <- ProteinCodingGene(
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
#'   RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA"),
#'   proteinID = "ENSP000001",
#'   proteinSEQ = AAString("MEEEIAALVIDNGSG")
#' )
#' 
#' pcg <- setGeneID(pcg, "ENST000002")
#' gene_id <- getGeneID(pcg)
#' print(gene_id)  
setGeneric("setGeneID" , function(geneobject, newid) standardGeneric("setGeneID"))



#' Get RNA ID Generic Function
#'
#' Retrieves the RNA ID from a Gene object
#'
#' @param geneobject An object of class Gene
#'
#' @return A character string representing the RNA ID
#'
#' @export
#' 
#' @importFrom methods new
#' @importFrom IRanges IRanges
#'
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#' 
#' pcg <- ProteinCodingGene(
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
#'   RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA"),
#'   proteinID = "ENSP000001",
#'   proteinSEQ = AAString("MEEEIAALVIDNGSG")
#' )
#' 
#' getRnaID(pcg)
setGeneric("getRnaID" , function(geneobject) standardGeneric("getRnaID"))




#' Get RNA Sequence Generic Function
#'
#' Retrieves the RNA sequence from a Gene object
#'
#' @param geneobject An object of class Gene
#'
#' @return An object of class RNAString representing the RNA Sequence
#'
#' @export
#' 
#' @importFrom methods new
#' @importFrom IRanges IRanges
#'
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#' 
#' pcg <- ProteinCodingGene(
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
#'   RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA"),
#'   proteinID = "ENSP000001",
#'   proteinSEQ = AAString("MEEEIAALVIDNGSG")
#' )
#' 
#' getRnaSEQ(pcg)
setGeneric("getRnaSEQ", function(geneobject) standardGeneric("getRnaSEQ"))



#' Length Product Generic Function
#'
#' Computes the length of the gene's main product, that varies based on the gene type:
#' - ProteinCodingGene: Number of amino acids in the protein sequence
#' - lncRNAGene: Length of the lncRNA sequence
#' - miRNAGene: Length of the miRNA sequence
#'
#' @param geneobject An object of class Gene or from one subclass derived from it
#'
#' @return An integer representing the length of the gene's main product
#'
#' @export
#' 
#' @importFrom methods new
#' @importFrom IRanges IRanges
#'
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#' 
#' pcg <- ProteinCodingGene(
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
#'   RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA"),
#'   proteinID = "ENSP000001",
#'   proteinSEQ = AAString("MEEEIAALVIDNGSG")
#' )
#' 
#' lengthProduct(pcg)
#'
#' lnc <- lncRNAGene(
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
#'   RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA")
#' )
#' 
#' lengthProduct(lnc)
#'
#' mir <- miRNAGene(
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
#'   RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA"),
#'   mirna_seed_SEQ= RNAString("UAGCUA")
#' )
#' 
#' lengthProduct(mir)
setGeneric("lengthProduct" , function(geneobject) standardGeneric("lengthProduct"))



#' Get Protein ID Generic Function
#'
#' Retrieves the Protein ID from a Protein Coding Gene object
#'
#' @param geneobject An object of class ProteinCodingGene
#'
#' @return A character string representing the Protein ID
#'
#' @export
#'
#' @importFrom methods new
#' @importFrom IRanges IRanges
#'
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#' 
#' pcg <- ProteinCodingGene(
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
#'   RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA"),
#'   proteinID = "ENSP000001",
#'   proteinSEQ = AAString("MEEEIAALVIDNGSG")
#' )
#' 
#' getProtID(pcg)
setGeneric("getProtID" , function(geneobject) standardGeneric("getProtID"))



#' Get Protein SEQ Generic Function
#'
#' Retrieves the Protein sequence from a Protein Coding Gene object
#'
#' @param geneobject An object of class ProteinCodingGene
#'
#' @return An object of class AAString representing the protein sequence
#'
#' @export
#' 
#' @importFrom methods new
#' @importFrom IRanges IRanges
#'
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#' 
#' pcg <- ProteinCodingGene(
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
#'   RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA"),
#'   proteinID = "ENSP000001",
#'   proteinSEQ = AAString("MEEEIAALVIDNGSG")
#' )
#' 
#' getProtSeq(pcg)
setGeneric("getProtSeq" , function(geneobject) standardGeneric("getProtSeq"))



#' Get miRNAs Seed Sequence Generic Function
#'
#' Retrieves the miRNA seed sequence from a miRNAGene object
#'
#' @param geneobject An object of class miRNAGene
#'
#' @return An object of class RNAString representing the miRNA seed sequence
#'
#' @export
#' 
#' @importFrom methods new
#' @importFrom IRanges IRanges
#'
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#' 
#' mir <- miRNAGene(
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
#'   RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA"),
#'   mirna_seed_SEQ= RNAString("UAGCUA")
#' )
#' 
#' getSeedSEQ(mir)
setGeneric("getSeedSEQ", function(geneobject) standardGeneric("getSeedSEQ"))



# =====================================================================
#   ***                     Class Methods                     ***
# =====================================================================


# • Gene Virtual Class methods

#' Get Gene ID Method for Gene Class
#'
#' Retrieves the Gene ID from a Gene object
#'
#' @param geneobject An object of class Gene
#'
#' @return A character string representing the Gene ID
#'
#' @exportMethod getGeneID
#' 
#' @importFrom methods new
#' @importFrom IRanges IRanges
#'
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#' 
#' pcg <- ProteinCodingGene(
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
#'   RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA"),
#'   proteinID = "ENSP000001",
#'   proteinSEQ = AAString("MEEEIAALVIDNGSG")
#' )
#' 
#' gene_id <- getGeneID(pcg)
#' print(gene_id)  
setMethod("getGeneID", signature = "Gene" , function(geneobject) {
  return(geneobject@geneID)
})  
  

#' Set Gene ID Method for Gene Class
#'
#' Sets a new Gene ID for a Gene object
#'
#' @param geneobject An object of class Gene
#' @param newid A character string representing the new Gene ID
#'
#' @return The updated Gene object with the new Gene ID
#'
#' @exportMethod setGeneID
#' 
#' @importFrom methods new
#' @importFrom IRanges IRanges
#'
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#' 
#' pcg <- ProteinCodingGene(
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
#'   RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA"),
#'   proteinID = "ENSP000001",
#'   proteinSEQ = AAString("MEEEIAALVIDNGSG")
#' )
#' 
#' pcg <- setGeneID(pcg, "ENST000002")
#' gene_id <- getGeneID(pcg)
#' print(gene_id)  
setMethod("setGeneID", signature = "Gene" , function(geneobject , newid) {
  if (!is.character(newid)) stop("The new ID you want to set has to be of class character.")
  geneobject@geneID <- newid
  geneobject
})  


#' Get RNA ID Method for Gene Class
#'
#' Retrieves the RNA ID from a Gene object
#'
#' @param geneobject An object of class Gene
#'
#' @return A character string representing the RNA ID
#'
#' @exportMethod getRnaID
#' 
#' @importFrom methods new
#' @importFrom IRanges IRanges
#'
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#' 
#' pcg <- ProteinCodingGene(
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
#'   RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA"),
#'   proteinID = "ENSP000001",
#'   proteinSEQ = AAString("MEEEIAALVIDNGSG")
#' )
#' 
#' rna_id <- getRnaID(pcg)
#' print(rna_id)
setMethod("getRnaID", signature = "Gene" , function(geneobject) {
  return(geneobject@RnaID)
}) 
  

#' Get RNA Sequence Method for Gene Class
#'
#' Retrieves the RNA sequence from a Gene object
#'
#' @param geneobject An object of class Gene
#'
#' @return An object of class RNAString representing the RNA sequence
#'
#' @exportMethod getRnaSEQ
#'
#' @importFrom methods new 
#' @importFrom IRanges IRanges
#'
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#' 
#' pcg <- ProteinCodingGene(
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
#'   RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA"),
#'   proteinID = "ENSP000001",
#'   proteinSEQ = AAString("MEEEIAALVIDNGSG")
#' )
#' 
#' getRnaSEQ(pcg)
setMethod("getRnaSEQ", signature = "Gene" , function(geneobject) {
  return(geneobject@RnaSEQ)
}) 






# • Protein Coding Genes Class methods 

#' Get Protein ID Method for ProteinCodingGene Class
#'
#' Retrieves the Protein ID from a ProteinCodingGene object
#'
#' @param geneobject An object of class ProteinCodingGene
#'
#' @return A character string representing the Protein ID
#'
#' @exportMethod getProtID
#' 
#' @importFrom methods new
#' @importFrom IRanges IRanges
#'
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#' 
#' pcg <- ProteinCodingGene(
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
#'   RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA"),
#'   proteinID = "ENSP000001",
#'   proteinSEQ = AAString("MEEEIAALVIDNGSG")
#' )
#' 
#' prot_id <- getProtID(pcg)
#' print(prot_id)  
setMethod("getProtID", signature = "ProteinCodingGene" , function(geneobject) {
  return(geneobject@proteinID)
}) 


#' Get Protein SEQ Method for ProteinCodingGene Class
#'
#' Retrieves the Protein Sequence from a ProteinCodingGene object
#'
#' @param geneobject An object of class ProteinCodingGene
#'
#' @return An object of class AAString representing the Protein sequence
#'
#' @exportMethod getProtSeq
#' 
#' @importFrom methods new
#' @importFrom IRanges IRanges
#'
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#' 
#' pcg <- ProteinCodingGene(
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
#'   RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA"),
#'   proteinID = "ENSP000001",
#'   proteinSEQ = AAString("MEEEIAALVIDNGSG")
#' )
#' 
#' getProtSeq(pcg)
setMethod("getProtSeq", signature = "ProteinCodingGene" , function(geneobject) {
  return(geneobject@proteinSEQ)
}) 



#' Length Product Method for ProteinCodingGene Class
#'
#' Computes the length of the protein sequence in a ProteinCodingGene object
#'
#' @param geneobject An object of class ProteinCodingGene
#'
#' @return An integer representing the number of amino acids in the protein sequence
#'
#' @exportMethod lengthProduct
#' 
#' @importFrom methods new
#' @importFrom IRanges IRanges
#'
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#' 
#' pcg <- ProteinCodingGene(
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
#'   RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA"),
#'   proteinID = "ENSP000001",
#'   proteinSEQ = AAString("MEEEIAALVIDNGSG")
#' )
#' 
#' len_pcg <- lengthProduct(pcg)
#' print(len_pcg) 
setMethod("lengthProduct", signature = "ProteinCodingGene" , function(geneobject) {
  return(length(geneobject@proteinSEQ))
}) 



# • Long non-Coding RNA Genes Class methods 

#' Length Product Method for lncRNAGene Class
#'
#' Computes the length of the RNA sequence in a lncRNAGene object
#'
#' @param geneobject An object of class lncRNAGene
#'
#' @return An integer representing the length of the RNA sequence
#'
#' @exportMethod lengthProduct
#' 
#' @importFrom methods new
#' @importFrom IRanges IRanges
#'
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#' 
#' lnc <- lncRNAGene(
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
#'   RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA")
#' )
#' 
#' len_lnc <- lengthProduct(lnc)
#' print(len_lnc)  
setMethod("lengthProduct", signature = "lncRNAGene" , function(geneobject) {
  return(length(geneobject@RnaSEQ))
}) 



# • micro RNA Genes Class methods 

#' Length Product Method for miRNAGene Class
#'
#' Computes the length of the RNA sequence in a miRNAGene object
#'
#' @param geneobject An object of class miRNAGene
#'
#' @return An integer representing the length of the RNA sequence
#'
#' @exportMethod lengthProduct
#' 
#' @importFrom methods new
#' @importFrom IRanges IRanges
#'
#'
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#' 
#' mir <- miRNAGene(
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
#'   RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA"),
#'   mirna_seed_SEQ= RNAString("UAGCUA")
#' )
#' 
#' len_mir <- lengthProduct(mir)
#' print(len_mir)  
setMethod("lengthProduct", signature = "miRNAGene" , function(geneobject) {
  return(length(geneobject@RnaSEQ))
})  

#' Get Seed Sequence Method for miRNAGene Class
#'
#' Retrieves the miRNA seed sequence from a miRNAGene object.
#'
#' @param geneobject An object of class miRNAGene
#'
#' @return An object of class RNAString representing the miRNA seed sequence
#'
#' @exportMethod getSeedSEQ
#'
#' @importFrom methods new
#' @importFrom IRanges IRanges
#'
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#' 
#' mir <- miRNAGene(
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
#'   RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA"),
#'   mirna_seed_SEQ= RNAString("UAGCUA")
#' )
#' seed_seq <- getSeedSEQ(mir)
#' print(seed_seq)  
setMethod("getSeedSEQ", signature = "miRNAGene" , function(geneobject) {
  return(geneobject@mirna_seed_SEQ)
}) 


