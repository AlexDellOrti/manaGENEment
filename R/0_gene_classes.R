
# Gene class 
#' Gene Class Initialization (Virtual class object)
#'
#' A virtual S4 class representing a generic Gene. It defines
#' slots common to various gene types.
#' 
#' @author Alexander Dell'Orti \cr Politecnico di Milano\cr Maintainer: Alexander
#' Dell'Orti \cr E-Mail: <alexander.dellorti@@mail.polimi.it>
#'
#' @slot geneID Character, Ensembl or NCBI Gene identifier
#' @slot h_symbol Character, Gene HUGO symbol
#' @slot geneName Character, Full gene name
#' @slot description Character, Gene description
#' @slot geneStructure GRanges, Genomic Ranges object representing the gene's structure
#' @slot RnaID Character, RNA identifier
#' @slot RnaSEQ RNAString, RNA sequence in RNAString format
#'
#' @exportClass Gene
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom Biostrings RNAString
#' @importFrom IRanges IRanges
#' @importFrom methods new
#'
#' @seealso \code{\linkS4class{ProteinCodingGene}}, \code{\linkS4class{lncRNAGene}}, \code{\linkS4class{miRNAGene}}
#'
#' @examples
#' # Gene class is virtual and objects cannot be created directly from the user.
#' # Attempting to create a Gene object will result in an error.
setClass("Gene" , 
         representation(
           geneID="character" , 
           h_symbol = "character" , 
           geneName="character" , 
           description="character" ,
           geneStructure = "GRanges",
           RnaID = "character" , 
           RnaSEQ = "RNAString", 
           "VIRTUAL"),
         prototype =list(
           geneID= NA_character_ , 
           h_symbol = NA_character_ , 
           geneName=NA_character_ , 
           description=NA_character_ ,
           geneStructure = GRanges(),
           RnaID = NA_character_ , 
           RnaSEQ = Biostrings::RNAString(""))
         )
         
         
         
# Protein coding RNA class (constructr)
#' ProteinCodingGene Class (S4 class object)
#'
#' An S4 class representing protein-coding genes, inheriting from the virtual Gene class.
#' It includes additional slots specific to protein products.
#' 
#' @author Alexander Dell'Orti \cr Politecnico di Milano\cr Maintainer: Alexander
#' Dell'Orti \cr E-Mail: <alexander.dellorti@@mail.polimi.it>
#'
#' @slot proteinID Character, Protein ID
#' @slot proteinSEQ AAString, Protein amino acid sequence in AAString format
#'
#' @exportClass ProteinCodingGene
#'
#' @importFrom Biostrings AAString
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom methods new
#'
#' @seealso \code{\linkS4class{Gene}}
#'
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#'
#' pcg_object <- new("ProteinCodingGene",
#'            geneID = "ENST000001",
#'            h_symbol = "ACTB",
#'            geneName = "Beta Actin",
#'            description = "A major component of the cytoskeleton",
#'            geneStructure = GenomicRanges::GRanges(seqnames = "chr7",
#'                                     ranges = IRanges::IRanges(start = 5527175, end = 5535161),
#'                                     strand = "+"),
#'            RnaID = "ENST_RNA_001",
#'            RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA"),
#'            proteinID = "ENSP000001",
#'            proteinSEQ = AAString("MEEEIAALVIDNGSG"))
#'
setClass("ProteinCodingGene" , 
         representation( 
           proteinID = "character" , 
           proteinSEQ = "AAString"),
         prototype = list(
           proteinID = NA_character_ , 
           proteinSEQ = Biostrings::AAString("")),
         contains = "Gene")




# long non coding RNA class 
#' lncRNAGene Class (S4 class object)
#'
#' An S4 class representing long non-coding RNA (lncRNA) genes, inheriting from the virtual Gene class.
#' It does not include additional slots beyond those inherited from Gene since the RNA is it's final product.
#'
#' @author Alexander Dell'Orti \cr Politecnico di Milano\cr Maintainer: Alexander
#' Dell'Orti \cr E-Mail: <alexander.dellorti@@mail.polimi.it>
#'
#' @exportClass lncRNAGene
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom Biostrings RNAString
#' @importFrom IRanges IRanges
#' @importFrom methods new
#'
#' @seealso \code{\linkS4class{Gene}}
#'
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#'
#' lncr_object <- new("lncRNAGene",
#'            geneID = "ENST000002",
#'            h_symbol = "LNC1",
#'            geneName = "Long Non-Coding RNA 1",
#'            description = "A long non-coding RNA gene involved in regulation",
#'            geneStructure = GenomicRanges::GRanges(seqnames = "chr2",
#'                                     ranges = IRanges::IRanges(start = 1500000, end = 1501000),
#'                                     strand = "-"),
#'            RnaID = "ENST_RNA_002",
#'            RnaSEQ = RNAString("AUCGAUCGAUGCA"))
#'
setClass("lncRNAGene",
         contains = "Gene")





# micro RNA class
#' miRNAGene Class (S4 class object)
#'
#' An S4 class representing microRNA (miRNA) genes, inheriting from the virtual Gene class.
#' It includes a slot specific to miRNA seed sequences.
#' 
#' @author Alexander Dell'Orti \cr Politecnico di Milano\cr Maintainer: Alexander
#' Dell'Orti \cr E-Mail: <alexander.dellorti@@mail.polimi.it>
#'
#' @slot mirna_seed_SEQ RNAString, miRNA seed sequence in RNAString format
#'
#' @exportClass miRNAGene
#'
#' @importFrom Biostrings RNAString
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom methods new
#'
#' @seealso \code{\linkS4class{Gene}}
#'
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' library(IRanges)
#'
#' mir_object <- new("miRNAGene",
#'            geneID = "ENST000003",
#'            h_symbol = "MIR1",
#'            geneName = "MicroRNA 1",
#'            description = "A microRNA gene involved in post-transcriptional regulation",
#'            geneStructure = GenomicRanges::GRanges(seqnames = "chr3",
#'                                     ranges = IRanges::IRanges(start = 2500000, end = 2500500),
#'                                     strand = "+"),
#'            RnaID = "ENST_RNA_003",
#'            RnaSEQ = RNAString("AUGCUAGCUA"),
#'            mirna_seed_SEQ = RNAString("AUCGUAA"))
#'
setClass("miRNAGene" , 
         representation( 
           mirna_seed_SEQ = "RNAString"),
         prototype = list(
           mirna_seed_SEQ = Biostrings::RNAString("")),
         contains = "Gene")



                    
            