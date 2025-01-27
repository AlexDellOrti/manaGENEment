library(testthat)
library(manaGENEment)  
library(GenomicRanges)
library(Biostrings)
library(IRanges)


# testing every class for complete set of input and for missing ones. 
# So one warning message for each missing attribute + default parameters specified in prototype during
# class definitions are expected




## protein coding gene class constructor

# test constructor with correct parameters
test_that("ProteinCodingGene constructor works with all parameters", {
  pcg <- ProteinCodingGene(
    geneID = "ENST000001",
    h_symbol = "ACTB",
    geneName = "Beta Actin",
    description = "A major component of the cytoskeleton",
    geneStructure = GenomicRanges::GRanges(
      seqnames = "chr7",
      ranges = IRanges::IRanges(start = 5527175, end = 5535161), strand = "+"),
    RnaID = "ENST_RNA_001",
    RnaSEQ = Biostrings::RNAString("AUGGCUAGCUAGCUAGCUA"),
    proteinID = "ENSP000001",
    proteinSEQ = Biostrings::AAString("MEEEIAALVIDNGSG")
  )
  
  expect_s4_class(pcg, "ProteinCodingGene")
  expect_equal(pcg@geneID, "ENST000001")
  expect_equal(pcg@h_symbol, "ACTB")
  expect_equal(pcg@geneName, "Beta Actin")
  expect_equal(pcg@description, "A major component of the cytoskeleton")
  expect_true(inherits(pcg@geneStructure, "GRanges"))
  expect_equal(pcg@RnaID, "ENST_RNA_001")
  expect_equal(as.character(pcg@RnaSEQ), "AUGGCUAGCUAGCUAGCUA")
  expect_equal(pcg@proteinID, "ENSP000001")
  expect_equal(as.character(pcg@proteinSEQ), "MEEEIAALVIDNGSG")
})



# test constructor with parameters missing:
# geneID              PRESENT
# h_symbol            PRESENT
# geneName            MISSING
# description         MISSING
# geneStructure       MISSING
# RnaID               MISSING
# RnaSEQ              MISSING
# proteinID           MISSING
# proteinSEQ          MISSING


test_that("ProteinCodingGene constructor assigns defaults and returns warnings for missing parameters", {
  
  # expected warnings
  expect_warning(
    pcg <- ProteinCodingGene(
      geneID = "ENST000001",
      h_symbol = "ACTB"
    ),
    regexp = "Warning! Missing parameter\\(s\\): geneName, description, geneStructure, RnaID, RnaSEQ, proteinID, proteinSEQ"
  )
  
  
  # expected attributes values
  expect_s4_class(pcg, "ProteinCodingGene")
  expect_equal(pcg@geneID, "ENST000001")
  expect_equal(pcg@h_symbol, "ACTB")
  expect_true(is.na(pcg@geneName))
  expect_true(is.na(pcg@description))
  expect_true(inherits(pcg@geneStructure, "GRanges"))
  expect_true(is.na(pcg@RnaID))
  expect_equal(as.character(pcg@RnaSEQ), "")
  expect_true(is.na(pcg@proteinID))
  expect_equal(as.character(pcg@proteinSEQ), "")
})








## long non-coding RNA gene class constructor

# test constructor with correct parameters
test_that("lncRNAGene constructor works with all parameters", {
  lnc <- lncRNAGene(
    geneID = "ENST000001",
    h_symbol = "ACTB",
    geneName = "Beta Actin",
    description = "A major component of the cytoskeleton",
    geneStructure = GenomicRanges::GRanges(
      seqnames = "chr7",
      ranges = IRanges::IRanges(start = 5527175, end = 5535161), strand = "+"),
    RnaID = "ENST_RNA_001",
    RnaSEQ = Biostrings::RNAString("AUGGCUAGCUAGCUAGCUA"),
  )
  
  expect_s4_class(lnc, "lncRNAGene")
  expect_equal(lnc@geneID, "ENST000001")
  expect_equal(lnc@h_symbol, "ACTB")
  expect_equal(lnc@geneName, "Beta Actin")
  expect_equal(lnc@description, "A major component of the cytoskeleton")
  expect_true(inherits(lnc@geneStructure, "GRanges"))
  expect_equal(lnc@RnaID, "ENST_RNA_001")
  expect_equal(as.character(lnc@RnaSEQ), "AUGGCUAGCUAGCUAGCUA")
})




# test constructor with parameters missing:
# geneID              PRESENT
# h_symbol            PRESENT
# geneName            MISSING
# description         MISSING
# geneStructure       MISSING
# RnaID               MISSING
# RnaSEQ              MISSING



test_that("lncRNAGene constructor assigns defaults and returns warnings for missing parameters", {
  
  # expected warnings
  expect_warning(
    lnc <- lncRNAGene(
      geneID = "ENST000001",
      h_symbol = "ACTB"
    ),
    regexp = "Warning! Missing parameter\\(s\\): geneName, description, geneStructure, RnaID, RnaSEQ"
  )
  
  
  # expected attributes values
  expect_s4_class(lnc, "lncRNAGene")
  expect_equal(lnc@geneID, "ENST000001")
  expect_equal(lnc@h_symbol, "ACTB")
  expect_true(is.na(lnc@geneName))
  expect_true(is.na(lnc@description))
  expect_true(inherits(lnc@geneStructure, "GRanges"))
  expect_true(is.na(lnc@RnaID))
  expect_equal(as.character(lnc@RnaSEQ), "")
})












## micro RNA gene class constructor

# test constructor with correct parameters
test_that("miRNAGene constructor works with all parameters", {
  mir <- miRNAGene(
    geneID = "ENST000001",
    h_symbol = "ACTB",
    geneName = "Beta Actin",
    description = "A major component of the cytoskeleton",
    geneStructure = GenomicRanges::GRanges(
      seqnames = "chr7",
      ranges = IRanges::IRanges(start = 5527175, end = 5535161), strand = "+"),
    RnaID = "ENST_RNA_001",
    RnaSEQ = Biostrings::RNAString("AUGGCUAGCUAGCUAGCUA"),
    mirna_seed_SEQ = Biostrings::RNAString("AUCGUAA")
  )
  
  expect_s4_class(mir, "miRNAGene")
  expect_equal(mir@geneID, "ENST000001")
  expect_equal(mir@h_symbol, "ACTB")
  expect_equal(mir@geneName, "Beta Actin")
  expect_equal(mir@description, "A major component of the cytoskeleton")
  expect_true(inherits(mir@geneStructure, "GRanges"))
  expect_equal(mir@RnaID, "ENST_RNA_001")
  expect_equal(as.character(mir@RnaSEQ), "AUGGCUAGCUAGCUAGCUA")
  expect_equal(as.character(mir@mirna_seed_SEQ), "AUCGUAA")
})


# test constructor with parameters missing:
# geneID              PRESENT
# h_symbol            PRESENT
# geneName            MISSING
# description         MISSING
# geneStructure       MISSING
# RnaID               MISSING
# RnaSEQ              MISSING
# mirna_seed_SEQ      MISSING


test_that("miRNAGene constructor assigns defaults and returns warnings for missing parameters", {
  
  # expected warnings
  expect_warning(
    mir <- miRNAGene(
      geneID = "ENST000001",
      h_symbol = "ACTB"
    ),
    regexp = "Warning! Missing parameter\\(s\\): geneName, description, geneStructure, RnaID, RnaSEQ, mirna_seed_SEQ"
  )
  
  
  # expected attributes values
  expect_s4_class(mir, "miRNAGene")
  expect_equal(mir@geneID, "ENST000001")
  expect_equal(mir@h_symbol, "ACTB")
  expect_true(is.na(mir@geneName))
  expect_true(is.na(mir@description))
  expect_true(inherits(mir@geneStructure, "GRanges"))
  expect_true(is.na(mir@RnaID))
  expect_equal(as.character(mir@RnaSEQ), "")
  expect_equal(as.character(mir@mirna_seed_SEQ), "")
})



