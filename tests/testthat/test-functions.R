
library(testthat)
library(manaGENEment)  
library(GenomicRanges)
library(Biostrings)
library(IRanges)



# the strategy to test all the functions is been organized like follows:
# - create a ProteinCodingGene object
# - create a lncRNAGene object
# - create a miRNAGene object
# 
# - test getGeneID on the ProteinCodingGene object
# - test setGeneID on the ProteinCodingGene object
# - test getRnaID on the ProteinCodingGene object
# - test getRnaSEQ on the ProteinCodingGene object
# - test getProtID on the ProteinCodingGene object
# - test getProtSeq on the ProteinCodingGene object
# - test lengthProduct on the ProteinCodingGene object
# 
# - test getGeneID on the lncRNAGene object
# - test setGeneID on the lncRNAGene object
# - test getRnaID on the lncRNAGene object
# - test getRnaSEQ on the lncRNAGene object
# - test lengthProduct on the lncRNAGene object
# 
# - test getGeneID on the miRNAGene object
# - test setGeneID on the miRNAGene object
# - test getRnaID on the miRNAGene object
# - test getRnaSEQ on the miRNAGene object
# - test getSeedSEQ on the miRNAGene object
# - test lengthProduct on the miRNAGene object

# - test errors when no parameters are given






# =====================================================================
#   ***                  Test Objects creation                  ***
# =====================================================================

# in order to test all the methods of each class, an object of each class is created, 
# on which all the test will be applied.

# Create a ProteinCodingGene object
pcg <- ProteinCodingGene(
  geneID = "ENST000001",
  h_symbol = "ACTB",
  geneName = "Beta Actin",
  description = "A major component of the cytoskeleton",
  geneStructure = GRanges(
    seqnames = "chr7",
    ranges = IRanges(start = 5527175, end = 5535161),
    strand = "+"
  ),
  RnaID = "ENST_RNA_001",
  RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA"),
  proteinID = "ENSP000001",
  proteinSEQ = AAString("MEEEIAALVIDNGSG")
)

# Create a lncRNAGene object
lnc <- lncRNAGene(
  geneID = "ENST000002",
  h_symbol = "ACTG",
  geneName = "Gamma Actin",
  description = "Gamma Actin gene description",
  geneStructure = GRanges(
    seqnames = "chr7",
    ranges = IRanges(start = 5527175, end = 5535161),
    strand = "+"
  ),
  RnaID = "ENST_RNA_002",
  RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA")
)

# Create a miRNAGene object
mir <- miRNAGene(
  geneID = "ENST000003",
  h_symbol = "ACTD",
  geneName = "Delta Actin",
  description = "Delta Actin gene description",
  geneStructure = GRanges(
    seqnames = "chr7",
    ranges = IRanges(start = 5527175, end = 5535161),
    strand = "+"
  ),
  RnaID = "ENST_RNA_003",
  RnaSEQ = RNAString("AUGGCUAGCUAGCUAGCUA"),
  mirna_seed_SEQ = RNAString("AUCGUAA")
)

# =====================================================================
#   ***               ProteinCodingGene Tests                  ***
# =====================================================================

test_that("getGeneID retrieves the correct Gene ID for ProteinCodingGene", {
  expect_equal(getGeneID(pcg), "ENST000001")
})

test_that("setGeneID updates the Gene ID for ProteinCodingGene", {
  updated_pcg <- setGeneID(pcg, "ENST000004")
  expect_equal(getGeneID(updated_pcg), "ENST000004")
})

test_that("getRnaID retrieves the correct RNA ID for ProteinCodingGene", {
  expect_equal(getRnaID(pcg), "ENST_RNA_001")
})

test_that("getRnaSEQ retrieves the correct RNA sequence for ProteinCodingGene", {
  expect_equal(as.character(getRnaSEQ(pcg)), "AUGGCUAGCUAGCUAGCUA")
})

test_that("getProtID retrieves the correct Protein ID for ProteinCodingGene", {
  expect_equal(getProtID(pcg), "ENSP000001")
})

test_that("getProtSeq retrieves the correct Protein Sequence for ProteinCodingGene", {
  expect_equal(as.character(getProtSeq(pcg)), "MEEEIAALVIDNGSG")
})

test_that("lengthProduct returns correct protein length for ProteinCodingGene", {
  expect_equal(lengthProduct(pcg), length(pcg@proteinSEQ))  # Should be 15
})

# =====================================================================
#   ***                lncRNAGene Tests                        ***
# =====================================================================

test_that("getGeneID retrieves the correct Gene ID for lncRNAGene", {
  expect_equal(getGeneID(lnc), "ENST000002")
})

test_that("setGeneID updates the Gene ID for lncRNAGene", {
  updated_lnc <- setGeneID(lnc, "ENST000006")
  expect_equal(getGeneID(updated_lnc), "ENST000006")
})

test_that("getRnaID retrieves the correct RNA ID for lncRNAGene", {
  expect_equal(getRnaID(lnc), "ENST_RNA_002")
})

test_that("getRnaSEQ retrieves the correct RNA sequence for lncRNAGene", {
  expect_equal(as.character(getRnaSEQ(lnc)), "AUGGCUAGCUAGCUAGCUA")
})

test_that("lengthProduct returns correct RNA length for lncRNAGene", {
  expect_equal(lengthProduct(lnc), length(lnc@RnaSEQ))  # Should be 18
})

# =====================================================================
#   ***                 miRNAGene Tests                          ***
# =====================================================================

test_that("getGeneID retrieves the correct Gene ID for miRNAGene", {
  expect_equal(getGeneID(mir), "ENST000003")
})

test_that("setGeneID updates the Gene ID for miRNAGene", {
  updated_mir <- setGeneID(mir, "ENST000007")
  expect_equal(getGeneID(updated_mir), "ENST000007")
})

test_that("getRnaID retrieves the correct RNA ID for miRNAGene", {
  expect_equal(getRnaID(mir), "ENST_RNA_003")
})

test_that("getRnaSEQ retrieves the correct RNA sequence for miRNAGene", {
  expect_equal(as.character(getRnaSEQ(mir)), "AUGGCUAGCUAGCUAGCUA")
})

test_that("getSeedSEQ retrieves the correct miRNA seed sequence for miRNAGene", {
  expect_equal(as.character(getSeedSEQ(mir)), "AUCGUAA")
})

test_that("lengthProduct returns correct RNA length for miRNAGene", {
  expect_equal(lengthProduct(mir), length(mir@RnaSEQ))  # Should be 18
})

# =====================================================================
#   ***           Setter Functions Error Handling                ***
# =====================================================================

test_that("setGeneID throws an error when newid is not a character string", {
  # ProteinCodingGene object
  expect_error(
    setGeneID(pcg, 12345),
    "The new ID you want to set has to be of class character."
  )
  
  # lncRNAGene object
  expect_error(
    setGeneID(lnc, TRUE),
    "The new ID you want to set has to be of class character."
  )
  
  # miRNAGene object
  expect_error(
    setGeneID(mir, list("ENST000008")),
    "The new ID you want to set has to be of class character."
  )
})


# =====================================================================
#   ***           Constructors Errors                          ***
# =====================================================================

test_that("All the Constructors throw an error when all parameters are missing", {
  expect_error(
    ProteinCodingGene(),
    "Insert at least one argument"
  )
  
  expect_error(
    lncRNAGene(),
    "Insert at least one argument"
  )
  
  expect_error(
    miRNAGene(),
    "Insert at least one argument"
  )
})