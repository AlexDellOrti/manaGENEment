## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=FALSE---------------------------------------------------------------
# Load the manaGENEment package
library(manaGENEment)

## ----echo=TRUE, collapse =FALSE-----------------------------------------------

# Create a ProteinCodingGene object
pcg <- ProteinCodingGene(
  geneID = "ENST000001",
  h_symbol = "ACTB",
  geneName = "Beta Actin",
  description = "One of six different actin isoforms and a major component of the cytoskeleton",
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


## ----echo=TRUE, collapse =FALSE-----------------------------------------------
# Create a lncRNAGene object
lnc <- lncRNAGene(
  geneID = "ENST000002",
  h_symbol = "HOTAIR",
  geneName = "HOX Transcript Antisense RNA",
  description = "Involved in epigenetic regulation",
  geneStructure = GenomicRanges::GRanges(
    seqnames = "chr12",
    ranges = IRanges::IRanges(start = 1150000, end = 1160000),
    strand = "-"
  ),
  RnaID = "ENST_RNA_002",
  RnaSEQ = Biostrings::RNAString("AUGGCUAGCUAGCUAGCUA")
)


## ----echo=TRUE, collapse =FALSE-----------------------------------------------
# Create a miRNAGene object
mir <- miRNAGene(
  geneID = "ENST000003",
  h_symbol = "MIR21",
  geneName = "miR-21",
  description = "Oncogenic miRNA involved in various cancers",
  geneStructure = GenomicRanges::GRanges(
    seqnames = "chr17",
    ranges = IRanges::IRanges(start = 43045629, end = 43045716),
    strand = "+"
  ),
  RnaID = "ENST_RNA_003",
  RnaSEQ = Biostrings::RNAString("AUGGCUAGCUAGCUAGCUA"),
  mirna_seed_SEQ = Biostrings::RNAString("UAGCUA")
)

## ----echo=TRUE, collapse =FALSE-----------------------------------------------
# Retrieve the Gene ID from the ProteinCodingGene object
gene_id_pcg <- getGeneID(pcg)
print(gene_id_pcg)

# Retrieve the Gene ID from the lncRNAGene object
gene_id_lnc <- getGeneID(lnc)
print(gene_id_lnc)

# Retrieve the Gene ID from the miRNAGene object
gene_id_mir <- getGeneID(mir)
print(gene_id_mir)


## ----echo=TRUE, collapse =FALSE-----------------------------------------------
# Retrieve the RNA ID from the ProteinCodingGene object
rna_id_pcg <- getRnaID(pcg)
print(rna_id_pcg)

# Retrieve the RNA ID from the lncRNAGene object
rna_id_lnc <- getRnaID(lnc)
print(rna_id_lnc)

# Retrieve the RNA ID from the miRNAGene object
rna_id_mir <- getRnaID(mir)
print(rna_id_mir)


## ----echo=TRUE, collapse =FALSE-----------------------------------------------
# Set a new Gene ID for the ProteinCodingGene object
pcg <- setGeneID(pcg, "ENST000004")
print(getGeneID(pcg))


## ----echo=TRUE, collapse =FALSE-----------------------------------------------
# Compute the length of the protein sequence
len_pcg <- lengthProduct(pcg)
print(len_pcg)  # Number of amino acids

## ----echo=TRUE, collapse =FALSE-----------------------------------------------
# Compute the length of the RNA sequence
len_lnc <- lengthProduct(lnc)
print(len_lnc)  # Length of RNA sequence

## ----echo=TRUE, collapse =FALSE-----------------------------------------------
# Compute the length of the RNA sequence
len_mir <- lengthProduct(mir)
print(len_mir)  # Length of RNA sequence

## ----echo=TRUE, collapse =FALSE-----------------------------------------------
# Retrieve the Protein ID and Sequence from the ProteinCodingGene object
prot_id <- getProtID(pcg)
prot_seq <- getProtSeq(pcg)
print(prot_id)
print(prot_seq)

## ----echo=TRUE, collapse =FALSE-----------------------------------------------
# Retrieve the miRNA seed sequence from the miRNAGene object
seed_seq <- getSeedSEQ(mir)
print(seed_seq)

