#!/usr/bin/env Rscript
# RCODE_process_cellsnp.r
# This script processes the output from cellsnp-lite for allele-specific expression.
# It reads in the sample barcode file, the AD, DP, and OTH matrices, and the VCF file,
# then computes the reference counts (DP - AD), pivots the matrices to long format,
# merges the data with SNP information from the VCF, and writes the final table.
#
# Usage:
#   Rscript RCODE_process_cellsnp.r <SAMPLE_NAME> <ADmtx> <DPmtx> <OTHmtx> <VCF> <OUT>

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(Matrix)  # For readMM()

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 6) {
  stop("Usage: Rscript RCODE_process_cellsnp.r <SAMPLE_NAME> <ADmtx> <DPmtx> <OTHmtx> <VCF> <OUT>")
}

SAMPLE_NAME <- args[1]
ADmtx      <- args[2]
DPmtx      <- args[3]
OTHmtx     <- args[4]
VCF        <- args[5]
OUT        <- args[6]

# Read the sample barcode file (without header)
sample_name <- read_tsv(SAMPLE_NAME, col_names = FALSE)

# Read the count matrices
AD  <- Matrix::readMM(ADmtx)
DP  <- Matrix::readMM(DPmtx)
OTH <- Matrix::readMM(OTHmtx)
# Compute the reference count matrix
RD  <- DP - AD

# Read the VCF file, skipping header lines starting with "#"
vcf <- read_tsv(VCF, comment = "#", col_names = FALSE)

# Build SNP information from the VCF.
# Creating a SNP ID as "chr:pos:REF:ALT" and selecting relevant columns.
SNP_INFO <- vcf %>%
  mutate(SNP_ID = str_c(X1, ":", X2, ":", X4, ":", X5)) %>%
  select(SNP_ID, X1, X2, X4, X5)
colnames(SNP_INFO) <- c("SNP_ID", "CHR", "POS", "REF", "ALT")

# Convert the matrices to data frames and attach SNP IDs.
alt_res <- AD %>% as.matrix() %>% as_tibble()
colnames(alt_res) <- sample_name$X1
alt_res$SNP <- SNP_INFO$SNP_ID

ref_res <- RD %>% as.matrix() %>% as_tibble()
colnames(ref_res) <- sample_name$X1
ref_res$SNP <- SNP_INFO$SNP_ID

oth_res <- OTH %>% as.matrix() %>% as_tibble()
colnames(oth_res) <- sample_name$X1
oth_res$SNP <- SNP_INFO$SNP_ID

# Pivot each matrix to long format and filter to keep only positive counts
alt_res <- alt_res %>%
  pivot_longer(-SNP, names_to = "cell_barcode", values_to = "ALTcount") %>%
  filter(ALTcount > 0)

ref_res <- ref_res %>%
  pivot_longer(-SNP, names_to = "cell_barcode", values_to = "REFcount") %>%
  filter(REFcount > 0)

oth_res <- oth_res %>%
  pivot_longer(-SNP, names_to = "cell_barcode", values_to = "OTHcount") %>%
  filter(OTHcount > 0)

# Merge the allele count data
result <- full_join(ref_res, alt_res, by = c("SNP", "cell_barcode")) %>%
  full_join(oth_res, by = c("SNP", "cell_barcode"))

# Merge in the SNP information using the SNP_ID
result <- left_join(SNP_INFO, result, by = c("SNP_ID" = "SNP"))
result[is.na(result)] <- 0
result <- distinct(result)

# Write the final processed table to the specified output file
write_tsv(result, OUT)
