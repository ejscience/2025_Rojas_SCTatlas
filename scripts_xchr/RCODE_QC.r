#!/usr/bin/env Rscript
# RCODE_QC_RNA_only.r
# This script performs QC on RNA SNP calls by merging raw SNP data with ANNOVAR annotation.
# It is modified for RNA-only data and does not filter based on allele frequency.
# The output file name is generated based on the input file name.

library(tidyverse)
library(viridis)
library(patchwork)
library(tools)  # for file_path_sans_ext()

# Get command-line arguments (expecting two arguments: raw RNA SNP summary and RNA multianno file)
argv <- commandArgs(trailingOnly = TRUE)
if(length(argv) < 2) {
  stop("Usage: Rscript RCODE_QC_RNA_only.r <RNA_SNPCELL> <RNA_ANNOVAR>")
}

RNA_SNPCELL <- argv[1]
RNA_ANNOVAR <- argv[2]

# Load raw SNP call data (should contain a pre-computed SNP_ID column)
rna <- read_tsv(RNA_SNPCELL)

# Function to load annotation data from ANNOVAR multianno file
load_anno <- function(ANNOVAR) {
  anno <- read_tsv(ANNOVAR) %>%
    mutate(SNP_ID = str_c(Chr, Start, Ref, Alt, sep=":"))
  
  # If the ALL.sites.2015_08 column is missing, create it with default 0
  if(!"ALL.sites.2015_08" %in% colnames(anno)) {
    anno <- anno %>% mutate(ALL.sites.2015_08 = 0)
  }
  
  # Select and rename relevant columns
  anno <- anno %>% dplyr::select(SNP_ID, Func.refGene, Gene.refGene, ALL.sites.2015_08)
  colnames(anno) <- c("SNP_ID", "Region", "Gene", "ALL_Freq")
  anno$ALL_Freq[is.na(anno$ALL_Freq)] <- 0
  return(anno)
}

rna_anno <- load_anno(RNA_ANNOVAR)

# Merge raw data with annotation data
merge_anno <- function(data, anno) {
  df <- left_join(data, anno, by = "SNP_ID") %>%
    # Remove variants with multiple gene assignments
    filter(!str_detect(Gene, pattern = ";")) %>%
    filter(!is.na(Gene)) %>%
    # Keep variants in regions of interest (adjust as needed)
    filter(Region %in% c("intergenic", "intronic", "UTR5", "UTR3", "exonic", 
                          "ncRNA_exonic", "ncRNA_intronic", "splicing"))
  return(df)
}

rna_df <- merge_anno(rna, rna_anno)

# For RNA-only, we no longer filter on allele frequency since ALL_Freq is 0 for most variants
rna_df_QC_pass <- rna_df

# Generate an output file name based on the RNA input file name
output_file <- paste0(file_path_sans_ext(RNA_SNPCELL), "_QC_passed_SNP_df.tsv.gz")

# Write the QC-passed SNP data to the output file
write_tsv(rna_df_QC_pass, output_file)
