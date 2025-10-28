#!/bin/bash
# List of samples to process
samples=("SCT02" "SCT03" "SCT06" "SCT09" "SCT10") # "SCT01_nuc")
    

# Path to the reference genome used in cellranger (update this path as needed)
refseq="/tank/data/ernesto/genomes/refdata-gex-GRCh38-2020-A/fasta/genome.fa"

# Number of threads to use
threads=24

# Base directories
bam_base="/tank/data/ernesto/data/scrna/cellranger_outs"
barcode_dir="/tank/data/ernesto/analysis/Xchr_SCT_AlleleSpecific"

for sample in "${samples[@]}"; do
    # BAM file is the X-chromosome specific file created earlier
    bam_file="${bam_base}/${sample}/outs/${sample}.X.bam"
    # Barcode file in our analysis folder (assumed not modified)
    barcode_file="${barcode_dir}/${sample}_RNA.barcode.list.tsv"
    # Output directory for cellsnp-lite results (will be created by cellsnp-lite)
    output_dir="${barcode_dir}/${sample}.chrX.snp.call"
    
    echo "Running cellsnp-lite for sample ${sample}"
    
    cellsnp-lite -s "$bam_file" \
        --minCOUNT 20 -b "$barcode_file" \
        --UMItag UB -p ${threads} \
        --minMAF 0.05 --minMAPQ 20 \
        --refseq "$refseq" \
        --chrom=chrX -O "$output_dir"
        
    echo "Finished cellsnp-lite for sample ${sample}"
done
