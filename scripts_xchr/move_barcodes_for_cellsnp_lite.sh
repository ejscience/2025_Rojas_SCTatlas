#!/bin/bash

samples=("SCT01_nuc" "SCT02" "SCT03" "SCT06" "SCT09" "SCT10")
base_dir="/tank/data/ernesto/data/scrna/cellranger_outs"
out_dir="/tank/data/ernesto/analysis/Xchr_SCT_AlleleSpecific"

for sample in "${samples[@]}"; do
    barcode_file="${base_dir}/${sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
    out_file="${out_dir}/${sample}_RNA.barcode.list.tsv"
    
    # Extract the barcodes without modifying them:
    # zcat "$barcode_file" > "$out_file"
    
    # Uncomment the following line if you need to change "-1" to "-RNA" in the barcodes:
    zcat "$barcode_file" | sed 's/-1/-RNA/g' > "$out_file"
    
    echo "Processed ${sample} barcodes into ${out_file}"
done


