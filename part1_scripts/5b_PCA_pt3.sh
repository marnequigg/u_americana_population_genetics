#!/bin/bash

SOURCE_DIR="/home/imh4101/all_marked_bams/rg_bams"

declare -A csv_to_folder=(
    ["cluster1_4n_final.csv"]="/home/imh4101/gatk_again/tetraploids"
    ["cluster2_other_final.csv"]="/home/imh4101/gatk_again/other"
    ["cluster3_2n_final.csv"]="/home/imh4101/gatk_again/diploids"
)

for csv_file in "${!csv_to_folder[@]}"; do
    dest_folder="${csv_to_folder[$csv_file]}"
    mkdir -p "$dest_folder"

    while IFS= read -r sample_id; do
        # Use find to locate matching BAM files
        found_files=$(find "$SOURCE_DIR" -type f -name "${sample_id}_b*_rg.bam")
        
        if [[ -z "$found_files" ]]; then
            echo "Warning: No BAM file found for sample ID $sample_id"
        else
            while IFS= read -r bam_file; do
                cp "$bam_file" "$dest_folder"
                echo "Copied $bam_file to $dest_folder"
            done <<< "$found_files"
        fi
    done < "$csv_file"
done
