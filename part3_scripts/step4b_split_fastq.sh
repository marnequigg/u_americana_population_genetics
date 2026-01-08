#!/bin/bash

# Set your input and output directories
INPUT_DIR="/data/labs/Fant/Quigg/03e.hybpiper/03.hybphaser1/05_phasing/phased_reads/diploid_sub_fastq"
OUTPUT_DIR="/data/labs/Fant/Quigg/03e.hybpiper/03.hybphaser1/05_phasing/phased_reads/diploid_sub_fastq_split"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all .fastq files in the input directory
for file in "$INPUT_DIR"/*.fastq; do
    # Extract sample name without extension
    base=$(basename "$file" .fastq)

    echo "🔄 Splitting $base..."

    # Run BBMap's reformat.sh to split interleaved FASTQ
    reformat.sh in="$file" \
                out1="$OUTPUT_DIR/${base}_R1.fastq" \
                out2="$OUTPUT_DIR/${base}_R2.fastq"

    echo "✅ Finished: $base → R1 and R2 saved"
done

echo "🎉 All interleaved FASTQs have been split."
