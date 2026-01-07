#this starts with bwa alignment and ends with primerclip

#download the softwares
conda install bioconda::bwa #version 0.7.17
conda install bioconda::samtools #version 1.5
conda install bioconda::primerclip #version 0.3.8

#Dont forgot to index
#only need to do this once
bwa index /data/labs/Fant/Quigg/genomes/masked_genome_elm.fa

#start your ENGINESSSSSSSS
nano pipeline.sh

#!/bin/bash

### Step 1: BWA MEM Alignment ###

FASTQ_DIR="/data/labs/Fant/Quigg/02.trimming"
REFERENCE="/data/labs/Fant/Quigg/genomes/masked_genome_elm.fa"
ALIGNMENT_DIR="/data/labs/Fant/Quigg/03.alignment_clip"

# Create output directory if needed
mkdir -p "$ALIGNMENT_DIR"

# Align all paired-end FASTQ files with BWA MEM
for R1 in "$FASTQ_DIR"/*_R1_paired.fastq.gz; do
    SAMPLE=$(basename "$R1" _R1_paired.fastq.gz)
    R2="$FASTQ_DIR/${SAMPLE}_R2_paired.fastq.gz"
    SAM_OUTPUT="$ALIGNMENT_DIR/${SAMPLE}.sam"

    echo "Aligning $SAMPLE..."
    bwa mem -M -t 20 "$REFERENCE" "$R1" "$R2" > "$SAM_OUTPUT"
    echo "Alignment complete for $SAMPLE!"
done

echo "All alignments completed!"

### Step 2: Convert SAM to BAM ###

BAM_DIR="/data/labs/Fant/Quigg/03.alignment_clip/bams"
mkdir -p "$BAM_DIR"

for SAM_FILE in "$ALIGNMENT_DIR"/*.sam; do
    SAMPLE=$(basename "$SAM_FILE" .sam)
    BAM_OUTPUT="$BAM_DIR/${SAMPLE}.bam"

    echo "Converting $SAM_FILE to BAM..."
    samtools view -b "$SAM_FILE" -o "$BAM_OUTPUT"
    echo "Conversion completed for $SAMPLE!"
done

### Step 3: Sort BAM & Generate Flagstats ###
SORTED_DIR="/data/labs/Fant/Quigg/03.alignment_clip/sorted_bams/"
FLAGSTAT_DIR="/data/labs/Fant/Quigg/03.alignment_clip/flagstats/"
mkdir -p "$SORTED_DIR" "$FLAGSTAT_DIR"

for BAM_FILE in "$BAM_DIR"/*.bam; do
    SAMPLE=$(basename "$BAM_FILE" .bam)
    SORTED_BAM="$SORTED_DIR/${SAMPLE}_sorted.bam"
    
    echo "Sorting BAM file for $SAMPLE..."
    samtools sort -o "$SORTED_BAM" "$BAM_FILE"

    echo "Generating flagstat for $SAMPLE..."
    samtools flagstat "$SORTED_BAM" > "$FLAGSTAT_DIR/${SAMPLE}_flagstat.txt"
done

echo "Sorting and flagstat analysis completed!"

### Step 5: Sort SAM Files by Name ###
INDS=($(for i in "$ALIGNMENT_DIR"/*.sam; do echo $(basename "$i" .sam); done))

for IND in "${INDS[@]}"; do
    samtools sort -n "$ALIGNMENT_DIR/${IND}.sam" -o "$ALIGNMENT_DIR/${IND}_nsort.sam"
done

echo "All SAM sorting completed!"

### Step 6: Run PrimerClip ###

PRIMERFILE="/data/labs/Fant/Quigg/genomes/cp927_Masterfile.txt"

for IND in "${INDS[@]}"; do
    echo "Running PrimerClip on $IND..."
    primerclip "$PRIMERFILE" "$ALIGNMENT_DIR/${IND}_nsort.sam" "$ALIGNMENT_DIR/${IND}_clipped.sam"
done

echo "Primer clipping completed!"
#exit

chmod +x pipeline.sh

#LETS GOOOOOOOOOOOOOOOOOO
screen -L ./pipeline.sh

#to disconnect the screen (keep it running in the background)
#type this into the screen itself
ctrl + a
d
