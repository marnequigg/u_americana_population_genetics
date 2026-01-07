#!/bin/bash
set -euo pipefail

# 📁 Reference genome
REF_GENOME="/home/imh4101/genomes/masked_elm/masked_genome_elm.fa"

# 📁 Input and output directories
RG_DIR="/home/imh4101/gatk_again/diploids/01.rg_data"
GVCF_DIR="/home/imh4101/gatk_again/diploids/02.gvcfs"
STATUS_LOG="${GVCF_DIR}/gatk2n_pipeline_status.log"

mkdir -p "$GVCF_DIR"

# 📝 Initialize status log
echo -e "SampleID\tBatch\tStatus\tMessage" > "$STATUS_LOG"

# 📋 Logging function
log_status() {
  local sample=$1
  local batch=$2
  local status=$3
  local message=$4
  echo -e "${sample}\t${batch}\t${status}\t${message}" >> "$STATUS_LOG"
}
export -f log_status

# 🧬 Sample processing function
process_sample() {
  BAM="$1"
  filename=$(basename "$BAM")
  SAMPLE_ID=$(echo "$filename" | cut -d'_' -f1)
  BATCH=$(echo "$filename" | grep -o 'b[23]')
  OUTPUT_VCF="${GVCF_DIR}/${SAMPLE_ID}.g.vcf.gz"

  # 🔍 Index BAM
BAI="${BAM}.bai"
if [ ! -f "$BAI" ]; then
  if ! samtools index "$BAM" -@ 15; then
    log_status "$SAMPLE_ID" "$BATCH" "FAILED" "Indexing step"
    exit 2
  fi
fi

  # 🧬 Run GATK HaplotypeCaller
  if ! gatk HaplotypeCaller \
    -R "$REF_GENOME" \
    -I "$BAM" \
    -O "$OUTPUT_VCF" \
    -ERC GVCF \
    --sample-ploidy 2 \
    --native-pair-hmm-threads 15; then
      log_status "$SAMPLE_ID" "$BATCH" "FAILED" "GATK HaplotypeCaller"
      exit 3
  fi

  # ✅ Success
  log_status "$SAMPLE_ID" "$BATCH" "SUCCESS" "GATK completed"
}
export -f process_sample
export REF_GENOME GVCF_DIR STATUS_LOG

# 🚀 Run in parallel (2 jobs at a time)
find "$RG_DIR" -maxdepth 1 -name "*_b[23]_rg.bam" | parallel -j 2 process_sample {}
