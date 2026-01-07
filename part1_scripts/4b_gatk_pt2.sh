#!/bin/bash
set -euo pipefail

# Define paths
INPUT_DIR="/home/imh4101/gvcfs/"
REFERENCE="/home/imh4101/genomes/masked_elm/masked_genome_elm.fa"
COHORT_GVCF="${INPUT_DIR}/cohort.g.vcf.gz"
FINAL_VCF="${INPUT_DIR}/final.vcf.gz"

# Step 1: Combine GVCFs
echo "🔗 Combining GVCFs from $INPUT_DIR..."

GVCF_LIST=""
for gvcf in "$INPUT_DIR"/*.g.vcf.gz; do
  [[ "$gvcf" == "$COHORT_GVCF" || "$gvcf" == "$FINAL_VCF" ]] && continue
  GVCF_LIST+=" --variant $gvcf"
done

gatk --java-options "-Xmx16g" CombineGVCFs \
  -R "$REFERENCE" \
  $GVCF_LIST \
  -O "$COHORT_GVCF"

echo "✅ GVCFs combined into: $COHORT_GVCF"

# Step 2: Genotype GVCFs
echo "🧬 Running GenotypeGVCFs..."

gatk --java-options "-Xmx24g" GenotypeGVCFs \
  -R "$REFERENCE" \
  -V "$COHORT_GVCF" \
  -O "$FINAL_VCF"

echo "✅ Final VCF saved to: $FINAL_VCF"
