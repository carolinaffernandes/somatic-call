#!/bin/bash
#SBATCH --job-name=MergeVCFs
#SBATCH --output=logs/MergeVCFs_%j.out
#SBATCH --error=logs/MergeVCFs_%j.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00

CONFIG="$PWD/config.yaml"
SAMPLESHEET="$PWD/samplesheet.tsv"

get_yaml() { python3 get_yaml.py "$CONFIG" "$1"; }

FILTERED_DIR=$(get_yaml filtered_dir)
REFERENCE=$(get_yaml reference)
SINGULARITY_IMG=$(get_yaml singularity_img)
OUTPUT_MERGED="${FILTERED_DIR}/all_samples_merged.vcf.gz"

VCF_LIST=$(ls "$FILTERED_DIR"/*_filtered.vcf.gz | tr '\n' ' ')

echo ">> MergeVCFs - Todos os filtros aplicados"
singularity exec "$SINGULARITY_IMG" java -jar /gatk/picard.jar MergeVcfs \
    I=$VCF_LIST \
    O="$OUTPUT_MERGED" \
    &> logs/MergeVCFs.log

if [[ $? -ne 0 ]]; then
    echo "ERRO no MergeVCFs" | tee -a logs/MergeVCFs.log >&2
    exit 1
fi

echo "MergeVCFs concluído: $OUTPUT_MERGED"

