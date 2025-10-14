#!/bin/bash
#SBATCH --job-name=FilterMutectCalls
#SBATCH --output=logs/FilterMutectCalls_%A_%a.out
#SBATCH --error=logs/FilterMutectCalls_%A_%a.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00

CONFIG="$PWD/config.yaml"
SAMPLESHEET="$PWD/samplesheet.tsv"
CHROM_FILE="$PWD/chrom_list.txt"

get_yaml() { python3 get_yaml.py "$CONFIG" "$1"; }

OUTPUT_DIR=$(get_yaml output_dir)
FILTERED_DIR=$(get_yaml filtered_dir)
REFERENCE=$(get_yaml reference)
SINGULARITY_IMG=$(get_yaml singularity_img)

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$SAMPLESHEET")
SAMPLE=$(echo "$LINE" | cut -f1)
CHROM=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$CHROM_FILE")

MUTECT2_OUT="${OUTPUT_DIR}/${SAMPLE}_${CHROM}.vcf.gz"
FILTERED_OUT="${FILTERED_DIR}/${SAMPLE}_${CHROM}_filtered.vcf.gz"
LOG_FILE="logs/${SAMPLE}_FilterMutectCalls_${CHROM}.log"

mkdir -p "$FILTERED_DIR" "logs"

echo ">> $SAMPLE - FilterMutectCalls - $CHROM" | tee -a "$LOG_FILE"
singularity exec "$SINGULARITY_IMG" gatk FilterMutectCalls \
    -R "$REFERENCE" \
    -V "$MUTECT2_OUT" \
    -O "$FILTERED_OUT" \
    &>> "$LOG_FILE"

if [[ $? -ne 0 ]]; then
    echo "ERRO no FilterMutectCalls para $SAMPLE - $CHROM" | tee -a "$LOG_FILE" >&2
    exit 1
fi

