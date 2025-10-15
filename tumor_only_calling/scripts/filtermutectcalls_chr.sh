#!/usr/bin/env bash
#SBATCH --job-name=FilterMutectCalls
#SBATCH --output=logs/FilterMutectCalls_%A_%a.out
#SBATCH --error=logs/FilterMutectCalls_%A_%a.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00

set -euo pipefail

CONFIG="$PWD/config.yaml"
TASK_LIST="$PWD/task_list.txt"

get_yaml() { python3 get_yaml.py "$CONFIG" "$1"; }

OUTPUT_DIR=$(get_yaml output_dir)
FILTERED_DIR=$(get_yaml filtered_dir)
REFERENCE=$(get_yaml reference)
SINGULARITY_IMG=$(get_yaml singularity_img)

TASK_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASK_LIST")
SAMPLE=$(echo "$TASK_LINE" | cut -f1)
CHROM=$(echo "$TASK_LINE" | cut -f3)

MUTECT2_OUT="${OUTPUT_DIR}/${SAMPLE}_${CHROM}.vcf.gz"
FILTERED_OUT="${FILTERED_DIR}/${SAMPLE}_${CHROM}_filtered.vcf.gz"
LOG_DIR="logs"
LOG_FILE="${LOG_DIR}/${SAMPLE}_FilterMutectCalls_${CHROM}.log"

mkdir -p "$FILTERED_DIR" "$LOG_DIR"

# verifica se já filtrado
if [[ -s "$FILTERED_OUT" ]]; then
  echo "[$(date)] SKIP: ${FILTERED_OUT} already exists — skipping FilterMutectCalls." | tee -a "$LOG_FILE"
  exit 0
fi

# verifica se o VCF de entrada existe e está pronto
if [[ ! -s "$MUTECT2_OUT" ]]; then
  echo "[$(date)] ERROR: input VCF ${MUTECT2_OUT} not found or empty — cannot run FilterMutectCalls." | tee -a "$LOG_FILE" >&2
  exit 2
fi

echo ">> $SAMPLE - FilterMutectCalls - $CHROM" | tee -a "$LOG_FILE"
singularity exec "$SINGULARITY_IMG" gatk FilterMutectCalls \
    -R "$REFERENCE" \
    -V "$MUTECT2_OUT" \
    -O "$FILTERED_OUT" \
    &>> "$LOG_FILE"

if [[ $? -ne 0 ]]; then
    echo "[$(date)] ERRO no FilterMutectCalls para $SAMPLE - $CHROM" | tee -a "$LOG_FILE" >&2
    exit 1
fi

echo "DONE FILTERMUTECTCALLS" >> "$LOG_FILE"
echo "[$(date)] FilterMutectCalls finished for ${SAMPLE}_${CHROM}" | tee -a "$LOG_FILE"
