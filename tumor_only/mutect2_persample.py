#!/bin/bash

CONFIG_FILE="config.yaml"

if ! command -v yq &> /dev/null; then
    echo "Erro: yq não encontrado. Instale com 'conda install -c conda-forge yq' ou 'mamba install yq'."
    exit 1
fi

CONDA_ENV=$(yq '.conda_env' "$CONFIG_FILE")
REFERENCE=$(yq '.reference' "$CONFIG_FILE")
GERMLINE_RESOURCE=$(yq '.germline_resource' "$CONFIG_FILE")
OUT_DIR=$(yq '.output_dir' "$CONFIG_FILE")
TMP_DIR=$(yq '.tmp_dir' "$CONFIG_FILE")
SAMPLE_SHEET=$(yq '.sample_sheet' "$CONFIG_FILE")

mkdir -p "$OUT_DIR/raw" "$OUT_DIR/filtered" "$TMP_DIR" logs

source ~/miniconda3/etc/profile.d/conda.sh
conda activate "$CONDA_ENV"

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLE_SHEET")
SAMPLE=$(echo "$LINE" | cut -f1)
BAM=$(echo "$LINE" | cut -f2)

RAW_VCF="${OUT_DIR}/raw/${SAMPLE}.vcf.gz"
FILTERED_VCF="${OUT_DIR}/filtered/${SAMPLE}.filtered.vcf.gz"

gatk Mutect2 \
    -R "$REFERENCE" \
    -I "$BAM" \
    --germline-resource "$GERMLINE_RESOURCE" \
    -O "$RAW_VCF" \
    --tmp-dir "$TMP_DIR"

gatk FilterMutectCalls \
    -V "$RAW_VCF" \
    -O "$FILTERED_VCF" \
    --tmp-dir "$TMP_DIR"

echo "Pipeline concluído para $SAMPLE"
