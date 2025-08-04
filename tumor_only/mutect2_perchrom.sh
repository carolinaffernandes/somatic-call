#!/bin/bash

set -euo pipefail

CONFIG_FILE="config.yaml"

read_config() {
  local key="$1"
  python3 -c "
import yaml
with open('$CONFIG_FILE') as f:
    cfg = yaml.safe_load(f)
keys = '$key'.split('.')
val = cfg
for k in keys:
    val = val.get(k, None)
    if val is None:
        break
print(val if val else '')
" 2>/dev/null
}

ENV_NAME=$(read_config conda_env)

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

REFERENCE=$(read_config reference)
GERMLINE_RESOURCE=$(read_config germline_resource)
OUT_DIR=$(read_config output_dir)
TMP_DIR=$(read_config tmp_dir)
SAMPLE_SHEET=$(read_config sample_sheet)

mkdir -p "$OUT_DIR/raw" "$OUT_DIR/filtered" "$TMP_DIR" logs

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLE_SHEET")
SAMPLE=$(echo "$LINE" | cut -f1)
BAM=$(echo "$LINE" | cut -f2)

RAW_VCF="${OUT_DIR}/raw/${SAMPLE}.vcf.gz"
FILTERED_VCF="${OUT_DIR}/filtered/${SAMPLE}.filtered.vcf.gz"

echo "Processando $SAMPLE"
echo "BAM: $BAM"

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

echo "Finalizado $SAMPLE"

