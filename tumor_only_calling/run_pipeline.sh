#!/usr/bin/env bash
#SBATCH --job-name=run_pipeline
#SBATCH --output=logs/run_pipeline_%j.out
#SBATCH --error=logs/run_pipeline_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

set -euo pipefail

CONFIG="$PWD/config.yaml"
SAMPLESHEET="$PWD/samplesheet.tsv"
CHROM_FILE="$PWD/chrom_list.txt"
TASK_LIST="$PWD/task_list.txt"

OUTPUT_DIR=$(python3 get_yaml.py "$CONFIG" "output_dir")
FILTERED_DIR=$(python3 get_yaml.py "$CONFIG" "filtered_dir")
TMP_DIR=$(python3 get_yaml.py "$CONFIG" "tmp_dir")
REFERENCE=$(python3 get_yaml.py "$CONFIG" "reference")
GERMLINE=$(python3 get_yaml.py "$CONFIG" "germline_resource")
PON=$(python3 get_yaml.py "$CONFIG" "pon")
SINGULARITY_IMG=$(python3 get_yaml.py "$CONFIG" "singularity_img")
MAX_JOBS=$(python3 get_yaml.py "$CONFIG" "max_jobs")

mkdir -p "$OUTPUT_DIR" "$FILTERED_DIR" "$TMP_DIR" "logs"

echo "=== CONFIGURAÇÕES CARREGADAS ==="
echo "OUTPUT_DIR: $OUTPUT_DIR"
echo "FILTERED_DIR: $FILTERED_DIR"
echo "TMP_DIR: $TMP_DIR"
echo "REFERENCE: $REFERENCE"
echo "GERMLINE: $GERMLINE"
echo "PON: $PON"
echo "SINGULARITY_IMG: $SINGULARITY_IMG"
echo "MAX_JOBS: $MAX_JOBS"
echo "==============================="

echo ">> Gerando lista de tarefas (amostra × cromossomo)..."
> "$TASK_LIST"
tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r SAMPLE BAM_PATH; do
  for CHROM in $(cat "$CHROM_FILE"); do
    echo -e "${SAMPLE}\t${BAM_PATH}\t${CHROM}" >> "$TASK_LIST"
  done
done

TOTAL_TASKS=$(wc -l < "$TASK_LIST")
echo "Arquivo $TASK_LIST criado com $TOTAL_TASKS combinações."
echo ">> Submetendo Mutect2 em blocos de ${MAX_JOBS}..."

# Função para rodar sbatch em blocos
submit_in_chunks() {
  local script="$1"
  local total="$2"
  local step="$3"
  local dep="$4"
  local chunk_size="$MAX_JOBS"

  local start=1
  while [[ $start -le $total ]]; do
    local end=$((start + chunk_size - 1))
    [[ $end -gt $total ]] && end=$total

    if [[ -z "$dep" ]]; then
      JOBID=$(sbatch --parsable --array=${start}-${end}%${chunk_size} "$script")
    else
      JOBID=$(sbatch --parsable --dependency=afterok:${dep} --array=${start}-${end}%${chunk_size} "$script")
    fi

    echo "  -> Submetido bloco ${start}-${end} ($step) JOBID=$JOBID"
    dep=$JOBID

    # Aguarda os jobs do bloco terminarem antes do próximo
    echo "Aguardando bloco ${start}-${end} de ${step} finalizar..."
    squeue --job $JOBID > /dev/null 2>&1
    while squeue -j "$JOBID" 2>/dev/null | grep -q "$USER"; do
      sleep 60
    done

    start=$((end + 1))
  done

  echo "$dep"
}

# 1️⃣ Mutect2
JOBID_MUTECT2=$(submit_in_chunks "scripts/mutect2_chr.sh" "$TOTAL_TASKS" "Mutect2" "")

# 2️⃣ FilterMutectCalls (dependente)
JOBID_FMC=$(submit_in_chunks "scripts/filtermutectcalls_chr.sh" "$TOTAL_TASKS" "FilterMutectCalls" "$JOBID_MUTECT2")

# 3️⃣ MergeVCFs (depende de FilterMutectCalls)
JOBID_MERGE=$(sbatch --parsable --dependency=afterok:${JOBID_FMC} scripts/merge_vcfs.sh)

echo "=============================================="
echo "Pipeline submetida com sucesso!"
echo "Mutect2 JOBID final: $JOBID_MUTECT2"
echo "FilterMutectCalls JOBID final: $JOBID_FMC"
echo "MergeVCFs JOBID: $JOBID_MERGE"
echo "=============================================="

