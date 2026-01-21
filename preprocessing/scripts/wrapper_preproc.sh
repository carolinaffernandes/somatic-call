#!/bin/bash
#SBATCH --job-name=Worker_PreProc
#SBATCH --output=logs/worker_%A_%a.out
#SBATCH --error=logs/worker_%A_%a.err

set -e 

TASK_FILE=$1
CONFIG_FILE=$2

# Funcao leitura YAML
get_config() {
    python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE'))['$1'])"
}

# Identifica tarefa
LINE_NUM=$SLURM_ARRAY_TASK_ID
TASK_DATA=$(sed "${LINE_NUM}q;d" "$TASK_FILE")

SAMPLE=$(echo "$TASK_DATA" | awk '{print $1}')
R1=$(echo "$TASK_DATA" | awk '{print $2}')
R2=$(echo "$TASK_DATA" | awk '{print $3}')

echo "=== INICIANDO PRE-PROCESSAMENTO ==="
echo "Amostra: $SAMPLE"
echo "Task ID: $SLURM_ARRAY_TASK_ID"

# Carrega configs
REF=$(get_config reference_genome)
SIF=$(get_config singularity_image)
DBSNP=$(get_config known_sites_snp)
INDELS=$(get_config known_sites_indels)
THREADS=$SLURM_CPUS_PER_TASK
RESULTS_DIR=$(get_config results_dir)

# Define Diretorios
DIR_BWA="$RESULTS_DIR/bwa"
DIR_MARKDUP="$RESULTS_DIR/markduplicates"
DIR_BQSR="$RESULTS_DIR/bqsr"

# 1. BWA
echo "[Step 1] Alinhamento BWA..."
singularity exec -B /mnt:/mnt "$SIF" /bin/bash -c \
    "bwa mem -t $THREADS -R '@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA' $REF $R1 $R2 | \
    samtools sort -@ $THREADS -o $DIR_BWA/$SAMPLE.sorted.bam -"

# 2. MarkDuplicates
echo "[Step 2] MarkDuplicates..."
singularity exec -B /mnt:/mnt "$SIF" gatk MarkDuplicates \
    -I "$DIR_BWA/$SAMPLE.sorted.bam" \
    -O "$DIR_MARKDUP/$SAMPLE.dedup.bam" \
    -M "$DIR_MARKDUP/$SAMPLE.metrics.txt" \
    --CREATE_INDEX true

# 3. BaseRecalibrator
echo "[Step 3] BaseRecalibrator..."
singularity exec -B /mnt:/mnt "$SIF" gatk BaseRecalibrator \
    -I "$DIR_MARKDUP/$SAMPLE.dedup.bam" -R "$REF" \
    --known-sites "$DBSNP" --known-sites "$INDELS" \
    -O "$DIR_BQSR/$SAMPLE.recal_data.table"

# 4. ApplyBQSR
echo "[Step 4] ApplyBQSR..."
singularity exec -B /mnt:/mnt "$SIF" gatk ApplyBQSR \
    -I "$DIR_MARKDUP/$SAMPLE.dedup.bam" -R "$REF" \
    --bqsr-recal-file "$DIR_BQSR/$SAMPLE.recal_data.table" \
    -O "$DIR_BQSR/$SAMPLE.final.bam"

echo "SUCESSO: Amostra $SAMPLE finalizada."
