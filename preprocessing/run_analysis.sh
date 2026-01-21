#!/bin/bash
#SBATCH --job-name=Manager_PreProc
#SBATCH --output=logs/manager_%j.out
#SBATCH --error=logs/manager_%j.err
#SBATCH --time=5-00:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1

# CONFIGURACOES DO GERENTE
CHUNK_SIZE=10      # Quantas amostras enviar por vez
MAX_PARALLEL=10    # Quantas rodam simultaneamente dentro do lote
CONFIG_FILE="config.yaml"
TASK_LIST="preproc_task_list.txt"

echo "=== GERENTE DE PIPELINE (PRE-PROCESSING) ==="
echo "Data: $(date)"

mkdir -p logs

echo ">> [Manager] Gerando lista de tarefas..."
cd scripts
python3 create_tasklist.py
cd ..

if [ ! -f "$TASK_LIST" ]; then echo "ERRO: Lista de tarefas nao encontrada"; exit 1; fi

TOTAL_TASKS=$(wc -l < "$TASK_LIST")
echo "Total de amostras: $TOTAL_TASKS"

# Cria estrutura de pastas
RESULTS_DIR=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE'))['results_dir'])")
mkdir -p "$RESULTS_DIR/bwa" "$RESULTS_DIR/markduplicates" "$RESULTS_DIR/bqsr"

submit_in_chunks() {
    local script="$1"
    local total="$2"
    local start=1
    
    while [[ $start -le $total ]]; do
        local end=$((start + CHUNK_SIZE - 1))
        [[ $end -gt $total ]] && end=$total
        
        echo "---------------------------------------------------"
        echo ">> [Manager] Submetendo lote: $start ate $end..."
        
        # Pega recursos do config.yaml
        TIME=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE'))['slurm']['time'])")
        MEM=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE'))['slurm']['mem'])")
        CPUS=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE'))['slurm']['cpus_per_task'])")

        # Submete Array
        JOBID=$(sbatch --parsable \
            --array=${start}-${end}%${MAX_PARALLEL} \
            --time=$TIME \
            --mem=$MEM \
            --cpus-per-task=$CPUS \
            scripts/wrapper_preproc.sh "$TASK_LIST" "$CONFIG_FILE")
            
        echo "   Job ID lancado: $JOBID"
        echo "   [Manager] Aguardando lote terminar..."
        
        # Loop de espera (checa a cada 60s)
        while squeue -j "$JOBID" 2>/dev/null | grep -q "$USER"; do
            sleep 60
        done
        
        echo "   Lote finalizado."
        start=$((end + 1))
    done
}

# Executa
chmod +x scripts/wrapper_preproc.sh
submit_in_chunks "scripts/wrapper_preproc.sh" "$TOTAL_TASKS"

echo "PIPELINE FINALIZADO."
