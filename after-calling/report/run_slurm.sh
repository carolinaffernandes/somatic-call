#!/usr/bin/env bash

#SBATCH --job-name=Report
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=log/report_%j.out  
#SBATCH --error=log/report_%j.err

LOG_DIR="log"
mkdir -p "$LOG_DIR"

# === EXECUÇÃO ===
echo "Iniciando em: $(date)"
python3 report.py
echo "Finalizado em: $(date)"
