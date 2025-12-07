#!/bin/bash
#SBATCH --job-name=sanity
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH -o log/qc_%j.out
#SBATCH -e log/qc_%j.err

LOG_DIR="log"
mkdir -p "$LOG_DIR"

python qc_panel.py

