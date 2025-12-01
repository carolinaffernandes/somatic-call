#!/bin/bash
#SBATCH --job-name=sanity
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH -o qc_%j.out
#SBATCH -e qc_%j.err

python qc_panel.py

