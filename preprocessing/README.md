# Pipeline de Pré-processamento de Exoma (WES) 🧬

Este diretório contém os scripts para automação do pré-processamento de dados de sequenciamento de exoma (WES), seguindo as **GATK Best Practices**. O pipeline converte arquivos brutos **FastQ** em arquivos **BAM prontos para análise** (Analysis-Ready BAMs), adequados para chamada de variantes somáticas (Mutect2).

## 📋 Visão Geral do Fluxo

O pipeline executa os seguintes passos sequenciais para cada amostra:

1.  **Alinhamento (BWA-MEM):** Mapeamento das reads contra o genoma de referência (hg38).
2.  **Sort & Index:** Ordenação do BAM por coordenadas genômicas.
3.  **MarkDuplicates (GATK):** Identificação e marcação de duplicatas de PCR/Ópticas.
4.  **Base Quality Score Recalibration (BQSR - GATK):** * `BaseRecalibrator`: Calcula erros sistemáticos nos scores de qualidade.
    * `ApplyBQSR`: Aplica a correção no BAM.

## 🛠️ Pré-requisitos

* **Ambiente:** Cluster HPC com gerenciador de tarefas **SLURM**.
* **Container:** Imagem Singularity com GATK 4+ e BWA instalados (arquivo `.sif`).
* **Referências:** * Genoma de Referência (FASTA + índices `.fai`, `.dict`, `.bwt`, etc.).
    * Banco de dados de variantes conhecidas (dbSNP, Mills/Indels) para o BQSR.

## 📂 Estrutura de Arquivos

```text
preprocessing/
├── config.yaml                # Arquivo de configuração principal
├── run_analysis.sh            # Script "Gerente" (Submissão SLURM)
├── scripts/
│   ├── create_tasklist.py     # Gera lista de amostras (R1/R2)
│   └── wrapper_preproc.sh     # Script "Operário" (Executa BWA/GATK)
├── logs/                      # Logs de execução (Output/Error)
└── README.md                  # Este arquivo

## Observação:

Os arquivos de referência a serem utilizados podem ser obtidos através do GATK Resource Bundle.
