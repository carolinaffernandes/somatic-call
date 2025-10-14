------------------------------------------------------------------------

# 🧬 Pipeline de Chamada e Filtragem de Variantes

(Mutect2 + FilterMutectCalls + MergeVCFs)

Este repositório contém um pipeline modular em **Bash** para realizar chamadas de variantes somáticas utilizando o **GATK Mutect2**, seguido de filtragem com **FilterMutectCalls** e junção dos resultados com **Picard MergeVCFs**.\
O pipeline foi projetado para ambientes **HPC com SLURM** e execução via **Singularity** (ou Docker, se adaptado).

------------------------------------------------------------------------

## Estrutura do Projeto

```         
project/
.
├── chrom_list.txt
├── config.yaml
├── environment_variant-calling.yml
├── generate_task_list.sh
├── get_yaml.py
├── README.md
├── results
├── run_pipeline.sh
├── samples
├── samplesheet.tsv
└── scripts
    ├── download_ref.sh
    ├── filtermutectcalls_chr.sh
    ├── merge_vcf.sh
    └── mutect2_chr.sh

```

------------------------------------------------------------------------

## 1. Configuração

### Controle do ambiente com conda

```         
conda env create -f environment_variant-calling.yml
```

### Conferência do `config.yaml`

Defina os caminhos e parâmetros principais do pipeline

### Arquivo `samplesheet.tsv`

O arquivo deve conter **duas colunas separadas por tabulação**:

```         
sample_id   tumor_bam
Sample01    /path/tumor.bam
Sample02    /path/tumor.bam
```

### Download dos arquivos de referência

> Para esta etapa, você precisa ter o `gsutil` configurado (via `gcloud auth login` ou `gsutil config`).

Os arquivos para o genoma de referência serão transferidos para o diretório ref. O download será feito do Resource Bundle do GATK. Acessar para maiores informações.

Executar:

```         
scripts/download_ref.sh
```

------------------------------------------------------------------------

## 🚀 2. Execução do Pipeline

### Submissão via SLURM

O pipeline foi projetado para rodar em um **HPC com SLURM**.\
Para iniciar, execute:

``` bash
sbatch run_pipeline.sh
```

O script `run_pipeline.sh`:

\- Lê `samplesheet.tsv`

\- Envia jobs SLURM para cada amostra

\- Dentro de cada job, roda `mutect2_chr.sh` por cromossomo (em paralelo)

\- Em seguida, chama `filtermutectcalls_chr.sh`

\- Por fim, o `merge_vcfs.sh` junta todos os resultados em um único arquivo por amostra.

### Controle de Dependências no SLURM

Os scripts utilizam dependências (`--dependency=afterok:<JOBID>`) para garantir que: - `FilterMutectCalls` só inicie após todas as etapas do `Mutect2` estarem completas; - `MergeVCFs` só rode após a filtragem.

------------------------------------------------------------------------

## 3. Modularidade

Cada etapa do pipeline pode ser executada separadamente:

``` bash
bash download_ref.sh	# Download de arquivos de referência
bash mutect2_chr.sh   # Chamada de variantes
bash filtermutectcalls_chr.sh  # Filtragem
bash merge_vcfs.sh  # Junção final
```

Isso facilita o reprocessamento parcial sem repetir todo o fluxo.

------------------------------------------------------------------------

## 4. Observações Importantes

-   Todos os diretórios de trabalho (`output_dir`, `temp_dir`) devem existir e ser acessíveis.\
-   O pipeline não utiliza `--bind` no Singularity, pois assume que `/scratch` e `/temporario2` já estão montados.\
-   As imagens Docker devem estar disponíveis localmente ou acessíveis via cache.

------------------------------------------------------------------------

## 5. Resultados Esperados

Ao final da execução, o diretório de saída (`results/`) conterá:

```         
results/
├── Sample01/
│   ├── mutect2/
│   │   ├── Sample01_chr1.vcf.gz
│   │   ├── ...
│   ├── filtered/
│   │   ├── Sample01_chr1_filtered.vcf.gz
│   │   ├── ...
│   └── merged/
│       └── Sample01_merged.vcf.gz
```

