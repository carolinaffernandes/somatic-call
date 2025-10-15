------------------------------------------------------------------------

# 🧬 Filtragem e Anotação de Variantes exclusivas (hard-filter)

(Bcftools + ANNOVAR)

Este repositório contém um pipeline modular em **Bash** para filtragem e anotação de variantes somáticas, exclusivas para cada amostra, utilizando o bcftools e ANNOVAR. O pipeline foi projetado para ambientes **HPC com SLURM** e execução via **Singularity** (ou Docker, se adaptado).

> É esperado que sejam realizados os processos anteriores (../tumor_only_calling/) e que os resultados (tumor_only_calling/results) obtidos, sejam utilizados como samples nesta etapa.

------------------------------------------------------------------------

## 1. Configuração

### Controle do ambiente com conda

```         
conda env create -f environment_variant-calling.yml
```

### Arquivo `samplesheet.tsv`

O arquivo deve conter **duas colunas separadas por tabulação**:

-   As amostras podem ser isoladas no diretório (samples/), mas também podem ser referenciadas em seu path absoluto no diretório (results/) do tumor_only_calling.

-   Para a anotação de variantes e conversão de .vcf.gz em .vcf e .txt "anotados", utilizaremos o ANNOVAR que deve ser transferido mediante registro conforme instruído no [site oficial do ANNOVAR](https://annovar.openbioinformatics.org/en/latest/user-guide/download)
