### ⌛ - Construindo a Pipeline:

Para utilizar este repositório, recomenda-se o uso do arquivo de definição (`.def`) para a construção de um container **Apptainer/Singularity** otimizado para análise de variantes somáticas.

O ambiente foi desenhado para seguir as **GATK Best Practices**, garantindo reprodutibilidade e portabilidade entre estações de trabalho locais e ambientes de computação de alto desempenho (HPC).

Este ambiente será utilizado principalmente para as etapas de `pre-processing` e `tumor-only-calling`, mas também ressalta-se que, caso seja para uso de módulos específicos, também disponibilizam-se YMLs para configurar ambientes conda.

#### Sobre a receita `(.def)`:

Todas as ferramentas foram compiladas a partir do código-fonte ou obtidas de imagens oficiais para garantir estabilidade:

| Ferramenta   | Versão    | Origem/Notas                           |
|:-------------|:----------|:---------------------------------------|
| **GATK**     | 4.6.1.0   | Base Image (Broad Institute)           |
| **BWA**      | 0.7.17    | Compilado (com *fix* para GCC moderno) |
| **Samtools** | 1.21      | Compilado (HTSlib 1.21)                |
| **Bcftools** | 1.21      | Compilado (HTSlib 1.21)                |
| **Picard**   | (Incluso) | Integrado ao GATK 4                    |

#### Pré-requisitos

Para construir a imagem, você precisa ter instalado em sua máquina local (Linux/WSL): **Apptainer** (ou SingularityCE) e acesso `root` (sudo) para a construção.

> **Nota:** Não é necessário ter acesso root para *executar* o container no HPC.

### Como Construir (Build)

Recomenda-se construir a imagem localmente para evitar problemas de permissão em clusters compartilhados.

1.  Clone este repositório:

    ``` bash
    git clone https://github.com/carolinaffernandes/somatic-call.git 
    cd somatic-call 
    ```

2.  Execute o build da imagem .sif:

``` bash
sudo singularity build gatk_4.6.1.sif gatk_versionado.def 
```

> Transfira o arquivo `gatk_4.6.1.sif` gerado para o seu cluster/HPC via SCP ou RSYNC.

3.  Como Usar:

O container permite executar todas as ferramentas encapsuladas sem necessidade de instalação ou configuração de *paths*.

Exemplo: Alinhamento com BWA

``` bash
singularity exec gatk_4.6.1.sif bwa mem \
    -t 16 \
    -R "@RG\tID:amostra\tSM:amostra\tPL:ILLUMINA" \
    ref.fa \
    read1.fq.gz read2.fq.gz > alinhamento.sam
```

Exemplo: Chamada de Variantes (Mutect2)

``` bash
singularity exec gatk_4.6.1.sif gatk Mutect2 \
    -R ref.fa \
    -I tumor.bam \
    -I normal.bam \
    -O somatic.vcf.gz
```

> Nos scripts disponíveis neste repositório, basta adicionar o *path* da sua instalação no 'singularity_img' do config.yaml

Este container foi desenvolvido aderindo aos princípios *FAIR,* servindo como um objeto digital persistente que encapsula o ambiente computacional exato utilizado nas análises.

-   **Arquivo de Definição:** `gatk_versionado.def`

-   **Imagem Imutável:** `gatk_4.6.1.sif`

Autoria: Maria Carolina F. C. Santos
