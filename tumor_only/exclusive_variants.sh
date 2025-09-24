#!/bin/bash
#SBATCH --job-name=isec_filter
#SBATCH --output=isec_filter_%j.out
#SBATCH --error=isec_filter_%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

set -euo pipefail

###### FUNÇÃO DE USO

usage() {
    echo "Uso: $0 -o <output_dir> [-i <input_dir> | VCF1 VCF2 ...] [-d <min_depth>] [-t <variant_type>]"
    echo
    echo "Parâmetros obrigatórios:"
    echo "  -o  Diretório de saída (será criado se não existir)"
    echo
    echo "Parâmetros opcionais:"
    echo "  -i  Diretório contendo arquivos VCF (.vcf.gz)."
    echo "      Se usado, o script processará automaticamente todos os VCFs encontrados no diretório."
    echo "  -d  Cobertura mínima (DP). Default: 10"
    echo "  -t  Tipo de variante: snps, indels ou both. Default: both"
    echo
    echo "Exemplo usando diretório:"
    echo "  sbatch $0 -o results -i /path/to/vcfs -d 20 -t snps"
    echo
    echo "Exemplo passando arquivos:"
    echo "  sbatch $0 -o results -d 20 -t snps sampleA.vcf.gz sampleB.vcf.gz sampleC.vcf.gz"
    exit 1
}

###### PARÂMETROS INICIAIS

INPUT_DIR=""
MIN_DEPTH=10
VARIANT_TYPE="both"
OUTPUT_DIR=""

while getopts "o:d:t:i:" opt; do
    case $opt in
        o) OUTPUT_DIR="$OPTARG" ;;
        d) MIN_DEPTH="$OPTARG" ;;
        t) VARIANT_TYPE="$OPTARG" ;;
        i) INPUT_DIR="$OPTARG" ;;
        *) usage ;;
    esac
done
shift $((OPTIND -1))

###### DEFININDO LISTA DE VCFs

if [[ -n "$INPUT_DIR" ]]; then
    if [[ ! -d "$INPUT_DIR" ]]; then
        echo "ERRO: Diretório $INPUT_DIR não existe."
        exit 1
    fi
    mapfile -t VCF_LIST < <(find "$INPUT_DIR" -maxdepth 1 -type f -name "*.vcf.gz" | sort)
else
    VCF_LIST=("$@")
fi

if [[ ${#VCF_LIST[@]} -lt 2 ]]; then
    echo "ERRO: É necessário fornecer pelo menos dois arquivos VCF."
    usage
fi

mkdir -p "$OUTPUT_DIR"

###### EXIBINDO PARÂMETROS
echo "=== PARÂMETROS DE EXECUÇÃO ==="
echo "Diretório de saída: $OUTPUT_DIR"
[[ -n "$INPUT_DIR" ]] && echo "Diretório de entrada: $INPUT_DIR"
echo "Cobertura mínima (DP): $MIN_DEPTH"
echo "Tipo de variante: $VARIANT_TYPE"
echo "VCFs de entrada: ${VCF_LIST[*]}"
echo "=============================="

###### PASSO 1 - Determinar parâmetro de tipo
TYPE_PARAM=""
if [[ "$VARIANT_TYPE" == "snps" ]]; then
    TYPE_PARAM="-v snps"
elif [[ "$VARIANT_TYPE" == "indels" ]]; then
    TYPE_PARAM="-v indels"
fi

###### PASSO 2 - Filtrar VCFs individuais
echo "[PASSO 2] Filtrando cada VCF (PASS + tipo + FORMAT/DP >= $MIN_DEPTH)..."

FILTERED_VCFS=()

for VCF in "${VCF_LIST[@]}"; do
    BASENAME=$(basename "$VCF" .vcf.gz)
    FILTERED_VCF="$OUTPUT_DIR/${BASENAME}_filtered.vcf.gz"

    echo "   -> Processando $VCF"
    bcftools view -f PASS $TYPE_PARAM "$VCF" -Ou | \
    bcftools filter -i "FORMAT/DP>=$MIN_DEPTH" -Oz -o "$FILTERED_VCF"

    tabix -p vcf "$FILTERED_VCF"
    FILTERED_VCFS+=("$FILTERED_VCF")
done

###### PASSO 3 - Rodar bcftools isec para variantes únicas e comuns
echo "[PASSO 3] Executando bcftools isec..."
bcftools isec -p "$OUTPUT_DIR/isec_output" -Oz "${FILTERED_VCFS[@]}"

###### PASSO 4 - Renomear arquivos de saída (únicos e comuns)
echo "[PASSO 4] Renomeando arquivos de saída..."

INDEX=0
for VCF in "${FILTERED_VCFS[@]}"; do
    SAMPLE_NAME=$(basename "$VCF" _filtered.vcf.gz)
    mv "$OUTPUT_DIR/isec_output/000${INDEX}.vcf.gz" "$OUTPUT_DIR/${SAMPLE_NAME}_unique.vcf.gz"
    tabix -p vcf "$OUTPUT_DIR/${SAMPLE_NAME}_unique.vcf.gz"
    INDEX=$((INDEX + 1))
done

# Arquivo com variantes comuns a todos
if [[ -f "$OUTPUT_DIR/isec_output/000${INDEX}.vcf.gz" ]]; then
    mv "$OUTPUT_DIR/isec_output/000${INDEX}.vcf.gz" "$OUTPUT_DIR/common_variants.vcf.gz"
    tabix -p vcf "$OUTPUT_DIR/common_variants.vcf.gz"
fi

###### PASSO 5 - Gerar a união real de todos os VCFs filtrados
echo "[PASSO 5] Criando combined_all.vcf.gz (união de todas as variantes filtradas)..."


bcftools merge -Oz -o "$OUTPUT_DIR/combined_all.vcf.gz" "${FILTERED_VCFS[@]}"
tabix -p vcf "$OUTPUT_DIR/combined_all.vcf.gz"

###### PASSO 6 - Anotar com ANNOVAR
echo "[PASSO 6] Anotando arquivos com ANNOVAR..."

ANNOVAR_BIN="/scratch/14196861/annovar"
ANNOVAR_DIR="$OUTPUT_DIR/Anotados"
mkdir -p "$ANNOVAR_DIR"

# Lista de arquivos para anotar
ANNOVAR_FILES=("${FILTERED_VCFS[@]}" "$OUTPUT_DIR"/*_unique.vcf.gz "$OUTPUT_DIR/common_variants.vcf.gz" "$OUTPUT_DIR/combined_all.vcf.gz")

for VCF in "${ANNOVAR_FILES[@]}"; do
    BASENAME=$(basename "$VCF" .vcf.gz)
    ANNOT_TYPE="$VARIANT_TYPE"
    ANNOVAR_PREFIX="${BASENAME}_${ANNOT_TYPE}"

    echo "   -> Executando table_annovar.pl para $VCF..."
    "$ANNOVAR_BIN/table_annovar.pl" "$VCF" "$ANNOVAR_BIN/humandb/" \
        -buildver hg38 \
        -out "$ANNOVAR_DIR/${ANNOVAR_PREFIX}" \
        -remove \
        -protocol refGene,avsnp150 \
        -operation g,f \
        -nastring . \
        -vcfinput
done

echo "[FINALIZADO] ANNOVAR"
echo "Arquivos gerados em: $ANNOVAR_DIR"

