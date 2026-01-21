import os
import glob
import yaml

# Define o caminho absoluto do arquivo de saida
TASK_LIST_PATH = "/scratch/14196861/doutorado-U251MG/somatic-call/preprocessing/preproc_task_list.txt"
CONFIG_PATH = "../config.yaml"

def load_config(path):
    abs_config_path = os.path.abspath(path)
    if not os.path.exists(abs_config_path):
        print(f"ERRO: Config nao encontrado em: {abs_config_path}")
        exit(1)
    with open(abs_config_path, 'r') as file:
        return yaml.safe_load(file)

def create_task_list():
    print("--- Gerador de Lista de Tarefas (Preprocessing) ---")
    
    config = load_config(CONFIG_PATH)
    
    # Se nao tiver no config, edite a linha abaixo com o caminho fixo dos FASTQs
    fastq_dir_raw = "/scratch/14196861/doutorado-U251MG/fastq-filtered" 
    
    if not os.path.exists(fastq_dir_raw):
        print(f"ERRO: Pasta de FASTQ nao existe: {fastq_dir_raw}")
        return

    print(f"Buscando FASTQs em: {fastq_dir_raw}")

    # Padrao de busca (ajuste se seus arquivos forem .fq.gz ou tiverem outro nome)
    # Assume padrao _R1.fastq.gz e _R2.fastq.gz
    r1_files = sorted(glob.glob(os.path.join(fastq_dir_raw, "*_R1.fastq.gz")))
    
    if not r1_files:
        print("AVISO: Nenhum arquivo _R1.fastq.gz encontrado.")
        return

    print(f"Encontrados {len(r1_files)} pares de arquivos.")

    with open(TASK_LIST_PATH, 'w') as f:
        for r1 in r1_files:
            # Determina o R2 baseado no R1
            r2 = r1.replace("_R1.fastq.gz", "_R2.fastq.gz")
            
            if not os.path.exists(r2):
                print(f"ERRO: Par R2 nao encontrado para: {r1}")
                continue
                
            # Define ID da amostra (nome do arquivo sem _R1...)
            filename = os.path.basename(r1)
            sample_id = filename.replace("_R1.fastq.gz", "")
            
            # Escreve: SAMPLE <TAB> R1 <TAB> R2
            f.write(f"{sample_id}\t{r1}\t{r2}\n")

    print(f"Lista salva em: {TASK_LIST_PATH}")

if __name__ == "__main__":
    create_task_list()
