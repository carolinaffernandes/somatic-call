#!/usr/bin/env python3
import os
import sys
import yaml
import pandas as pd
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup

# --- CONFIGURAÇÃO ---
def load_config(path="config.yaml"):
    if not os.path.exists(path): return {}
    with open(path, "r") as f: data = yaml.safe_load(f)
    return data.get("annotation_pipeline", data)

# --- FUNÇÕES AUXILIARES ---
def check_cosmic(value):
    val = str(value).strip()
    if val == '.' or val == '' or pd.isna(value): return "Not in COSMIC"
    return "In COSMIC"

def get_cosmic_details(df):
    if 'cosmic70' not in df.columns: return []
    hits = df[df['cosmic70'].astype(str) != '.'].copy()
    details = []
    for _, row in hits.iterrows():
        change = row.get('AAChange.refGene', '.')
        if change == '.': change = "N/A"
        details.append({
            'gene': row.get('Gene.refGene', 'Unknown'),
            'cosmic_id': row.get('cosmic70', 'Unknown'),
            'change': change,
            'func': row.get('ExonicFunc.refGene', row.get('Func.refGene', ''))
        })
    return details

def get_gene_counts(df):
    relevant_funcs = ['exonic', 'splicing', 'exonic;splicing', 'UTR3', 'UTR5']
    if 'Func.refGene' in df.columns:
        df_genes = df[df['Func.refGene'].isin(relevant_funcs)].copy()
    else: return {}

    all_genes = []
    for entry in df_genes['Gene.refGene'].dropna():
        genes = str(entry).split(';')
        all_genes.extend(genes)
    
    return pd.Series(all_genes).value_counts().head(20).to_dict()

# --- PARSE ANNOVAR ---
def parse_annovar(file_path):
    try:
        sep = ',' if file_path.endswith('.csv') else '\t'
        header = pd.read_csv(file_path, sep=sep, nrows=0).columns.tolist()
        
        cols_needed = ['Func.refGene', 'ExonicFunc.refGene', 'cosmic70', 'Gene.refGene', 'AAChange.refGene']
        use_cols = [c for c in cols_needed if c in header]
        
        if not use_cols:
            print(f"[AVISO] Colunas necessárias ausentes em {file_path}")
            return None

        df = pd.read_csv(file_path, sep=sep, usecols=use_cols, low_memory=False)
        
        # 1. Estatísticas
        counts_func = {}
        if 'Func.refGene' in df.columns:
            counts_func = df['Func.refGene'].value_counts().head(10).to_dict()

        counts_exonic = {}
        if 'ExonicFunc.refGene' in df.columns:
            counts_exonic = df[df['ExonicFunc.refGene'] != '.']['ExonicFunc.refGene'].value_counts().head(10).to_dict()

        # 2. COSMIC
        hits_cosmic = 0
        cosmic_details = []
        if 'cosmic70' in df.columns:
            df['COSMIC_Status'] = df['cosmic70'].apply(check_cosmic)
            hits_cosmic = df['COSMIC_Status'].value_counts().get("In COSMIC", 0)
            cosmic_details = get_cosmic_details(df)

        # 3. Genes
        gene_counts = {}
        if 'Gene.refGene' in df.columns:
            gene_counts = get_gene_counts(df)
        
        print(f"   -> {len(df)} vars | {hits_cosmic} COSMIC Hits")

        return {
            "total": len(df),
            "counts_func": counts_func,
            "counts_exonic": counts_exonic,
            "hits_cosmic": hits_cosmic,
            "cosmic_details": cosmic_details,
            "gene_counts": gene_counts
        }

    except Exception as e:
        print(f"[ERRO] {file_path}: {e}")
        return None

# --- OUTPUT FUNCTIONS ---

def save_tsv(data_dict, output_path, col_names=["Category", "Count"]):
    """Salva dicionário como TSV simples."""
    if not data_dict: return
    try:
        df = pd.DataFrame(list(data_dict.items()), columns=col_names)
        df.to_csv(output_path, sep='\t', index=False)
    except Exception as e:
        print(f"[ERRO] Falha ao salvar TSV {output_path}: {e}")

def save_bar_plot(data_dict, title, outpath, color='skyblue'):
    if not data_dict: return
    # Ajustei figsize para (8, 10) para ficar mais alto/estreito
    # já que agora serão 3 colunas, precisamos de altura para ler os labels
    plt.figure(figsize=(8, 10)) 
    sorted_items = sorted(data_dict.items(), key=lambda x: x[1], reverse=True)
    labels = [str(k) for k, v in sorted_items]
    values = [v for k, v in sorted_items]

    plt.barh(labels, values, color=color, edgecolor='black', alpha=0.8)
    plt.title(title, fontsize=16, fontweight='bold')
    plt.xlabel("Contagem", fontsize=14)
    plt.yticks(fontsize=12)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(outpath, dpi=100)
    plt.close()

# --- RELATÓRIO HTML ---
def generate_html_report(output_dir, results_list, report_title, filename):
    html_path = os.path.join(output_dir, filename)
    
    style = """
    <style>
        body { font-family: 'Segoe UI', sans-serif; padding: 20px; background-color: #f4f6f9; color: #333; }
        h1 { color: #2c3e50; border-bottom: 4px solid #3498db; padding-bottom: 15px; }
        h2 { color: #2c3e50; margin-top: 40px; border-left: 6px solid #3498db; padding-left: 15px; background: #eaf2f8; padding-top:10px; padding-bottom:10px;}
        
        /* Stats Cards */
        .stats-grid { display: flex; gap: 20px; margin-bottom: 20px; }
        .stat-card { background: white; padding: 15px; border-radius: 8px; flex: 1; text-align: center; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
        .stat-value { font-size: 28px; font-weight: bold; color: #2980b9; }
        .stat-label { font-size: 13px; color: #7f8c8d; text-transform: uppercase; }
        
        /* Grid de Plots: 3 colunas IGUAIS */
        .plot-container { 
            display: grid; 
            grid-template-columns: repeat(3, 1fr); /* 3 colunas iguais */
            gap: 15px; 
            width: 100%;
            margin-bottom: 30px;
        }
        .plot-box { 
            background: white; 
            padding: 10px; 
            border-radius: 8px; 
            box-shadow: 0 2px 5px rgba(0,0,0,0.05); 
            text-align: center;
        }
        .plot-img { width: 100%; height: auto; border: 1px solid #eee; }

        /* Tabela COSMIC */
        .cosmic-section { margin-top: 30px; background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 5px rgba(0,0,0,0.05); border-left: 5px solid #c0392b; }
        .cosmic-title { font-size: 18px; font-weight: bold; color: #c0392b; margin-bottom: 15px; }
        table { width: 100%; border-collapse: collapse; margin-top: 10px; font-size: 14px; }
        th, td { padding: 10px; text-align: left; border-bottom: 1px solid #ddd; }
        th { background-color: #f8f9fa; color: #555; }
        tr:hover { background-color: #f1f1f1; }
    </style>
    """

    soup = BeautifulSoup(f"<html><head>{style}</head><body><h1>{report_title}</h1></body></html>", "html.parser")

    for res in results_list:
        sample = res['sample']
        div_sample = soup.new_tag("div", attrs={"class": "sample-wrapper"})
        
        h2 = soup.new_tag("h2"); h2.string = f"Amostra: {sample}"; div_sample.append(h2)

        # Cards
        div_stats = soup.new_tag("div", attrs={"class": "stats-grid"})
        stats_data = [
            ("Total Variantes", f"{res['total']:,}"),
            ("Exonic Variants", sum(res['counts_exonic'].values())),
            ("COSMIC Hits", res['hits_cosmic'])
        ]
        for label, val in stats_data:
            card = soup.new_tag("div", attrs={"class": "stat-card"})
            v = soup.new_tag("div", attrs={"class": "stat-value"}); v.string = str(val)
            l = soup.new_tag("div", attrs={"class": "stat-label"}); l.string = label
            card.append(v); card.append(l)
            div_stats.append(card)
        div_sample.append(div_stats)

        # Plots (3 na mesma linha)
        div_plots = soup.new_tag("div", attrs={"class": "plot-container"})
        
        plot_files = [
            ("Func", "Distribuição Genômica"),
            ("Exonic", "Classificação Exônica"),
            ("Genes", "Top Genes Mutados")
        ]
        
        for p_type, p_desc in plot_files:
            img_name = f"{sample}_{p_type}.png"
            if os.path.exists(os.path.join(output_dir, "plots", img_name)):
                box = soup.new_tag("div", attrs={"class": "plot-box"})
                img = soup.new_tag("img", src=f"plots/{img_name}", attrs={"class": "plot-img"})
                lbl = soup.new_tag("div", attrs={"style": "font-weight:bold; margin-top:5px; color:#555;"})
                lbl.string = p_desc
                box.append(img); box.append(lbl)
                div_plots.append(box)
        
        div_sample.append(div_plots)

        # Tabela COSMIC (Única tabela mantida no HTML)
        if res['hits_cosmic'] > 0:
            div_cosmic = soup.new_tag("div", attrs={"class": "cosmic-section"})
            ctitle = soup.new_tag("div", attrs={"class": "cosmic-title"})
            ctitle.string = f"Achados Clínicos (COSMIC) - {res['hits_cosmic']} variantes"
            div_cosmic.append(ctitle)

            table = soup.new_tag("table")
            thead = soup.new_tag("thead"); tr_head = soup.new_tag("tr")
            for h in ["Gene", "COSMIC ID", "Alteração (AA)", "Tipo"]:
                th = soup.new_tag("th"); th.string = h; tr_head.append(th)
            thead.append(tr_head); table.append(thead)

            tbody = soup.new_tag("tbody")
            for det in res['cosmic_details']:
                tr = soup.new_tag("tr")
                for k in ['gene', 'cosmic_id', 'change', 'func']:
                    td = soup.new_tag("td"); td.string = str(det[k]); tr.append(td)
                tbody.append(tr)
            
            table.append(tbody)
            div_cosmic.append(table)
            div_sample.append(div_cosmic)

        soup.body.append(div_sample)

    with open(html_path, "w") as f: f.write(str(soup))
    print(f"HTML criado: {html_path}")

# --- MAIN ---
def main():
    cfg = load_config()
    html_filename = cfg.get("html_filename", "Relatorio_Clean.html")
    report_title = cfg.get("report_title", "Relatório de Variantes Somáticas")
    output_dir = cfg.get("output_dir", "results_report_clean")
    samplesheet = cfg.get("samplesheet", "samples.tsv")
    
    if not os.path.exists(samplesheet):
        print("[ERRO] Samplesheet 'samples.tsv' não encontrada."); return

    try:
        sep = '\t' if samplesheet.endswith('.tsv') else ','
        df = pd.read_csv(samplesheet, sep=sep)
    except: print("[ERRO] Erro ao ler samplesheet."); return

    # Cria pastas de saída
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, "plots"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "tables"), exist_ok=True) # Pasta nova para TSVs

    results = []
    for _, row in df.iterrows():
        sample = str(row['sample'])
        fpath = str(row['file'])
        
        print(f"Processando {sample}...")
        data = parse_annovar(fpath)
        if data:
            data['sample'] = sample
            results.append(data)
            
            # 1. Gera Plots
            save_bar_plot(data['counts_func'], f"Regiões - {sample}", 
                          os.path.join(output_dir, "plots", f"{sample}_Func.png"), "#27ae60")
            save_bar_plot(data['counts_exonic'], f"Tipos Exônicos - {sample}", 
                          os.path.join(output_dir, "plots", f"{sample}_Exonic.png"), "#2980b9")
            save_bar_plot(data['gene_counts'], f"Top Genes - {sample}", 
                          os.path.join(output_dir, "plots", f"{sample}_Genes.png"), "#8e44ad")

            # 2. Gera TSVs Separados
            save_tsv(data['counts_func'], os.path.join(output_dir, "tables", f"{sample}_regions.tsv"), ["Region", "Count"])
            save_tsv(data['counts_exonic'], os.path.join(output_dir, "tables", f"{sample}_exonic_types.tsv"), ["Type", "Count"])
            save_tsv(data['gene_counts'], os.path.join(output_dir, "tables", f"{sample}_top_genes.tsv"), ["Gene", "Count"])
            if data['cosmic_details']:
                # Salva cosmic também em TSV para garantir
                pd.DataFrame(data['cosmic_details']).to_csv(os.path.join(output_dir, "tables", f"{sample}_cosmic.tsv"), sep='\t', index=False)

    if results:
        generate_html_report(output_dir, results, report_title, html_filename)
        print("[OK] Relatório e Tabelas TSV gerados com sucesso!")

main()
