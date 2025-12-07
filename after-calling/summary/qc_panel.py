#!/usr/bin/env python3
import os
import sys
import yaml
import gzip
import pandas as pd
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup

def load_config(path="config.yaml"):
    if not os.path.exists(path):
        print("[ERRO] Arquivo config.yaml não encontrado.")
        sys.exit(1)
    with open(path, "r") as f:
        return yaml.safe_load(f)["qc_pipeline"]

def parse_vcf(vcf_path, dp_field="DP", af_field="AF"):
    snps = 0
    indels = 0
    transitions = 0
    transversions = 0
    dp_values = []
    af_values = []
    ti_set = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}

    # Lógica para abrir arquivo GZ ou Texto
    if vcf_path.endswith(".gz"):
        open_func = gzip.open
        mode = "rt"
    else:
        open_func = open
        mode = "r"

    try:
        with open_func(vcf_path, mode) as vcf:
            for line in vcf:
                if line.startswith("#"): continue
                cols = line.strip().split("\t")
                if len(cols) < 10: continue

                ref = cols[3].upper()
                alt = cols[4].upper()

                if len(ref) == 1 and len(alt) == 1:
                    snps += 1
                    if (ref, alt) in ti_set: transitions += 1
                    else: transversions += 1
                else:
                    indels += 1

                format_fields = cols[8].split(":")
                sample_fields = cols[9].split(":")
                fmt = dict(zip(format_fields, sample_fields))

                if dp_field in fmt and fmt[dp_field].isdigit():
                    dp_values.append(int(fmt[dp_field]))

                if af_field in fmt:
                    try:
                        val = fmt[af_field].split(",")[0] 
                        af_values.append(float(val))
                    except: pass
    except Exception as e:
        print(f"[ERRO] Falha ao ler {vcf_path}: {e}")
        return None

    titv_ratio = transitions / transversions if transversions > 0 else 0.0

    return {
        "snps": snps, "indels": indels, "titv_ratio": round(titv_ratio, 2),
        "dp_values": dp_values, "af_values": af_values
    }

def save_plot(values, title, outpath, xmax=None, color='skyblue'):
    if not values: return
    plt.figure(figsize=(6,4))
    histo_range = (0, xmax) if xmax else None
    plt.hist(values, bins=40, range=histo_range, color=color, edgecolor='black', alpha=0.8)
    plt.title(title)
    plt.xlabel("Valor")
    plt.ylabel("Frequência")
    plt.grid(axis='y', alpha=0.3)
    if xmax: plt.xlim(0, xmax)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()

def generate_html_report(output_dir, df, report_title="Painel de QC", filename="painel_QC.html"):
    html_path = os.path.join(output_dir, filename)
    
    style = """
    <style>
        body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; padding: 20px; background-color: #f9f9f9; }
        h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
        h2 { color: #34495e; margin-top: 30px; }
        table { border-collapse: collapse; width: 100%; background: white; box-shadow: 0 1px 3px rgba(0,0,0,0.2); }
        th, td { text-align: left; padding: 12px; border-bottom: 1px solid #ddd; }
        th { background-color: #3498db; color: white; }
        tr:hover { background-color: #f1f1f1; }
        .plot-container { display: flex; flex-wrap: wrap; gap: 20px; justify-content: center; }
        .plot-box { background: white; border: 1px solid #ddd; padding: 15px; border-radius: 8px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); text-align: center; }
        .plot-title { font-weight: bold; margin-bottom: 10px; color: #555; }
    </style>
    """

    soup = BeautifulSoup(f"<html><head>{style}</head><body><h1>{report_title}</h1></body></html>", "html.parser")

    # --- FORMATADORES PARA O HTML ---
    formatters = {
        'titv_ratio': '{:,.2f}'.format, 
        'mean_DP': '{:,.1f}'.format, 
        'median_DP': '{:,.1f}'.format,  # Adicionado Mediana DP
        'mean_AF': '{:,.2f}'.format,
        'median_AF': '{:,.2f}'.format   # Adicionado Mediana AF
    }
    
    # --- COLUNAS PARA EXIBIR ---
    cols_to_show = [
        'sample', 'snps', 'indels', 'titv_ratio', 
        'mean_DP', 'median_DP',  # Adicionado aqui
        'mean_AF', 'median_AF'   # Adicionado aqui
    ]
    cols_final = [c for c in cols_to_show if c in df.columns]
    
    table_html = df[cols_final].to_html(index=False, formatters=formatters, classes="sortable")
    
    soup.body.append(BeautifulSoup("<h2>Métricas por Amostra</h2>", "html.parser"))
    soup.body.append(BeautifulSoup(table_html, "html.parser"))

    soup.body.append(BeautifulSoup("<h2>Distribuições (DP e AF)</h2>", "html.parser"))
    plots_div = soup.new_tag("div", attrs={"class": "plot-container"})
    plots_dir = os.path.join(output_dir, "plots")
    
    for sample in df['sample']:
        sample_plots = sorted([p for p in os.listdir(plots_dir) if p.startswith(f"{sample}_")])
        for img_file in sample_plots:
            div_box = soup.new_tag("div", attrs={"class": "plot-box"})
            title_div = soup.new_tag("div", attrs={"class": "plot-title"})
            display_name = img_file.replace(".png", "").replace(f"{sample}_", "")
            title_div.string = f"{sample} - {display_name}"
            div_box.append(title_div)
            img_tag = soup.new_tag("img", src=f"plots/{img_file}", width="400")
            div_box.append(img_tag)
            plots_div.append(div_box)

    soup.body.append(plots_div)
    with open(html_path, "w") as f: f.write(str(soup))
    print(f"HTML final criado em: {html_path}")

def main():
    cfg = load_config()
    output_dir = cfg.get("output_dir", "qc_output")
    samplesheet_path = cfg.get("samplesheet", "samples.tsv")
    
    report_title = cfg.get("report_title", "Painel de QC")
    html_filename = cfg.get("html_filename", "painel_QC.html")
    dp_color = cfg.get("dp_color", "skyblue")
    af_color = cfg.get("af_color", "orange")
    
    dp_field = cfg.get("dp_field", "DP")
    af_field = cfg.get("af_field", "AF")

    if not os.path.exists(samplesheet_path):
        print(f"[ERRO] Samplesheet não encontrada: {samplesheet_path}")
        return

    print(f"Lendo samplesheet: {samplesheet_path}")
    try:
        if samplesheet_path.endswith('.tsv') or samplesheet_path.endswith('.txt'):
            df_samples = pd.read_csv(samplesheet_path, sep='\t')
        else:
            df_samples = pd.read_csv(samplesheet_path)
    except Exception as e:
        print(f"[ERRO] Falha ao ler a planilha: {e}"); return
    
    if "sample" not in df_samples.columns or "vcf" not in df_samples.columns:
        print("[ERRO] Samplesheet deve ter colunas 'sample' e 'vcf'"); return

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, "plots"), exist_ok=True)

    results = []
    for _, row in df_samples.iterrows():
        sample = str(row['sample'])
        vcf_path = str(row['vcf'])
        print(f"Processando {sample}...")
        
        metrics = parse_vcf(vcf_path, dp_field, af_field)
        if metrics is None: continue
        
        dp = metrics.pop("dp_values")
        af = metrics.pop("af_values")

        plot_dp = os.path.join(output_dir, "plots", f"{sample}_DP.png")
        plot_af = os.path.join(output_dir, "plots", f"{sample}_AF.png")

        if dp: save_plot(dp, f"DP — {sample}", plot_dp, xmax=100, color=dp_color)
        if af: save_plot(af, f"AF — {sample}", plot_af, color=af_color)

        metrics["sample"] = sample
        
        # --- CÁLCULO DE MÉDIA E MEDIANA ---
        if dp:
            series_dp = pd.Series(dp)
            metrics["mean_DP"] = series_dp.mean()
            metrics["median_DP"] = series_dp.median() # De volta!
        else:
            metrics["mean_DP"] = None
            metrics["median_DP"] = None

        if af:
            series_af = pd.Series(af)
            metrics["mean_AF"] = series_af.mean()
            metrics["median_AF"] = series_af.median() # De volta!
        else:
            metrics["mean_AF"] = None
            metrics["median_AF"] = None
            
        results.append(metrics)

    df_final = pd.DataFrame(results)
    if not df_final.empty:
        df_final.sort_values(by="sample", inplace=True)

    csv_path = os.path.join(output_dir, "qc_panel.csv")
    df_final.to_csv(csv_path, index=False)
    print(f"[OK] CSV salvo em {csv_path}")

    if not df_final.empty:
        generate_html_report(output_dir, df_final, report_title, html_filename)

main()
