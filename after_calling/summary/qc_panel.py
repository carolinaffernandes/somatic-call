#!/usr/bin/env python3
import os
import yaml
import pandas as pd
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup

def load_config(path="config.yaml"):
    if not os.path.exists(path):
        return {
            "qc_pipeline": {
                "input_dir": "vcfs",
                "output_dir": "qc_results",
                "dp_field": "DP",
                "af_field": "AF"
            }
        }
    with open(path, "r") as f:
        return yaml.safe_load(f)["qc_pipeline"]

def parse_vcf(vcf_path, dp_field="DP", af_field="AF"):
    snps = 0
    indels = 0
    transitions = 0
    transversions = 0
    # Removidos contadores de Het/Hom
    
    dp_values = []
    af_values = []

    ti_set = {
        ('A', 'G'), ('G', 'A'),
        ('C', 'T'), ('T', 'C')
    }

    with open(vcf_path, "r") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue

            cols = line.strip().split("\t")
            ref = cols[3].upper()
            alt = cols[4].upper()

            # --- SNP vs INDEL e Ti/Tv ---
            if len(ref) == 1 and len(alt) == 1:
                snps += 1
                if (ref, alt) in ti_set:
                    transitions += 1
                else:
                    transversions += 1
            else:
                indels += 1

            # --- FORMAT parsing ---
            format_fields = cols[8].split(":")
            sample_fields = cols[9].split(":")
            fmt = dict(zip(format_fields, sample_fields))

            # (Removida lógica de GT/Genotype aqui)

            # --- DP (Depth) ---
            if dp_field in fmt and fmt[dp_field].isdigit():
                dp_values.append(int(fmt[dp_field]))

            # --- AF (Allele Frequency) ---
            if af_field in fmt:
                try:
                    val = fmt[af_field].split(",")[0] 
                    af_values.append(float(val))
                except:
                    pass

    titv_ratio = transitions / transversions if transversions > 0 else 0.0

    return {
        "snps": snps,
        "indels": indels,
        "titv_ratio": round(titv_ratio, 2),
        "dp_values": dp_values,
        "af_values": af_values
    }

def save_plot(values, title, outpath, xmax=None):
    if not values:
        return
    
    plt.figure(figsize=(6,4))
    
    histo_range = (0, xmax) if xmax else None
    
    plt.hist(values, bins=40, range=histo_range, color='skyblue', edgecolor='black')
    
    plt.title(title)
    plt.xlabel("Valor")
    plt.ylabel("Frequência")
    
    if xmax:
        plt.xlim(0, xmax)
        
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()

def generate_html_report(output_dir, df):
    html_path = os.path.join(output_dir, "painel_QC.html")
    
    style = """
    <style>
        body { font-family: sans-serif; padding: 20px; }
        table { border-collapse: collapse; width: 100%; margin-bottom: 30px; }
        th, td { text-align: left; padding: 8px; border-bottom: 1px solid #ddd; }
        th { background-color: #f2f2f2; }
        tr:hover { background-color: #f5f5f5; }
        h1, h2 { color: #333; }
        .plot-container { display: flex; flex-wrap: wrap; gap: 20px; }
        .plot-box { border: 1px solid #eee; padding: 10px; border-radius: 5px; }
    </style>
    """

    soup = BeautifulSoup(f"<html><head>{style}</head><body><h1>Painel de QC de Variantes (Somatic/Mutect)</h1></body></html>", "html.parser")

    formatters = {
        'titv_ratio': '{:,.2f}'.format,
        'mean_DP': '{:,.1f}'.format,
        'mean_AF': '{:,.2f}'.format
    }
    table_html = df.to_html(index=False, formatters=formatters, classes="sortable")
    
    soup.body.append(BeautifulSoup("<h2>Métricas por Amostra</h2>", "html.parser"))
    soup.body.append(BeautifulSoup(table_html, "html.parser"))

    soup.body.append(BeautifulSoup("<h2>Distribuições (DP e AF)</h2>", "html.parser"))
    plots_div = soup.new_tag("div", attrs={"class": "plot-container"})
    
    plots_dir = os.path.join(output_dir, "plots")
    
    # Itera sobre o DF ordenado para manter a ordem dos plots
    for sample in df['sample']:
        sample_plots = sorted([p for p in os.listdir(plots_dir) if p.startswith(f"{sample}_")])
        
        for img_file in sample_plots:
            div_box = soup.new_tag("div", attrs={"class": "plot-box"})
            title = soup.new_tag("div")
            title.string = img_file.replace(".png", "")
            div_box.append(title)
            img_tag = soup.new_tag("img", src=f"plots/{img_file}", width="400")
            div_box.append(img_tag)
            plots_div.append(div_box)

    soup.body.append(plots_div)

    with open(html_path, "w") as f:
        f.write(str(soup))

    print(f"HTML final criado em: {html_path}")

def main():
    cfg = load_config()
    input_dir = cfg.get("input_dir", ".")
    output_dir = cfg.get("output_dir", "qc_output")
    dp_field = cfg.get("dp_field", "DP")
    af_field = cfg.get("af_field", "AF")

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, "plots"), exist_ok=True)

    results = []
    files = [f for f in os.listdir(input_dir) if f.endswith(".vcf")]
    
    if not files:
        print(f"Nenhum arquivo .vcf encontrado em {input_dir}")
        return

    for file in files:
        sample = file.replace(".vcf", "")
        vcf_path = os.path.join(input_dir, file)
        print(f"Processando {sample}...")

        metrics = parse_vcf(vcf_path, dp_field, af_field)
        
        dp = metrics.pop("dp_values")
        af = metrics.pop("af_values")

        plot_dp = os.path.join(output_dir, "plots", f"{sample}_DP.png")
        plot_af = os.path.join(output_dir, "plots", f"{sample}_AF.png")

        if dp:
            # DP limitado
            save_plot(dp, f"DP — {sample}", plot_dp, xmax=200)
        if af:
            save_plot(af, f"AF — {sample}", plot_af)

        row = {
            "sample": sample,
            "snps": metrics["snps"],
            "indels": metrics["indels"],
            "titv_ratio": metrics["titv_ratio"],
            "mean_DP": pd.Series(dp).mean() if dp else None,
            "median_DP": pd.Series(dp).median() if dp else None,
            "mean_AF": pd.Series(af).mean() if af else None,
            "median_AF": pd.Series(af).median() if af else None,
        }
        results.append(row)

    df = pd.DataFrame(results)

    if not df.empty:
        # Ordenação alfabética
        df.sort_values(by="sample", inplace=True)

    csv_path = os.path.join(output_dir, "qc_panel.csv")
    df.to_csv(csv_path, index=False)
    print(f"[OK] CSV salvo em {csv_path}")

    if not df.empty:
        generate_html_report(output_dir, df)

if __name__ == "__main__":
    main()
