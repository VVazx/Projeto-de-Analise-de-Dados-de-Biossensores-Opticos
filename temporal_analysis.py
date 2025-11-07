import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import logging

def run_temporal_analysis(excel_path):
    """
    Carrega o arquivo Excel e realiza uma análise temporal dos dados.
    Plota o gráfico de dispersão (scatterplot) do comprimento de onda ressonante ao longo do tempo.
    :param excel_path: Caminho do arquivo Excel.
    """
    logging.info("Iniciando análise temporal...")
    try:
        # Carrega os dados do Excel
        logging.info(f"Carregando o arquivo Excel: {excel_path}")
        df = pd.read_excel(excel_path)

        if "horario" not in df.columns or "comprimento_onda_ressonante" not in df.columns: #Verifica se as colunas existem
            logging.error("Colunas necessárias não encontradas no arquivo Excel.")
            raise KeyError("Colunas 'horario' ou 'comprimento_onda_ressonante' não encontradas.")

        df["Tempo (s)"] = (df["horario"] - df["horario"].min()).dt.total_seconds() # Normaliza o tempo em segundos

        df = df.sort_values(by="horario") # Ordena os dados pelo tempo

        logging.info("Gerando gráfico de análise temporal...") # Gera o gráfico de dispersão
        plt.figure(figsize=(12, 6))
        sns.scatterplot(
            data=df,
            x="Tempo (s)",
            y="comprimento_onda_ressonante",
            hue="amostra",
            style="amostra",
            s=100,
            palette="tab10"
        )
        plt.xlabel("Tempo (s)")
        plt.ylabel("Comprimento de Onda (nm)")
        plt.title("Comprimento de Onda ao Longo do Tempo")
        plt.grid(True, linestyle="--", linewidth=0.5)
        plt.legend(title="Soluções", bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.show()

        logging.info("Análise temporal concluída com sucesso.")

    except Exception as e:
        logging.error(f"Erro durante a análise temporal: {e}")
        raise

# --- Bloco de Teste ---
# Para testar este módulo isoladamente
if __name__ == "__main__":
    
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    try:
        import config as cfg
        excel_file_path = cfg.excel
    except (ImportError, AttributeError):
        logging.error("Arquivo config.py não encontrado ou não contém a variável 'excel'.")
        excel_file_path = "analisedois.xlsx" # Mude se necessário
        logging.warning(f"Usando caminho de fallback para teste: {excel_file_path}")
    logging.info(f"Testando o módulo temporal_analysis.py com o arquivo {excel_file_path}...")
    try:
        run_temporal_analysis(excel_path=excel_file_path)
    except Exception as e:
        logging.error(f"Falha no teste do módulo: {e}")