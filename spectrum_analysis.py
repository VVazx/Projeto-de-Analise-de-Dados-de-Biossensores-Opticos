import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from utils import get_file_creation_time
import config as cfg
from process_spectra import funcs
import custom_filters as cf


class SpectrumProcessor:
    """
    Classe para processar pastas contendo arquivos de espectro.
    Encapsula a lógica de carregamento, filtragem e extração de comprimento de onda ressonante.
    """

    def __init__(self, rp=3, dwl=2):
        """
        Inicializa o processador com parâmentros padrão de filtragem e carrega a configuração de filtro do config.py.
        :param rp: Resolução de proximidade (resolution_proximity).
        :param dwl: dwl (passado para get_approximate_valley).
        """
        self.rp = rp
        self.dwl = dwl
        try:
            self.filter_settings = cfg.filter_config 
            f_type = self.filter_settings.get("type", "desconhecido")
            logging.info(f"Objeto SpectrumProcessor criado. Filtro {f_type} ativo.")
                
        except AttributeError as e:
            logging.error(f"A variável 'filter_config' não foi encontrada no config.py. Erro: {e}")
            raise
    
    
    def _process_folder(self, folder_path, time_threshold_seconds, prominence):
        """
        Processa todos os arquivos em uma pasta.
        :param folder_path: Caminho da pasta contendo os arquivos de espectro.
        :param time_threshold_seconds:Tempo mínimo em segundos para considerar o dado.
        :param prominence: Prominência mínima para detecção de vales.
        :return: Lista de comprimentos de onda ressonantes extraídos.
        """
        logging.info(f"Processando pasta e subpastas: {folder_path}")
        if not os.path.isdir(folder_path):
            logging.error(f"Pasta não encontrada: {folder_path}")
            raise FileNotFoundError(f"Pasta não encontrada: {folder_path}")
        all_files = [] #lista para armazenar todos os arquivos .txt na pasta
        for dirpath, dirnames, filenames in os.walk(folder_path):
            for f in filenames:
                if f.lower().endswith('.txt'):
                    # Adiciona o caminho completo (ex: pasta/subpasta/arquivo.txt)
                    all_files.append(os.path.join(dirpath, f))
        if not all_files:
            logging.warning(f"Nenhum arquivo .txt encontrado em: {folder_path} (ou em suas subpastas)")
            return [], []
        logging.info(f"Encontrados {len(all_files)} arquivos .txt. Iniciando processamento...")
        spec_data = []
        timestamps_raw = []

        for f in all_files:
            try:
                creation_time = get_file_creation_time(f) # Obtém o tempo de criação em segundos
                if creation_time is None:
                    logging.warning(f"Arquivo com nome de formato inválido ignorado: {f}")
                    continue
                data = np.loadtxt(f, delimiter = ';')

                mask = (data[:,0] >= 1500) & (data[:,0] <= 1600) # Filtra o intervalo de interesse
                if np.any(mask):
                    timestamps_raw.append(creation_time)
                    spec_data.append(data[mask])
            
            except Exception as e:
                logging.warning(f"Erro ao processar o arquivo {f}: {e}")

        if not spec_data:
            logging.error(f"Nenhum dado válido encontrado na pasta.{folder_path}")

        timestamps_raw = np.array(timestamps_raw)
        timestamps_normalized = timestamps_raw - timestamps_raw.min()
        final_wavelengths = []
        final_timestamps = []

        for i, timestamp in enumerate(timestamps_normalized):
            if timestamp >= time_threshold_seconds:
                try:
                    # --- Roteador de Filtros (VERSÃO NOVA) ---
                    f_config = self.filter_settings
                    f_type = f_config["type"]
                    f_params = f_config["params"] # Dicionário de parâmetros
                    
                    spec_filtered_data = None
                    
                    try:
                        if f_type == "savitzky":
                            spec_filtered_data = cf.apply_savitzky(
                                spec_data[i], 
                                window=f_params["window"], 
                                order=f_params["order"]
                            )
                        
                        elif f_type == "media_movel":
                            spec_filtered_data = cf.apply_media_movel(
                                spec_data[i], 
                                window=f_params["window"]
                            )
                        
                        elif f_type == "mediana":
                            spec_filtered_data = cf.apply_mediana(
                                spec_data[i], 
                                window=f_params["window"]
                            )
                            
                        elif f_type == "nenhum":
                            spec_filtered_data = cf.apply_nenhum(spec_data[i])
                            
                        else:
                            logging.warning(f"Tipo de filtro '{f_type}' desconhecido. Nenhum filtro aplicado.")
                            spec_filtered_data = spec_data[i]

                    except Exception as e:
                        logging.error(f"Erro ao aplicar o filtro '{f_type}': {e}")
                        spec_filtered_data = spec_data[i]

                    info = {} # Dicionário para 'get_approximate_valley'
                    _, info = funcs.get_approximate_valley(
                        spec_filtered_data, 
                        info, 
                        prominence=prominence, 
                        resolution_proximity=self.rp, 
                        p0=None, 
                        dwl=self.dwl
                    )
                    
                    if 'resonant_wl' in info and info['resonant_wl'] is not None:
                        final_wavelengths.append(info['resonant_wl'])
                        final_timestamps.append(timestamp)
                
                except Exception as e:
                    logging.error(f"Erro ao processar espectro do arquivo {i} (tempo {timestamp}s): {e}")

        logging.info(f"Processamento da pasta concluído. {len(final_wavelengths)} pontos de dados extraídos.")
        return final_timestamps, final_wavelengths
    

    def run_antibody_bsa_analysis(self, folder_path):
        """
        Executa a análise do Anticorpo/BSA, processa os dados e gera o gráfico.
        :param folder_path: Caminho da pasta contendo os arquivos de espectro.
        """
        logging.info("Iniciando análise Anticorpo/BSA...")
        try:
            new_times, wavelengths = self._process_folder( # Chama o método privado para processar a pasta
                folder_path=folder_path, 
                time_threshold_seconds=400, 
                prominence=1
            )

            if not new_times:
                logging.warning("Nenhum dado válido extraído para análise Anticorpo/BSA.")
                return
            # Plotagem dos resultados
            plt.figure(figsize=(10, 6))
            sns.scatterplot(x=new_times, y=wavelengths, marker='.', label="Anticorpo/BSA")
            plt.xlabel('Tempo (s)')
            plt.ylabel('Comprimento de Onda (nm)')
            plt.title('Análise de Anticorpo / BSA')
            plt.grid(True, linestyle="--", linewidth=0.5)
            plt.legend()
            plt.tight_layout()
            plt.show()
            
            logging.info("Análise Anticorpo/BSA concluída com sucesso.")
        
        except Exception as e:
            logging.error(f"Erro durante a análise Anticorpo/BSA: {e}")
            raise


    def calculate_lod(self, folder_path):
        """
        Executa a análise do PBS (branco) e calcula o LOD a partir dos dados processados na pasta.
        :param folder_path: Caminho da pasta contendo os arquivos de espectro.
        :return: Valor do LOD calculado.
        """
        logging.info("Iniciando cálculo do LOD...")
        try:
            
            _, wavelengths = self._process_folder( # Chama o método privado para processar a pasta
                folder_path=folder_path,
                time_threshold_seconds=1,  
                prominence=0.5             
            )
            
            if not wavelengths:
                logging.error("Nenhum dado de PBS (branco) encontrado. Não é possível calcular o LOD.")
                return None
            # Cálculo do LOD
            xb1 = np.mean(wavelengths)
            sb1 = np.std(wavelengths)  
            lambda_LOD = xb1 + 3 * sb1

            logging.info(f"Cálculo do LOD: Média(Xb1) = {xb1:.4f} nm")
            logging.info(f"Cálculo do LOD: Desvio Padrão(sigma) = {sb1:.4f} nm")
            logging.info(f"Limite de Detecção (LOD): {lambda_LOD:.4f} nm")
            logging.info("Cálculo do LOD Concluído com sucesso.")
            return lambda_LOD
        
        except Exception as e:
            logging.error(f"Erro durante o cálculo do LOD: {e}")
            raise

# --- Bloco de Teste ---
# Para testar este módulo isoladamente
if __name__ == "__main__":
    
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


    try:
        import config as cfg
        pbs_folder = cfg.pbs_path
        antibody_folder = cfg.antibodies_path
    except (ImportError, AttributeError):
        logging.error("Arquivo config.py não encontrado ou não contém 'pbs_path' e 'antibodies_path'.")
        pbs_folder = "C:/Users/Victor/OneDrive/Área de Trabalho/IL1_04-11-25/04_11_25/PBS"
        antibody_folder = "C:/Users/Victor/OneDrive/Área de Trabalho/IL1_04-11-25/04_11_25/ANT"
        logging.warning(f"Usando caminhos de fallback para teste: {pbs_folder}")

    logging.info("Testando o módulo spectrum_analysis.py...")
    
    try:
        processor = SpectrumProcessor(rp=3, dwl=2)
        processor.run_antibody_bsa_analysis(folder_path=antibody_folder)
        lod_calculado = processor.calculate_lod(folder_path=pbs_folder)
        if lod_calculado:
            logging.info(f"LOD calculado no teste: {lod_calculado:.4f} nm")
            
    except Exception as e:
        logging.error(f"Falha no teste do módulo: {e}")