import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit, root_scalar
from sklearn.metrics import r2_score, mean_squared_error
import re
import logging
import config as cfg 

class CalibrationAnalyzer:
    """
    Classe para analizar as curvas de calibração
    """
    def __init__(self, excel_path, sample_prefix):
        """
        Inicializa a classe com o caminho do arquivo Excel.
        :param excel_path: Caminho do arquivo Excel.
        :param sample_prefix: Amostra usada no experimento.
        """
        logging.info(f"Objeto CalibrationAnalyzer criado. Excel: {excel_path}")
        self.excel_path = excel_path 
        self.sample_prefix = sample_prefix # Amostra utilizada no experimento
        self.df_raw = None # Dataframe brutos carregados do Excel
        self.df_processed = None # Dataframe após limpeza e agrupamento
        self.x_raw_pg = None # Dados de concentração (pg/mL) - para PLOTS
        self.x_data_ng = None # Dados de concentração (ng/mL) - para FIT
        self.y_data = None # Dados de comprimento de onda (nm)
        self.y_error = None # Desvio padrão dos dados de comprimento de onda
        self.models_to_fit = []
        self.fit_results = {}
        self.best_model_name = None
        self.best_model_details = None


    def _convert_concentration(self, concentration):
        """
        Converte a concentração de string para float em pg/mL.
        :param concentration: Concentração em string (ex: "10ng").
        :return: Concentração em float (pg/mL).
        """
        try:
            number = float(concentration[:-2].replace(',', '.'))
            unit = concentration[-2:]
            factor = {'pg': 1, 'ng': 1e3, 'ug': 1e6}.get(unit, 1)
            return number * factor
        except Exception as e:
            logging.warning(f"Não foi possível converter concentração: {concentration}. Erro: {e}")
            return None


    def _logistic_inv(self, x, A, B, C, D):
        """
        Modelo logístico de 4 parâmetros.
        :param x: Concentração.
        :param A, B, C, D: Parâmetros do modelo.
        :return: Valor calculado pelo modelo.
        """
        return A / (1 + np.exp(B * (x - C))) + D
    

    def _model_equation(self, name, params):
        """
        Formata a equação do modelo como uma string para o gráfico.
        :param name: Nome do modelo.
        :param params: Parâmetros do modelo.
        :return: Equação formatada como string.
        """
        try: 
            if name == 'Logístico':
                return f'y = {params[0]:.2f} / (1 + np.exp({params[1]:.2f} * (x - {params[2]:.2f}))) + {params[3]:.2f}'
            elif name == 'Exponencial':
                return f'y = {params[0]:.2f} * np.exp({params[1]:.2f} * x) + {params[2]:.2f}'
            elif name == 'Polinomial 2ª ordem':
                return f'y = {params[0]:.2f} * x² + {params[1]:.2f} * x + {params[2]:.2f}'
            elif name == 'Polinomial 3ª ordem':
                return f'y = {params[0]:.2f} * x³ + {params[1]:.2f} * x² + {params[2]:.2f} * x + {params[3]:.2f}'
            else:
                return 'Equação não disponível'
        except Exception as e:
            logging.warning(f"Não foi possível formatar a equação para {name}: {e}")
            return "Erro na formatação da equação"


    def load_and_prepare_data(self):
        """
        Carrega os dados do arquivo Excel e prepara os dados para análise.
        """
        logging.info(f"Carregando o arquivo Excel: {self.excel_path}...")
        
        try:
            self.df_raw = pd.read_excel(self.excel_path, decimal=',')
            # Construção de um dataframe com os dados do excel.
            df = self.df_raw.copy()
            logging.info(f"Usando a amostra: {self.sample_prefix} para extrair as concentrações.")

            # Extração da concentração da amostra usando regex
            regex_pattern = rf'{self.sample_prefix}_(\d+[,\.]?\d*(?:ng|pg|ug))'
            def extract_conc(sample_name):
                match = re.findall(regex_pattern, sample_name)
                return match[0] if match else None
            df['Concentração_str'] = df['amostra'].apply(extract_conc)
            df = df[df['Concentração_str'].notna()]

            # Conversão de pg/mL para float
            df['Concentração'] = df['Concentração_str'].apply(self._convert_concentration)
            # Agrupamento dos dados por concentração
            calib_m = df.groupby('Concentração')['comprimento_onda_ressonante'].mean()
            calib_s = df.groupby('Concentração')['comprimento_onda_ressonante'].std()

            self.df_processed= calib_m.to_frame('Mean_WL')
            self.df_processed['Std_WL'] = calib_s

            self.x_raw_pg= calib_m.index.values # Concentração em pg/mL (para PLOT)
            self.x_data_ng = self.x_raw_pg / 1000 # Concentração em ng/mL (para FIT)
            self.y_data= calib_m.values # Comprimento de onda em nm
            self.y_error= calib_s.values # Desvio padrão do comprimento de onda
            
            logging.info("Dados carregados e preparados com sucesso.")
        except FileNotFoundError:
            logging.error(f"Arquivo não encontrado: {self.excel_path}")
            raise
        except Exception as e:
            logging.error(f"Erro ao carregar ou preparar os dados: {e}")
            raise


    def define_models(self):
        """
        Define os modelos a serem ajustados.
        """
        if self.y_data is None or self.x_data_ng is None: # Verificação se os dados foram carregados
            logging.error("Dados não carregados. Execute 'load_and_prepare_data()' primeiro.")
            raise ValueError("Dados não carregados.")
        
        self.models_to_fit = [ # Lista de modelos a serem ajustados
            {
                'name': 'Logístico',
                'function': self._logistic_inv,
                'p0': [max(self.y_data)-min(self.y_data), 1, np.median(self.x_data_ng), min(self.y_data)]
            }
        ] # Adicione outros modelos conforme necessário
        logging.info(f"{len(self.models_to_fit)} modelo(s) definido(s) para ajuste.")

    
    def fit_models(self):
        """
        Executa o 'curve_fit' para todos os modelos definidos.
        """
        if not self.models_to_fit:
            logging.error("Nenhum modelo definido para ajuste. Execute 'define_models()' primeiro.")
            raise ValueError("Nenhum modelo definido para ajuste.")
        
        logging.info("Iniciando ajuste de curvas...")
        # Dicionário para armazenar os resultados dos ajustes
        for model in self.models_to_fit:
            try:
                params, pcov= curve_fit(model['function'], self.x_data_ng, self.y_data, p0=model['p0'], maxfev=10000)
                y_fit= model['function'](self.x_data_ng, *params)
                r2= r2_score(self.y_data, y_fit)
                mse = mean_squared_error(self.y_data, y_fit)
                self.fit_results[model['name']] = {
                    'params': params,
                    'r2': r2,
                    'mse': mse,
                    'function': model['function']
                }
                logging.info(f"Ajuste concluído para o modelo {model['name']}: R²={r2:.4f}, MSE={mse:.4f}")
            
            except Exception as e:
                logging.error(f"Erro ao ajustar o modelo {model['name']}: {str(e)}")
                self.fit_results[model['name']] = None

        
        logging.info("\n--- Comparação dos Modelos de Ajuste ---") # Exibição dos resultados
        logging.info(f"{'Modelo':<20} {'R²':<10} {'MSE':<10}")
        for name, result in self.fit_results.items():
            if result is not None:
                logging.info(f"{name:<20} {result['r2']:<10.4f} {result['mse']:<15.4f}")

       
        try:
            self.best_model_name = max( # Encontra o melhor modelo com base no R²
                (name for name in self.fit_results if self.fit_results[name] is not None),
                key=lambda name: self.fit_results[name]['r2']
                )
            self.best_model_details = self.fit_results[self.best_model_name]
            logging.info(f"Melhor modelo selecionado: {self.best_model_name} com R²={self.best_model_details['r2']:.4f}")   
        except ValueError:
            logging.error("Nenhum modelo foi ajustado com sucesso. Não é possível determinar o melhor modelo.")
            raise Exception("Falha em todos os ajustes de modelo.")


    def plot_results(self):
        """
        Plota os resultados dos ajustes (comparação).
        """
        if not self.fit_results:
            logging.error("Nenhum resultado de ajuste disponível. Execute 'fit_models()' primeiro.")
            return 
            
        logging.info("Plotando a comparação de todos os modelos...")
        # Gráfico de comparação de todos os modelos
        plt.figure(figsize=(10, 6))
        plt.errorbar(self.x_raw_pg, self.y_data, yerr=2*self.y_error, fmt='ok', label='Dados Experimentais')
        x_fit_pg = np.logspace(np.log10(min(self.x_raw_pg)*0.8), np.log10(max(self.x_raw_pg)*1.2), 300)
        x_fit_ng = x_fit_pg / 1e3 # Conversão para ng/mL
        colors = ['r-','b-','g-','m-']
        for (name, result), color in zip(self.fit_results.items(), colors):
            if result is not None:
                plt.plot(x_fit_pg, result['function'](x_fit_ng, *result['params']), color, 
                         label=f'{name} (R²={result["r2"]:.3f})')
        plt.xscale('log')
        plt.xlabel('Concentração (pg/mL)')
        plt.ylabel('Comprimento de Onda (nm)')
        plt.title('Comparação de Modelos de Ajuste')
        plt.grid(True, which='both', ls='--', lw=0.5)
        plt.legend()
        plt.tight_layout()
        plt.show()


    def plot_best_fit(self):
        """
        Plota o melhor ajuste do modelo.
        """
        if not self.best_model_name:
            logging.warning("Nenhum melhor modelo disponível. Não é possível plotar.")
            return
        
        logging.info(f"Plotando o melhor ajuste: {self.best_model_name}...")

        plt.figure(figsize=(8, 5))
        plt.errorbar(self.x_raw_pg, self.y_data, yerr=2*self.y_error, fmt='ok', label='Dados Experimentais')
        #Curva do melhor ajuste
        best_func = self.best_model_details['function']
        best_params = self.best_model_details['params']
        best_r2 = self.best_model_details['r2']
        equation_str = self._model_equation(self.best_model_name, best_params)

        x_fit_pg = np.logspace(np.log10(min(self.x_raw_pg)*0.8), np.log10(max(self.x_raw_pg)*1.2), 300)
        x_fit_ng = x_fit_pg / 1e3 # Conversão para ng/mL

        plt.plot(x_fit_pg, best_func(x_fit_ng, *best_params), 
                 'r-', label=f'{self.best_model_name}\n{equation_str}\nR² = {best_r2 :.4f}')
        plt.xscale('log')
        plt.xlabel('Concentração (pg/mL)')
        plt.ylabel('Comprimento de Onda (nm)')
        plt.title(f'Melhor Ajuste: {self.best_model_name}')
        plt.grid(True, which='both', ls='--', lw=0.5)
        plt.legend()
        plt.tight_layout()
        plt.show()


    def predict_concentration(self, y_input, x_bounds=(1e-5, 1e3)):
        """
        Prevê a concentração com base no comprimento de onda usando o melhor modelo.
        :param y_input: Comprimento de onda medido.
        :param x_bounds: Limites (em ng/mL) para o solver procurar a raiz.
        :return: Concentração prevista (em pg/mL).
        """
        if not self.best_model_name:
            logging.error("Nenhum melhor modelo disponível. Execute 'fit_models()' primeiro.")
            return None
            
        best_func = self.best_model_details['function']
        best_params = self.best_model_details['params']
        
        # Função para prever a concentração (solver busca a raiz de f(x) = 0)
        def to_solve(x_ng):
            return best_func(x_ng, *best_params) - y_input
            
        try:
            sol = root_scalar(to_solve, bracket=x_bounds, method='brentq')
            if sol.converged:
                conc_ng = sol.root
                logging.info(f"Concentração prevista para {y_input} nm: {conc_ng*1e3:.2f} pg/mL")
                return conc_ng * 1e3 # Converter de ng/mL para pg/mL
            else:
                logging.warning(f"A solução não convergiu para o comprimento de onda {y_input} nm.")
                return None
        except Exception as e:
            logging.error(f"Erro ao prever concentração para {y_input} nm (verifique 'x_bounds'): {e}")
            return None


    def run_full_analysis(self):
        """
        Executa toda a análise de calibração.
        """
        logging.info("Iniciando Análise de Calibração Completa...")
        try:
            self.load_and_prepare_data()
            self.define_models()
            self.fit_models()
            self.plot_results()
            self.plot_best_fit()
            logging.info("Análise de Calibração Concluída com Sucesso.")
        
        except Exception as e:
            logging.error(f"ERRO na análise de calibração: {e}")
            raise e
        
# --- Bloco de Teste ---
# Apenas para testar o módulo independentemente
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    try:
        excel_file_path = cfg.excel
        sample_prefix = cfg.sample_prefix
    except (NameError, AttributeError):
        logging.error("Arquivo config.py não encontrado ou não contém a variável 'excel'.")
        exit()

    logging.info(f"Testando o módulo calibration_analysis.py com o arquivo {excel_file_path}...")
    try:
        analyzer = CalibrationAnalyzer(excel_path=excel_file_path, sample_prefix=sample_prefix)
        analyzer.run_full_analysis()
        wl_teste = 1550.5 # Mude este valor para testar
        conc_prevista = analyzer.predict_concentration(wl_teste) 
        if conc_prevista is not None:
            logging.info(f"Previsão para {wl_teste} nm: {conc_prevista:.2f} pg/mL")
            
    except Exception as e:
        logging.error(f"Falha no teste do módulo: {e}")