import logging
import sys
try:
    import config as cfg
except ImportError:
    logging.error("Arquivo config.py não encontrado.")
    sys.exit(1)
try:
    from calibration_analysis import CalibrationAnalyzer
    from spectrum_analysis import SpectrumProcessor
    from temporal_analysis import run_temporal_analysis
except ImportError as e:
    logging.error(f"Erro ao importar um módulo de análise: {e}")
    print("Certifique-se de que todos os arquivos necessários estão na pasta usada.")
    sys.exit(1)

def main():
    try:
        analyzer = CalibrationAnalyzer(excel_path=cfg.excel, sample_prefix=cfg.sample_prefix)
        analyzer.run_full_analysis()
    except Exception as e:
        logging.error(f"Erro durante a análise de calibração: {e}")
        logging.error("A rotina será interrompida.")
        return
    try:
       
        run_temporal_analysis(excel_path=cfg.excel)
        
    except Exception as e:
        logging.error(f"Erro durante a análise temporal: {e}")
    try:
        logging.info("Iniciando o processo de análise de espectro...")
        processor = SpectrumProcessor()
        processor.run_antibody_bsa_analysis(folder_path=cfg.antibodies_path)
        logging.info("Análise de espectro concluída com sucesso.")
    except Exception as e:
        logging.error(f"Erro durante a análise de espectro: {e}")
    logging.info("Processo de análise completo.")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, 
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        handlers=[
                            logging.FileHandler("analise.log"),
                            logging.StreamHandler(sys.stdout)
                                ])

    main()

       