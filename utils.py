import os
from datetime import datetime

def get_file_creation_time(file_path):
    """
    Extrai o tempo em segundos a partir de um nome de arquivo no formato "..._HH_MM_SS_MMM.txt"
    :param file_path: Caminho do arquivo.
    :return: Data de criação do arquivo.
    """
    try:
        # Pega o nome do arquivo sem a extensão
        nome_arq = os.path.splitext(os.path.basename(file_path))[0]
        # Divide o nome do arquivo em partes
        partes = nome_arq.split('_')
        # Extrai as partes de hora, minuto, segundo e milissegundo
        hora_str, minuto_str, segundo_str, milissegundo_str = partes[-4:]
        # Constrói a string de hora no formato "HH:MM:SS:MMM"
        hora = f"{hora_str}:{minuto_str}:{segundo_str}.{milissegundo_str}"
        # Converte para datetime
        data_criacao = datetime.strptime(hora, "%H:%M:%S.%f").time()
        # Calcula o total de segundos
        total_segundos = data_criacao.hour * 3600 + data_criacao.minute * 60 + data_criacao.second + data_criacao.microsecond / 1_000_000

        return total_segundos

    except (IndexError, ValueError) as e:
        print(f"[AVISO] Não foi possível processar o tempo do arquivo: {file_path}")
        print(f"  > Erro: {e}")
        return None