# Nome do arquivo em excel
excel = "analise0411.xlsx"
# Pasta com o arquivo PBS para calcular o LOD
# Exemplo: "C:/Users/Usuario/Documents/PBS/"
pbs_path = "C:/Users/Victor/OneDrive/Área de Trabalho/IL1_04-11-25/04_11_25/PBS"
# Pasta com o anticorpo e o BSA
# Exemplo: "C:/Users/Usuario/Documents/Anticorpo/"
antibodies_path = "C:/Users/Victor/OneDrive/Área de Trabalho/IL1_04-11-25/04_11_25/ANT"
# Prefixo da amostra usada no experimento
# Exemplo: "IL1", "IL10", etc.
sample_prefix = "IL1"
# Configuração do filtro de espectro
filter_config = {
    "type": "savitzky", # Tipo de filtro: "savitzky", "media_movel", "mediana",...
    "params": {
        "window": 53, # Esse é o tamanho da janela para o filtro escolhido, quanto maior mais suave fica o sinal
        "order": 2 # Ordem do polinômio para o filtro escolhido, quanto maior mais próximo do sinal original fica o sinal
    }
}

