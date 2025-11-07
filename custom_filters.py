import numpy as np
from scipy.signal import medfilt
from process_spectra import funcs # Este arquivo ainda precisa do especialista original

"""
Este módulo centraliza todas as implementações de filtros de espectro.
Ele atua como uma "biblioteca" interna para o SpectrumProcessor.
"""

def apply_savitzky(spec_data, window, order):
    """
    Aplica o filtro Savitzky-Golay (da biblioteca 'process_spectra').
    """
    try:
        filtered_data = funcs.filter_spectrum(spec_data, None, window, order, quiet=True)[0]
        return filtered_data
    except Exception as e:
        print(f"Erro no filtro Savitzky-Golay: {e}")
        return spec_data 

def apply_media_movel(spec_data, window):
    """
    Aplica um filtro de Média Móvel simples.
    """
    try:
        y_data = spec_data[:, 1]
        y_filtered = np.convolve(y_data, np.ones(window)/window, mode='same')
        filtered_spec = np.copy(spec_data)
        filtered_spec[:, 1] = y_filtered
        return filtered_spec
    except Exception as e:
        print(f"Erro no filtro Média Móvel: {e}")
        return spec_data

def apply_mediana(spec_data, window):
    """
    Aplica um filtro de Mediana (ótimo para "spikes").
    """
    try:
        if window % 2 == 0:
            window += 1 
        y_data = spec_data[:, 1]
        y_filtered = medfilt(y_data, kernel_size=window)
        filtered_spec = np.copy(spec_data)
        filtered_spec[:, 1] = y_filtered
        return filtered_spec
    except Exception as e:
        print(f"Erro no filtro Mediana: {e}")
        return spec_data

def apply_nenhum(spec_data):
    """
    Não faz nada, apenas retorna os dados originais.
    """
    return spec_data

# Caso necessário, mais filtros podem ser adicionados aqui seguindo o mesmo padrão.
