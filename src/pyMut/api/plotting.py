"""
API de gráficos para pyMut.

Este módulo proporciona funciones de alto nivel para crear visualizaciones
de resumen a partir de datos de mutaciones genéticas.
"""

import pandas as pd
import matplotlib.pyplot as plt
from typing import List, Dict, Union, Optional, Tuple

from ..utils import data_processing
from ..visualizations.summary import create_summary_plot

def plot_summary(data: Union[pd.DataFrame, str],
               figsize: Tuple[int, int] = (12, 10),
               title: str = "Resumen de mutaciones") -> plt.Figure:
    """
    Crea un gráfico de resumen con múltiples visualizaciones de los datos de mutaciones.
    
    Args:
        data: DataFrame con datos de mutaciones o ruta a un archivo TSV.
        figsize: Tamaño de la figura.
        title: Título del gráfico.
        
    Returns:
        Figura con las visualizaciones de resumen.
    """
    # Si data es una cadena, asumimos que es una ruta a un archivo TSV
    if isinstance(data, str):
        data = data_processing.read_tsv(data)
    
    # Preprocesar los datos para asegurar que tenemos las columnas necesarias
    processed_data = data_processing.extract_variant_classifications(
        data, 
        variant_column="Variant_Classification",
        funcotation_column="FUNCOTATION"
    )
    
    processed_data = data_processing.extract_variant_types(
        processed_data,
        variant_column="Variant_Type",
        funcotation_column="FUNCOTATION"
    )
    
    # Crear el gráfico de resumen
    return create_summary_plot(processed_data, figsize, title) 