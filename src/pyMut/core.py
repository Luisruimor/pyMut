import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Dict, Union, Optional, Tuple

class PyMutation:
    """
    Clase principal para visualizar mutaciones genéticas a partir de datos en formato TSV.
    
    Esta clase proporciona métodos para generar visualizaciones de resumen de datos
    de mutaciones genéticas, mostrando estadísticas generales como distribuciones
    de tipos de variantes, clasificaciones y cambios de nucleótidos.
    
    Attributes:
        data (pd.DataFrame): DataFrame que contiene los datos de mutaciones.
    """
    
    def __init__(self, data: pd.DataFrame):
        """
        Inicializa un objeto PyMutation con un DataFrame de pandas.
        
        Args:
            data (pd.DataFrame): DataFrame conteniendo datos de mutaciones.
        """
        self.data = data
    
    def summary_plot(self, 
                   figsize: Tuple[int, int] = (12, 10),
                   title: str = "Resumen de mutaciones") -> plt.Figure:
        """
        Genera un gráfico de resumen con estadísticas generales de las mutaciones.
        
        Esta visualización incluye múltiples gráficos:
        - Variant Classification: Distribución de clasificaciones de variantes
        - Variant Type: Distribución de tipos de variantes (SNP, INS, DEL, etc.)
        - SNV Class: Distribución de clases de SNV (cambios nucleotídicos como A>G, C>T, etc.)
        
        Args:
            figsize: Tamaño de la figura.
            title: Título del gráfico.
            
        Returns:
            Figura de matplotlib con el gráfico de resumen.
        """
        from .visualizations.summary import create_summary_plot
        from .utils.data_processing import extract_variant_classifications, extract_variant_types
        
        # Preprocesar los datos para asegurar que tenemos las columnas necesarias
        processed_data = extract_variant_classifications(
            self.data, 
            variant_column="Variant_Classification",
            funcotation_column="FUNCOTATION"
        )
        
        processed_data = extract_variant_types(
            processed_data,
            variant_column="Variant_Type",
            funcotation_column="FUNCOTATION"
        )
        
        # Generar el gráfico de resumen
        return create_summary_plot(processed_data, figsize, title)
