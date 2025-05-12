from typing import List, Optional
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Dict, Union, Optional, Tuple

class MutationMetadata:
    """
    Clase para almacenar metadatos de mutaciones.

    Atributos:
        source_format (str): Formato de origen (VCF, MAF, etc.).
        file_path (str): Ruta del archivo de origen.
        loaded_at (datetime): Fecha y hora de carga.
        filters (List[str]): Filtros aplicados al archivo.
        fasta (str): Ruta del archivo FASTA.
        notes (Optional[str]): Notas adicionales.
    """
    def __init__(self, source_format: str, file_path: str,
                 filters: List[str],fasta: str, notes: Optional[str] = None):
        self.source_format = source_format
        self.file_path = file_path
        self.loaded_at = datetime.now()
        self.filters = filters
        self.notes = notes
        self.fasta = fasta


class PyMutation:
    def __init__(self, data: pd.DataFrame, metadata: MutationMetadata):
        self.data = data
        self.metadata = metadata
    
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
    
    def variant_classification_plot(self,
                                  figsize: Tuple[int, int] = (8, 6),
                                  title: str = "Variant Classification") -> plt.Figure:
        """
        Genera un gráfico de barras horizontal mostrando la distribución de clasificaciones de variantes.
        
        Args:
            figsize: Tamaño de la figura.
            title: Título del gráfico.
            
        Returns:
            Figura de matplotlib con el gráfico de clasificación de variantes.
        """
        from .visualizations.summary import create_variant_classification_plot
        from .utils.data_processing import extract_variant_classifications
        
        # Preprocesar los datos para asegurar que tenemos la columna necesaria
        processed_data = extract_variant_classifications(
            self.data, 
            variant_column="Variant_Classification",
            funcotation_column="FUNCOTATION"
        )
        
        # Crear figura y ejes
        fig, ax = plt.subplots(figsize=figsize)
        
        # Generar el gráfico
        create_variant_classification_plot(processed_data, ax=ax)
        
        # Configurar título
        if title:
            fig.suptitle(title, fontsize=16)
        
        plt.tight_layout()
        return fig
    
    def variant_type_plot(self,
                        figsize: Tuple[int, int] = (8, 6),
                        title: str = "Variant Type") -> plt.Figure:
        """
        Genera un gráfico de barras horizontal mostrando la distribución de tipos de variantes.
        
        Args:
            figsize: Tamaño de la figura.
            title: Título del gráfico.
            
        Returns:
            Figura de matplotlib con el gráfico de tipos de variantes.
        """
        from .visualizations.summary import create_variant_type_plot
        from .utils.data_processing import extract_variant_types
        
        # Preprocesar los datos para asegurar que tenemos la columna necesaria
        processed_data = extract_variant_types(
            self.data,
            variant_column="Variant_Type",
            funcotation_column="FUNCOTATION"
        )
        
        # Crear figura y ejes
        fig, ax = plt.subplots(figsize=figsize)
        
        # Generar el gráfico
        create_variant_type_plot(processed_data, ax=ax)
        
        # Configurar título
        if title:
            fig.suptitle(title, fontsize=16)
        
        plt.tight_layout()
        return fig
    
    def snv_class_plot(self,
                      figsize: Tuple[int, int] = (8, 6),
                      title: str = "SNV Class",
                      ref_column: str = "REF",
                      alt_column: str = "ALT") -> plt.Figure:
        """
        Genera un gráfico de barras horizontal mostrando la distribución de clases de SNV.
        
        Args:
            figsize: Tamaño de la figura.
            title: Título del gráfico.
            ref_column: Nombre de la columna que contiene el alelo de referencia.
            alt_column: Nombre de la columna que contiene el alelo alternativo.
            
        Returns:
            Figura de matplotlib con el gráfico de clases de SNV.
        """
        from .visualizations.summary import create_snv_class_plot
        
        # Crear figura y ejes
        fig, ax = plt.subplots(figsize=figsize)
        
        # Generar el gráfico
        create_snv_class_plot(
            self.data, 
            ref_column=ref_column,
            alt_column=alt_column,
            ax=ax
        )
        
        # Configurar título
        if title:
            fig.suptitle(title, fontsize=16)
        
        plt.tight_layout()
        return fig
