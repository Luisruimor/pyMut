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
                   title: str = "Resumen de mutaciones",
                   show_interactive: bool = False) -> plt.Figure:
        """
        Genera un gráfico de resumen con estadísticas generales de las mutaciones.
        
        Esta visualización incluye múltiples gráficos:
        - Variant Classification: Distribución de clasificaciones de variantes
        - Variant Type: Distribución de tipos de variantes (SNP, INS, DEL, etc.)
        - SNV Class: Distribución de clases de SNV (cambios nucleotídicos como A>G, C>T, etc.)
        - Variants per Sample: Distribución de variantes por muestra y mediana (TMB)
        
        Args:
            figsize: Tamaño de la figura.
            title: Título del gráfico.
            show_interactive: Si es True, muestra la visualización en modo interactivo.
            
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
        fig = create_summary_plot(processed_data, figsize, title)
        
        # Si se solicita mostrar interactivamente
        if show_interactive:
            plt.show()
        
        return fig
    
    def variant_classification_plot(self,
                                    figsize: Tuple[int, int] = (8, 6),
                                    title: str = "Variant Classification",
                                    show_interactive: bool = False) -> plt.Figure:
        """
        Genera un gráfico de barras horizontal mostrando la distribución de clasificaciones de variantes.
        
        Args:
            figsize: Tamaño de la figura.
            title: Título del gráfico.
            show_interactive: Si es True, muestra la visualización en modo interactivo.
            
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
        
        # Si se solicita mostrar interactivamente
        if show_interactive:
            plt.show()
        
        return fig
    
    def variant_type_plot(self,
                          figsize: Tuple[int, int] = (8, 6),
                          title: str = "Variant Type",
                          show_interactive: bool = False) -> plt.Figure:
        """
        Genera un gráfico de barras horizontal mostrando la distribución de tipos de variantes.
        
        Args:
            figsize: Tamaño de la figura.
            title: Título del gráfico.
            show_interactive: Si es True, muestra la visualización en modo interactivo.
            
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
        
        # Si se solicita mostrar interactivamente
        if show_interactive:
            plt.show()
        
        return fig
    
    def snv_class_plot(self,
                        figsize: Tuple[int, int] = (8, 6),
                        title: str = "SNV Class",
                        ref_column: str = "REF",
                        alt_column: str = "ALT",
                        show_interactive: bool = False) -> plt.Figure:
        """
        Genera un gráfico de barras horizontal mostrando la distribución de clases de SNV.
        
        Args:
            figsize: Tamaño de la figura.
            title: Título del gráfico.
            ref_column: Nombre de la columna que contiene el alelo de referencia.
            alt_column: Nombre de la columna que contiene el alelo alternativo.
            show_interactive: Si es True, muestra la visualización en modo interactivo.
            
        Returns:
            Figura de matplotlib con el gráfico de clases de SNV.
        """
        from .visualizations.summary import create_snv_class_plot
        
        fig, ax = plt.subplots(figsize=figsize)
        create_snv_class_plot(
            self.data, 
            ref_column=ref_column,
            alt_column=alt_column,
            ax=ax
        )
        if title:
            fig.suptitle(title, fontsize=16)
        plt.tight_layout()
        
        # Si se solicita mostrar interactivamente
        if show_interactive:
            plt.show()
        
        return fig

    def variants_per_sample_plot(self,
                                 figsize: Tuple[int, int] = (10, 6),
                                 title: str = "Variants per Sample",
                                 variant_column: str = "Variant_Classification",
                                 sample_column: str = "Tumor_Sample_Barcode",
                                 show_interactive: bool = False) -> plt.Figure:
        """
        Genera un gráfico de barras apiladas mostrando el número de variantes por muestra (TMB)
        y su composición por tipo de variante.

        Args:
            figsize: Tamaño de la figura.
            title: Título del gráfico.
            variant_column: Nombre de la columna que contiene la clasificación de variante.
            sample_column: Nombre de la columna que contiene el identificador de la muestra,
                          o string que se usará para identificar columnas de muestra si las
                          muestras están como columnas.
            show_interactive: Si es True, muestra la visualización en modo interactivo.
            
        Returns:
            Figura de matplotlib con el gráfico de variantes por muestra.
        """
        from .visualizations.summary import create_variants_per_sample_plot
        from .utils.data_processing import extract_variant_classifications

        # Si variant_column no está en las columnas, intentar normalizarlo
        if variant_column not in self.data.columns:
            # Comprobar si hay una versión con diferente capitalización
            column_lower = variant_column.lower()
            for col in self.data.columns:
                if col.lower() == column_lower:
                    variant_column = col
                    break
        
        # Asegurar que la columna de clasificación de variantes existe o se extrae
        processed_data = extract_variant_classifications(
            self.data, 
            variant_column=variant_column,
            funcotation_column="FUNCOTATION"
        )
        
        fig, ax = plt.subplots(figsize=figsize)
        create_variants_per_sample_plot(
            processed_data, 
            variant_column=variant_column,
            sample_column=sample_column,
            ax=ax
        )
        
        # Si se proporciona un título personalizado diferente al título por defecto con mediana
        if title and not title.startswith("Variants per Sample (Median:"):
            fig.suptitle(title, fontsize=16, y=1.02)
        
        plt.tight_layout()
        
        # Si se solicita mostrar interactivamente
        if show_interactive:
            plt.show()
        
        return fig
