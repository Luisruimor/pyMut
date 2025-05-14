"""
Módulo para la generación de gráficos de resumen.

Este módulo contiene funciones para crear visualizaciones de resumen
que muestran diferentes estadísticas de los datos de mutaciones.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import re
from typing import List, Dict, Union, Optional, Tuple, Any


def create_variant_classification_plot(data: pd.DataFrame,
                                     variant_column: str = "Variant_Classification",
                                     ax: Optional[plt.Axes] = None) -> plt.Axes:
    """
    Crea un diagrama de barras horizontal mostrando el recuento por cada tipo de clasificación de variante.
    
    Args:
        data: DataFrame con los datos de mutaciones.
        variant_column: Nombre de la columna que contiene la clasificación de variante.
        ax: Eje de matplotlib donde dibujar. Si es None, se crea uno nuevo.
        
    Returns:
        Eje de matplotlib con la visualización.
    """
    # Contar las variantes por clasificación
    variant_counts = data[variant_column].value_counts().to_dict()
    
    # Ordenar los recuentos
    variant_counts = dict(sorted(variant_counts.items(), key=lambda item: item[1]))
    
    # Crear el eje si no se proporcionó uno
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 6))
    
    # Crear la gráfica de barras horizontal
    cmap = plt.colormaps['tab20']  # Obtén el colormap
    colores = cmap(range(len(variant_counts)))  # Genera los colores necesarios
    
    bars = ax.barh(list(variant_counts.keys()), list(variant_counts.values()), color=colores)
    
    # Ajustar título y etiquetas
    ax.set_title("Variant Classification", fontsize=14)
    ax.set_ylabel("")  # Se elimina el título del eje Y
    
    # Añadir etiquetas con el recuento en cada barra
    for bar in bars:
        ax.text(bar.get_width() + 10,
                 bar.get_y() + bar.get_height()/2,
                 f'{int(bar.get_width())}',
                 va='center', fontsize=10)
    
    # Ajustar márgenes y quitar recuadros innecesarios
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    return ax


def create_variant_type_plot(data: pd.DataFrame,
                           variant_column: str = "Variant_Type",
                           ax: Optional[plt.Axes] = None) -> plt.Axes:
    """
    Crea un diagrama de barras horizontal mostrando el recuento por cada tipo de variante.
    
    Args:
        data: DataFrame con los datos de mutaciones.
        variant_column: Nombre de la columna que contiene el tipo de variante.
        ax: Eje de matplotlib donde dibujar. Si es None, se crea uno nuevo.
        
    Returns:
        Eje de matplotlib con la visualización.
    """
    # Contar las variantes por tipo
    variant_counts = data[variant_column].value_counts().to_dict()
    
    # Ordenar los recuentos
    variant_counts = dict(sorted(variant_counts.items(), key=lambda item: item[1]))
    
    # Crear el eje si no se proporcionó uno
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 6))
    
    # Colores específicos (si el número de categorías excede la lista, se utilizará un colormap)
    colors = ['#D3C4E7', '#FFFACD', '#87CEEB']
    if len(variant_counts) > len(colors):
        cmap = plt.colormaps['tab20']
        colors = cmap(range(len(variant_counts)))
    else:
        colors = colors[:len(variant_counts)]
    
    # Crear la gráfica de barras horizontal
    bars = ax.barh(list(variant_counts.keys()), list(variant_counts.values()), color=colors)
    
    # Ajustar título y etiquetas
    ax.set_title("Variant Type", fontsize=14)
    ax.set_ylabel("")  # Se elimina el título del eje Y
    
    # Configuración de las etiquetas del eje Y en cursiva
    for tick in ax.get_yticklabels():
        tick.set_fontstyle('italic')
    
    # Añadir etiquetas con el recuento en cada barra
    for bar in bars:
        ax.text(bar.get_width() + 10,
                 bar.get_y() + bar.get_height()/2,
                 f'{int(bar.get_width())}',
                 va='center', fontsize=10)
    
    # Ajustar márgenes y quitar recuadros innecesarios
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    return ax


def create_snv_class_plot(data: pd.DataFrame,
                        ref_column: str = "REF",
                        alt_column: str = "ALT",
                        ax: Optional[plt.Axes] = None) -> plt.Axes:
    """
    Crea un diagrama de barras horizontal mostrando el recuento por cada clase de SNV.
    
    Args:
        data: DataFrame con los datos de mutaciones.
        ref_column: Nombre de la columna que contiene el alelo de referencia.
        alt_column: Nombre de la columna que contiene el alelo alternativo (tumoral).
        ax: Eje de matplotlib donde dibujar. Si es None, se crea uno nuevo.
        
    Returns:
        Eje de matplotlib con la visualización.
    """
    # Crear una copia del DataFrame para no modificar el original
    df_copy = data.copy()
    
    # Verificar si podemos generar la información de SNV Class
    if ref_column in df_copy.columns and alt_column in df_copy.columns:
        # Generar la SNV Class combinando el alelo de referencia y el alelo alternativo
        df_copy['SNV_Class'] = df_copy[ref_column] + '>' + df_copy[alt_column]
        # Filtrar para mantener solo los cambios de nucléotido único (un solo carácter a cada lado)
        df_copy = df_copy[df_copy['SNV_Class'].str.match(r'^[A-Z]>[A-Z]$')]
    else:
        # Si no podemos generar la información, devolver un eje vacío
        print(f"No se pueden generar las clases de SNV: faltan las columnas {ref_column} y/o {alt_column}.")
        if ax is None:
            _, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "No hay datos disponibles para SNV Class", 
               ha='center', va='center', fontsize=12)
        ax.set_title("SNV Class", fontsize=14)
        ax.axis('off')
        return ax
    
    # Eliminar filas sin información de SNV Class
    df_copy = df_copy[df_copy['SNV_Class'].notnull()]
    
    # Contar las variantes por clase de SNV
    snv_counts = df_copy['SNV_Class'].value_counts().to_dict()
    
    # Si no hay datos, mostrar un mensaje
    if not snv_counts:
        if ax is None:
            _, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "No hay datos disponibles para SNV Class", 
               ha='center', va='center', fontsize=12)
        ax.set_title("SNV Class", fontsize=14)
        ax.axis('off')
        return ax
    
    # Ordenar los recuentos
    snv_counts = dict(sorted(snv_counts.items(), key=lambda item: item[1]))
    
    # Crear el eje si no se proporcionó uno
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 6))
    
    # Definir colores específicos para cada clase
    colors = ['#FF8C00', '#9ACD32', '#FFD700', '#FF4500', '#4169E1', '#1E90FF']
    if len(snv_counts) > len(colors):
        cmap = plt.colormaps['tab20']
        colors = cmap(range(len(snv_counts)))
    else:
        colors = colors[:len(snv_counts)]
    
    # Crear la gráfica de barras horizontal
    bars = ax.barh(list(snv_counts.keys()), list(snv_counts.values()), color=colors)
    
    # Ajustar título y etiquetas
    ax.set_title("SNV Class", fontsize=14)
    ax.set_ylabel("")  # Se elimina el título del eje Y
    
    # Añadir etiquetas con el recuento en cada barra
    for bar in bars:
        ax.text(bar.get_width() + 10,
                 bar.get_y() + bar.get_height()/2,
                 f'{int(bar.get_width())}',
                 va='center', fontsize=10, fontstyle='italic')
    
    # Ajustar márgenes y quitar recuadros innecesarios
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    return ax


def create_summary_plot(data: pd.DataFrame,
                      figsize: Tuple[int, int] = (12, 10),
                      title: str = "Resumen de mutaciones") -> plt.Figure:
    """
    Crea un gráfico de resumen con múltiples visualizaciones de los datos de mutaciones.
    
    Args:
        data: DataFrame con los datos de mutaciones.
        figsize: Tamaño de la figura.
        title: Título del gráfico.
        
    Returns:
        Figura con las visualizaciones de resumen.
    """
    # from ..utils.data_processing import extract_variant_classifications, extract_variant_types
    
    # Crear una figura con múltiples subplots
    fig, axs = plt.subplots(2, 3, figsize=figsize)
    fig.suptitle(title, fontsize=16)
    
    # Preprocesar los datos para asegurar que tenemos las columnas necesarias
    # processed_data = extract_variant_classifications(
    # data, 
    # variant_column="Variant_Classification",
    # funcotation_column="FUNCOTATION"
    # )
    
    # processed_data = extract_variant_types(
    # processed_data,
    # variant_column="Variant_Type",
    # funcotation_column="FUNCOTATION"
    # )
    
    # Crear el gráfico de clasificación de variantes
    create_variant_classification_plot(data, ax=axs[0, 0])
    
    # Crear el gráfico de tipos de variantes
    create_variant_type_plot(data, ax=axs[0, 1])
    
    # Crear el gráfico de clases de SNV
    create_snv_class_plot(data, 
                         ref_column="REF",
                         alt_column="ALT",
                         ax=axs[0, 2])
    
    # TODO: Implementar los demás gráficos de resumen
    # Por ejemplo:
    # - Distribución de mutaciones por cromosoma
    # - Distribución de mutaciones por muestra
    # - Distribución de efectos en proteínas
    
    # Ajustar espaciado entre subplots
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)  # Ajustar para dejar espacio para el título principal
    
    return fig 