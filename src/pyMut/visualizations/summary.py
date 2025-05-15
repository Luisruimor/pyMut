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
                                     ax: Optional[plt.Axes] = None,
                                     color_map: Optional[Dict] = None) -> plt.Axes:
    """
    Crea un diagrama de barras horizontal mostrando el recuento por cada tipo de clasificación de variante.
    
    Args:
        data: DataFrame con los datos de mutaciones.
        variant_column: Nombre de la columna que contiene la clasificación de variante.
        ax: Eje de matplotlib donde dibujar. Si es None, se crea uno nuevo.
        color_map: Diccionario opcional que mapea las clasificaciones de variantes a colores.
        
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
    if color_map:
        colors = [color_map.get(variant, cm.get_cmap('tab20')(i % 20)) for i, variant in enumerate(variant_counts.keys())]
    else:
        cmap = plt.colormaps['tab20']  # Obtén el colormap
        colors = cmap(range(len(variant_counts)))  # Genera los colores necesarios
    
    bars = ax.barh(list(variant_counts.keys()), list(variant_counts.values()), color=colors)
    
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
                 va='center', fontstyle='italic')
    
    # Ajustar márgenes y quitar recuadros innecesarios
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    return ax


def create_variant_classification_summary_plot(data: pd.DataFrame,
                                             variant_column: str = "Variant_Classification",
                                             sample_column: str = "Tumor_Sample_Barcode",
                                             ax: Optional[plt.Axes] = None,
                                             color_map: Optional[Dict] = None) -> plt.Axes:
    """
    Crea un diagrama de cajas y bigotes (boxplot) que resume, para cada clasificación de variantes,
    la distribución (entre las muestras) del número de alelos alternativos detectados.

    Args:
        data: DataFrame con los datos de mutaciones.
        variant_column: Nombre de la columna que contiene la clasificación de variante.
        sample_column: Nombre de la columna que contiene el identificador de la muestra.
                       Si no existe, se asume que las muestras son columnas (formato ancho).
        ax: Eje de matplotlib donde dibujar. Si es None, se crea uno nuevo.
        color_map: Diccionario opcional que mapea las clasificaciones de variantes a colores.
        
    Returns:
        Eje de matplotlib con la visualización.
    """
    # Verificar si tenemos la columna de clasificación
    if variant_column not in data.columns:
        print(f"No se encontró la columna: {variant_column}")
        if ax is None:
            _, ax = plt.subplots(figsize=(12, 6))
        ax.text(0.5, 0.5, f"No hay datos disponibles para Variant Classification Summary\nFalta columna: {variant_column}", 
               ha='center', va='center', fontsize=12)
        ax.set_title("Variant Classification Summary", fontsize=14)
        ax.axis('off')
        return ax
    
    # Detectar si tenemos formato largo o ancho
    samples_as_columns = sample_column not in data.columns
    
    if samples_as_columns:
        # Formato ancho: las muestras son columnas
        # Encontrar columnas que podrían ser muestras (IDs de TCGA, etc.)
        potential_sample_cols = [col for col in data.columns if col.startswith('TCGA-') or 
                                (isinstance(col, str) and col.count('-') >= 2)]
        
        if not potential_sample_cols:
            print(f"No se encontraron columnas de muestra que parezcan identificadores")
            if ax is None:
                _, ax = plt.subplots(figsize=(12, 6))
            ax.text(0.5, 0.5, "No hay datos disponibles para Variant Classification Summary\nNo se detectaron columnas de muestra", 
                  ha='center', va='center', fontsize=12)
            ax.set_title("Variant Classification Summary", fontsize=14)
            ax.axis('off')
            return ax
            
        print(f"Detectadas {len(potential_sample_cols)} columnas de muestra")
        
        # Acumular los conteos por muestra y clasificación
        sample_variant_counts = {}
        
        # Obtener valores únicos de variantes para agrupar
        unique_variants = data[variant_column].unique()
        unique_variants = [v for v in unique_variants if pd.notna(v) and v != "Unknown"]
        
        # Primero determinar el formato de los valores en las columnas de muestra
        # Revisamos algunas filas para ver si son del tipo "A|B"
        sample_col = potential_sample_cols[0]
        sample_format = "unknown"
        sample_values = data[sample_col].dropna().values[:10]  # Tomar algunas muestras
        
        if len(sample_values) > 0:
            # Verificar si hay valores en formato "A|B" o "A/B"
            if any('|' in str(x) for x in sample_values):
                sample_format = "pipe_separated"
            elif any('/' in str(x) for x in sample_values):
                sample_format = "slash_separated"
            else:
                sample_format = "other"
                
        print(f"Formato detectado de muestras: {sample_format}")
        
        # Para cada clasificación de variante, contar cuántas muestras tienen esa variante
        for variant_class in unique_variants:
            # Filtrar filas con esta clasificación de variante
            variant_subset = data[data[variant_column] == variant_class]
            
            # Para cada muestra, contar cuántas variantes de este tipo tiene
            for sample in potential_sample_cols:
                # Contar según el formato detectado
                if sample_format == "pipe_separated":
                    # Para formato "A|B", contar cuando A y B son diferentes (alelo variante)
                    sample_count = variant_subset[sample].apply(
                        lambda x: 1 if (isinstance(x, str) and '|' in x and 
                                       x.split('|')[0] != x.split('|')[1]) else 0
                    ).sum()
                elif sample_format == "slash_separated":
                    # Similar para formato "A/B"
                    sample_count = variant_subset[sample].apply(
                        lambda x: 1 if (isinstance(x, str) and '/' in x and 
                                       x.split('/')[0] != x.split('/')[1]) else 0
                    ).sum()
                else:
                    # Para otros formatos, asumimos que valores distintos de "0", "0/0" u "0|0" indican variante
                    sample_count = variant_subset[sample].apply(
                        lambda x: 1 if (x != 0 and x != '0' and x != '0/0' and 
                                       x != '0|0' and not pd.isnull(x)) else 0
                    ).sum()
                
                if sample not in sample_variant_counts:
                    sample_variant_counts[sample] = {}
                
                if sample_count > 0:
                    sample_variant_counts[sample][variant_class] = sample_count
    else:
        # Formato largo: hay una columna específica para la muestra
        # Acumular los conteos por muestra y clasificación
        sample_variant_counts = {}
        
        # Procesar el DataFrame para contar por muestra y variante
        for _, row in data.iterrows():
            sample = row[sample_column]
            variant_class = row[variant_column]
            
            if pd.isnull(variant_class) or variant_class == "Unknown":
                continue
                
            if sample not in sample_variant_counts:
                sample_variant_counts[sample] = {}
                
            if variant_class not in sample_variant_counts[sample]:
                sample_variant_counts[sample][variant_class] = 0
                
            # Incrementar el contador (asumiendo 1 alelo alternativo por variante)
            sample_variant_counts[sample][variant_class] += 1
    
    # Convertir el diccionario a DataFrame (filas: muestras, columnas: variant classifications)
    df_data = {sample: dict(variants) for sample, variants in sample_variant_counts.items()}
    df = pd.DataFrame.from_dict(df_data, orient='index').fillna(0)
    
    if df.empty:
        print("No se encontraron datos para el análisis después del procesamiento.")
        if ax is None:
            _, ax = plt.subplots(figsize=(12, 6))
        ax.text(0.5, 0.5, "No hay datos disponibles para Variant Classification Summary\nNo se encontraron datos para analizar", 
               ha='center', va='center', fontsize=12)
        ax.set_title("Variant Classification Summary", fontsize=14)
        ax.axis('off')
        return ax
    
    # (Opcional) Reordenar las columnas según la suma total (de mayor a menor)
    col_order = df.sum(axis=0).sort_values(ascending=False).index.tolist()
    df = df[col_order]
    
    # Eliminar columnas con suma total 0
    df = df.loc[:, df.sum() > 0]
    
    # Preparar los datos para el boxplot: cada columna (variant classification) es una serie de conteos por muestra
    variant_types = df.columns.tolist()
    data_to_plot = [df[vt].values for vt in variant_types]
    
    # Crear el eje si no se proporcionó uno
    if ax is None:
        _, ax = plt.subplots(figsize=(12, 6))
    
    # Dibujar el boxplot con configuración mejorada
    bp = ax.boxplot(data_to_plot, 
                    patch_artist=True,  # Rellenar las cajas con color
                    medianprops=dict(color="red", linewidth=1.5),
                    boxprops=dict(linewidth=1.5),
                    whiskerprops=dict(linewidth=1.5),
                    capprops=dict(linewidth=1.5),
                    flierprops=dict(marker='o', markerfacecolor='gray', markersize=4, alpha=0.5),
                    showfliers=True,  # Mostrar outliers
                    widths=0.7)
    
    # Colorear cada caja usando una paleta de colores o el mapa de colores proporcionado
    if color_map:
        colors = [color_map.get(vt, plt.colormaps['tab20'](i % 20)) for i, vt in enumerate(variant_types)]
    else:
        # Usar un colormap de alta diferenciación para las cajas
        cmap = plt.colormaps['tab20']
        colors = [cmap(i % 20) for i in range(len(variant_types))]
    
    # Aplicar colores a cada caja
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    # Configurar las etiquetas del eje X con mejor rotación y formato
    ax.set_xticklabels(variant_types, rotation=45, ha='right', fontsize=10)
    
    # Ajustar los límites del eje Y
    ymin = 0
    ymax = max(max(d) if len(d) > 0 else 0 for d in data_to_plot) * 1.1
    if ymax == 0:  # Por si no hay datos positivos
        ymax = 1
    ax.set_ylim(ymin, ymax)
    
    # Mejorar la presentación visual
    # Añadir título más descriptivo
    ax.set_title("Variant Classification Summary", fontsize=14, fontweight='bold')
    ax.set_xlabel("Variant Classification", fontsize=12, fontweight='bold')
    ax.set_ylabel("Número de alelos alternativos por muestra", fontsize=12, fontweight='bold')
    
    # Añadir cuadrícula para mejor legibilidad
    ax.yaxis.grid(True, linestyle='--', alpha=0.3)
    
    # Eliminar algunos bordes de los ejes
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Ajustar espacio para las etiquetas del eje X
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    
    return ax


def create_summary_plot(data: pd.DataFrame,
                      figsize: Tuple[int, int] = (16, 12),
                      title: str = "Resumen de Mutaciones") -> plt.Figure:
    """
    Crea un gráfico de resumen con múltiples visualizaciones de los datos de mutaciones.
    
    Args:
        data: DataFrame con los datos de mutaciones.
        figsize: Tamaño de la figura.
        title: Título del gráfico.
        
    Returns:
        Figura con las visualizaciones de resumen.
    """
    # Crear una figura con múltiples subplots, haciendo más anchos los gráficos
    fig, axs = plt.subplots(2, 3, figsize=figsize, gridspec_kw={'width_ratios': [1.5, 1.5, 1.5], 'height_ratios': [1, 1.2]})
    fig.suptitle(title, fontsize=16)
    
    # Detectar nombres de columnas respetando capitalización original
    variant_classification_col = "Variant_Classification"
    sample_column = "Tumor_Sample_Barcode"
    
    # Buscar variantes de capitalización para variant_classification
    if variant_classification_col not in data.columns:
        for col in data.columns:
            if col.lower() == variant_classification_col.lower():
                variant_classification_col = col
                break
    
    # Generar un mapa de colores coherente para todas las clasificaciones de variantes
    unique_variants = data[variant_classification_col].unique()
    # Usar un colormap fijo para asegurar colores consistentes
    cmap = plt.colormaps['tab20']  # Colormap con buena variedad de colores
    variant_color_map = {variant: cmap(i % 20) for i, variant in enumerate(unique_variants) if pd.notna(variant)}
    
    # Crear el gráfico de clasificación de variantes usando el mapa de colores predefinido
    var_class_ax = create_variant_classification_plot(
        data, 
        variant_column=variant_classification_col, 
        ax=axs[0, 0],
        color_map=variant_color_map  # Pasar el mapa de colores
    )
    
    # Crear el gráfico de tipos de variantes
    create_variant_type_plot(data, ax=axs[0, 1])
    
    # Crear el gráfico de clases de SNV
    create_snv_class_plot(data, 
                         ref_column="REF",
                         alt_column="ALT",
                         ax=axs[0, 2])
    
    # Crear el gráfico de variantes por muestra (TMB)
    # Pasando el mismo mapa de colores que usamos para el gráfico de clasificación
    variants_ax = create_variants_per_sample_plot(
        data,
        variant_column=variant_classification_col,
        sample_column=sample_column,
        ax=axs[1, 0],
        variant_colors=variant_color_map  # Usar el mismo mapa de colores
    )
    
    # Crear el gráfico de variant classification summary
    var_boxplot_ax = create_variant_classification_summary_plot(
        data,
        variant_column=variant_classification_col,
        sample_column=sample_column,
        ax=axs[1, 1],
        color_map=variant_color_map  # Usar el mismo mapa de colores
    )
    
    # El último panel puede quedarse vacío o agregar información adicional
    axs[1, 2].axis('off')  # Desactivar el último panel por ahora
    
    # Quitar las leyendas individuales para evitar duplicación
    if var_class_ax.get_legend() is not None:
        var_class_ax.get_legend().remove()
    
    if variants_ax.get_legend() is not None:
        variants_ax.get_legend().remove()
        
    if var_boxplot_ax.get_legend() is not None:
        var_boxplot_ax.get_legend().remove()
    
    # Crear una leyenda común para los gráficos y colocarla en la parte inferior
    # Crear handles y labels manualmente para la leyenda global
    handles = []
    labels = []
    for variant, color in variant_color_map.items():
        if pd.isnull(variant) or variant == "Unknown":
            continue
        patch = plt.Rectangle((0,0), 1, 1, fc=color)
        handles.append(patch)
        labels.append(variant)
    
    # Añadir la leyenda común
    fig.legend(handles, labels, loc='lower center', ncol=min(len(labels), 5), 
              title=variant_classification_col, bbox_to_anchor=(0.5, 0.01))
    
    # Ajustar espaciado entre subplots
    plt.tight_layout()
    plt.subplots_adjust(top=0.9, bottom=0.18)  # Ajustar para leyenda común, dándole más espacio
    
    return fig


def create_variants_per_sample_plot(data: pd.DataFrame,
                                   variant_column: str = "Variant_Classification",
                                   sample_column: str = "Tumor_Sample_Barcode",
                                   ax: Optional[plt.Axes] = None,
                                   variant_colors: Optional[Dict] = None) -> plt.Axes:
    """
    Crea un gráfico de barras apiladas mostrando el número de variantes por muestra (TMB)
    y su composición por tipo de variante.

    Args:
        data: DataFrame con los datos de mutaciones. Puede tener las muestras como columnas
              o una columna específica con los identificadores de muestra.
        variant_column: Nombre de la columna que contiene la clasificación de variante.
        sample_column: Nombre de la columna que contiene el identificador de la muestra,
                      o string que se usará para identificar columnas de muestra si las 
                      muestras están como columnas.
        ax: Eje de matplotlib donde dibujar. Si es None, se crea uno nuevo.
        variant_colors: Diccionario opcional que mapea las clasificaciones de variantes a colores.
        
    Returns:
        Eje de matplotlib con la visualización.
    """
    # Verificar si tenemos las columnas necesarias
    if variant_column not in data.columns:
        print(f"No se encontró la columna: {variant_column}")
        if ax is None:
            _, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, f"No hay datos disponibles para Variants per Sample\nFalta columna: {variant_column}", 
               ha='center', va='center', fontsize=12)
        ax.set_title("Variants per Sample", fontsize=14)
        ax.axis('off')
        return ax
    
    # Detectar el formato de los datos
    # Si sample_column existe como columna, asumimos formato "largo"
    # Si no, asumimos formato "ancho" donde las muestras son columnas
    samples_as_columns = sample_column not in data.columns
    
    if samples_as_columns:
        # Encontrar columnas que podrían ser muestras (IDs de TCGA, etc.)
        potential_sample_cols = [col for col in data.columns if col.startswith('TCGA-') or '|' in str(data[col].iloc[0])]
        
        if not potential_sample_cols:
            print(f"No se encontraron columnas de muestra que empiecen con TCGA- o contengan el carácter |")
            if ax is None:
                _, ax = plt.subplots(figsize=(10, 6))
            ax.text(0.5, 0.5, "No hay datos disponibles para Variants per Sample\nNo se detectaron columnas de muestra", 
                  ha='center', va='center', fontsize=12)
            ax.set_title("Variants per Sample", fontsize=14)
            ax.axis('off')
            return ax
        
        # Convertir el formato ancho a largo para conteo
        # Primero, contar variantes por tipo para cada muestra
        variant_counts = {}
        
        # Agrupar por Variant_Classification
        for variant_type in data[variant_column].unique():
            variant_subset = data[data[variant_column] == variant_type]
            
            for sample_col in potential_sample_cols:
                # Si la columna de muestra contiene los genotipos en formato REF|ALT
                sample_variants = variant_subset[sample_col].apply(
                    lambda x: 1 if '|' in str(x) and str(x).split('|')[0] != str(x).split('|')[1] else 0
                ).sum()
                
                if sample_col not in variant_counts:
                    variant_counts[sample_col] = {}
                
                if sample_variants > 0:
                    variant_counts[sample_col][variant_type] = sample_variants
        
        # Convertir el diccionario a DataFrame
        if not variant_counts:
            print("No se encontraron variantes para ninguna muestra")
            if ax is None:
                _, ax = plt.subplots(figsize=(10, 6))
            ax.text(0.5, 0.5, "No hay datos disponibles para Variants per Sample\nNo se detectaron variantes", 
                  ha='center', va='center', fontsize=12)
            ax.set_title("Variants per Sample", fontsize=14)
            ax.axis('off')
            return ax
        
        # Crear DataFrame a partir del diccionario
        samples_df = []
        for sample, variants in variant_counts.items():
            for var_type, count in variants.items():
                samples_df.append({'Sample': sample, 'Variant_Classification': var_type, 'Count': count})
        
        if not samples_df:
            print("No se pudieron procesar las variantes por muestra")
            if ax is None:
                _, ax = plt.subplots(figsize=(10, 6))
            ax.text(0.5, 0.5, "No hay datos disponibles para Variants per Sample\nNo se pudieron procesar las variantes", 
                  ha='center', va='center', fontsize=12)
            ax.set_title("Variants per Sample", fontsize=14)
            ax.axis('off')
            return ax
        
        processed_df = pd.DataFrame(samples_df)
        variant_counts = processed_df.pivot(index='Sample', columns='Variant_Classification', values='Count').fillna(0)
        
    else:
        # Formato "largo" donde hay una columna con la muestra
        # Contar variantes por muestra y por clasificación
        variant_counts = data.groupby([sample_column, variant_column]).size().unstack(fill_value=0)
    
    # Calcular el total de variantes por muestra y ordenar de mayor a menor
    variant_counts['total'] = variant_counts.sum(axis=1)
    variant_counts = variant_counts.sort_values('total', ascending=False)
    
    # Calcular la mediana de variantes por muestra
    median_tmb = variant_counts['total'].median()
    
    # Eliminar la columna de total para la visualización
    if 'total' in variant_counts.columns:
        variant_counts = variant_counts.drop('total', axis=1)
    
    # Crear el eje si no se proporcionó uno
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))
    
    # Generar colores para las diferentes clasificaciones de variantes
    if variant_colors is not None:
        # Si se proporcionan colores específicos, usarlos para las variantes correspondientes
        colors = []
        for variant in variant_counts.columns:
            # Buscar el color en el mapa de colores para esta variante exacta
            if variant in variant_colors:
                colors.append(variant_colors[variant])
            else:
                # Si no se encuentra el color exacto, usar uno predeterminado
                cmap = cm.get_cmap('tab20', len(variant_counts.columns))
                colors.append(cmap(len(colors) % cmap.N))
    else:
        # Usar el mapa de colores predeterminado
        cmap = cm.get_cmap('tab20', len(variant_counts.columns))
        colors = [cmap(i) for i in range(len(variant_counts.columns))]
    
    # Crear la gráfica de barras apiladas
    variant_counts.plot(kind='bar', stacked=True, ax=ax, color=colors, width=0.8)
    
    # Añadir una línea horizontal para la mediana
    ax.axhline(y=median_tmb, color='red', linestyle='--', linewidth=1)
    
    # Configurar etiquetas y título (incluir la mediana en el título)
    ax.set_title(f"Variants per Sample\nMedian: {median_tmb:.1f}", fontsize=14)
    ax.set_xlabel("Samples", fontsize=12)
    ax.set_ylabel("Number of Variants", fontsize=12)
    
    # Quitar las etiquetas del eje X para mayor claridad cuando hay muchas muestras
    ax.set_xticklabels([])
    ax.tick_params(axis='x', which='both', bottom=False)
    
    # Ajustar leyenda
    ax.legend(title=variant_column, bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Ajustar márgenes y quitar recuadros innecesarios
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    return ax 