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
    
    # Colorear cada barra usando una paleta de colores
    if color_map:
        colors = [color_map.get(variant, plt.colormaps['tab20'](i % 20)) for i, variant in enumerate(variant_counts.keys())]
    else:
        # Usar el colormap 'tab20'
        cmap = plt.colormaps['tab20']  # En lugar de cm.get_cmap('tab20')
        colors = [cmap(i % 20) for i in range(len(variant_counts))]
    
    bars = ax.barh(list(variant_counts.keys()), list(variant_counts.values()), color=colors)
    
    # Ajustar título y etiquetas
    ax.set_title("Variant Classification", fontsize=14, fontweight='bold')
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
    ax.set_title("Variant Type", fontsize=14, fontweight='bold')
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
        ax.set_title("SNV Class", fontsize=14, fontweight='bold')
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
        ax.set_title("SNV Class", fontsize=14, fontweight='bold')
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
    ax.set_title("SNV Class", fontsize=14, fontweight='bold')
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
                                             color_map: Optional[Dict] = None,
                                             show_labels: bool = True) -> plt.Axes:
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
        show_labels: Si True, muestra las etiquetas de las clasificaciones en el eje X.
        
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
        ax.set_title("Variant Classification Summary", fontsize=14, fontweight='bold')
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
            ax.set_title("Variant Classification Summary", fontsize=14, fontweight='bold')
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
        ax.set_title("Variant Classification Summary", fontsize=14, fontweight='bold')
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
    if show_labels:
        ax.set_xticklabels(variant_types, rotation=45, ha='right', fontsize=10)
    else:
        ax.set_xticklabels([])  # No mostrar etiquetas
    
    # Ajustar los límites del eje Y
    ymin = 0
    ymax = max(max(d) if len(d) > 0 else 0 for d in data_to_plot) * 1.1
    if ymax == 0:  # Por si no hay datos positivos
        ymax = 1
    ax.set_ylim(ymin, ymax)
    
    # Mejorar la presentación visual
    # Añadir título más descriptivo
    ax.set_title("Variant Classification Summary", fontsize=14, fontweight='bold')
    
    # Añadir cuadrícula para mejor legibilidad
    ax.yaxis.grid(True, linestyle='--', alpha=0.3)
    
    # Eliminar algunos bordes de los ejes
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
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
    gene_column = "Hugo_Symbol"
    
    # Buscar variantes de capitalización para variant_classification
    if variant_classification_col not in data.columns:
        for col in data.columns:
            if col.lower() == variant_classification_col.lower():
                variant_classification_col = col
                break
    
    # Buscar variantes de capitalización para gene_column
    if gene_column not in data.columns:
        for col in data.columns:
            if col.lower() == gene_column.lower():
                gene_column = col
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
        color_map=variant_color_map  # Usar el mismo mapa de colores
    )
    
    # Crear el gráfico de variant classification summary
    var_boxplot_ax = create_variant_classification_summary_plot(
        data,
        variant_column=variant_classification_col,
        sample_column=sample_column,
        ax=axs[1, 1],
        color_map=variant_color_map,  # Usar el mismo mapa de colores
        show_labels=False  # No mostrar etiquetas en el summary plot
    )
    
    # Crear el gráfico de top mutated genes
    top_genes_ax = create_top_mutated_genes_plot(
        data,
        variant_column=variant_classification_col,
        gene_column=gene_column,
        sample_column=sample_column,
        mode="variants",
        count=10,
        ax=axs[1, 2],
        color_map=variant_color_map  # Usar el mismo mapa de colores
    )
    
    # Quitar las leyendas individuales para evitar duplicación
    if var_class_ax.get_legend() is not None:
        var_class_ax.get_legend().remove()
    
    if variants_ax.get_legend() is not None:
        variants_ax.get_legend().remove()
        
    if var_boxplot_ax.get_legend() is not None:
        var_boxplot_ax.get_legend().remove()
        
    if top_genes_ax.get_legend() is not None:
        top_genes_ax.get_legend().remove()
    
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
                                   color_map: Optional[Dict] = None) -> plt.Axes:
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
        color_map: Diccionario opcional que mapea las clasificaciones de variantes a colores.
        
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
        ax.set_title("Variants per Sample", fontsize=14, fontweight='bold')
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
            ax.set_title("Variants per Sample", fontsize=14, fontweight='bold')
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
            ax.set_title("Variants per Sample", fontsize=14, fontweight='bold')
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
    if color_map is not None:
        # Si se proporcionan colores específicos, usarlos para las variantes correspondientes
        colors = []
        for variant in variant_counts.columns:
            # Buscar el color en el mapa de colores para esta variante exacta
            if variant in color_map:
                colors.append(color_map[variant])
            else:
                # Si no se encuentra el color exacto, usar uno predeterminado
                # Usar plt.colormaps en lugar de cm.get_cmap
                colors.append(plt.colormaps['tab20'](len(colors) % 20))
    else:
        # Usar el mapa de colores predeterminado
        # Usar plt.colormaps en lugar de cm.get_cmap
        cmap = plt.colormaps['tab20']
        colors = [cmap(i % cmap.N) for i in range(len(variant_counts.columns))]
    
    # Crear la gráfica de barras apiladas
    variant_counts.plot(kind='bar', stacked=True, ax=ax, color=colors, width=0.8)
    
    # Añadir una línea horizontal para la mediana
    ax.axhline(y=median_tmb, color='red', linestyle='--', linewidth=1)
    
    # Configurar etiquetas y título (incluir la mediana en el título pero sin negrita)
    ax.set_title("Variants per Sample", fontsize=14, fontweight='bold')
    # Añadir la información de la mediana como subtítulo sin negrita
    ax.text(0.5, 0.92, f"Median: {median_tmb:.1f}", transform=ax.transAxes, ha='center', fontsize=12)
    ax.set_xlabel("")  # Eliminar la etiqueta del eje X "Samples"
    
    # Quitar las etiquetas del eje X para mayor claridad cuando hay muchas muestras
    ax.set_xticklabels([])
    ax.tick_params(axis='x', which='both', bottom=False)
    
    # Ajustar leyenda
    ax.legend(title=variant_column, bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Ajustar márgenes y quitar recuadros innecesarios
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    return ax


def create_top_mutated_genes_plot(data: pd.DataFrame,
                               mode: str = "variants",
                               variant_column: str = "Variant_Classification",
                               gene_column: str = "Hugo_Symbol",
                               sample_column: str = "Tumor_Sample_Barcode",
                               count: int = 10,
                               ax: Optional[plt.Axes] = None,
                               color_map: Optional[Dict] = None) -> plt.Axes:
    """
    Crea un diagrama de barras horizontal mostrando los genes más mutados y la distribución
    de variantes según su clasificación.

    Args:
        data: DataFrame con los datos de mutaciones.
        mode: Modo de conteo de mutaciones: "variants" (cuenta número total de variantes)
              o "samples" (cuenta número de muestras afectadas).
        variant_column: Nombre de la columna que contiene la clasificación de variante.
        gene_column: Nombre de la columna que contiene el símbolo del gen.
        sample_column: Nombre de la columna que contiene el identificador de la muestra,
                      o prefijo que se usará para identificar columnas de muestra si las
                      muestras están como columnas.
        count: Número de genes principales a mostrar.
        ax: Eje de matplotlib donde dibujar. Si es None, se crea uno nuevo.
        color_map: Diccionario opcional que mapea las clasificaciones de variantes a colores.
        
    Returns:
        Eje de matplotlib con la visualización.
    """
    # Verificar si tenemos las columnas necesarias
    if gene_column not in data.columns:
        print(f"No se encontró la columna: {gene_column}")
        if ax is None:
            _, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, f"No hay datos disponibles para Top Mutated Genes\nFalta columna: {gene_column}", 
               ha='center', va='center', fontsize=12)
        ax.set_title("Top Mutated Genes", fontsize=14, fontweight='bold')
        ax.axis('off')
        return ax
        
    if variant_column not in data.columns:
        print(f"No se encontró la columna: {variant_column}")
        if ax is None:
            _, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, f"No hay datos disponibles para Top Mutated Genes\nFalta columna: {variant_column}", 
               ha='center', va='center', fontsize=12)
        ax.set_title("Top Mutated Genes", fontsize=14, fontweight='bold')
        ax.axis('off')
        return ax
    
    # Detectar si las muestras están como columnas (formato ancho)
    # Buscar columnas que puedan ser muestras (por ejemplo, empiezan con TCGA-)
    sample_cols = [col for col in data.columns if str(col).startswith("TCGA-")]
    samples_as_columns = len(sample_cols) > 0
    
    if not samples_as_columns and sample_column not in data.columns:
        print(f"No se encontró la columna de muestras: {sample_column}")
        print("Tampoco se detectaron columnas de muestra con formato TCGA-*")
        if ax is None:
            _, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, f"No hay datos disponibles para Top Mutated Genes\nNo se detectaron muestras", 
               ha='center', va='center', fontsize=12)
        ax.set_title("Top Mutated Genes", fontsize=14, fontweight='bold')
        ax.axis('off')
        return ax
    
    # Crear el eje si no se proporcionó uno
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 8))
    
    # Filtrar filas con valores faltantes o "Unknown"
    data_filtered = data[(data[gene_column].notna()) & (data[variant_column].notna())]
    data_filtered = data_filtered[(data_filtered[gene_column] != "Unknown") & 
                                   (data_filtered[variant_column] != "Unknown")]
    
    if mode == "variants":
        # MODO "variants" - Contar el número total de variantes
        # Contar variantes por gen y tipo de variante
        gene_variant_counts = data_filtered.groupby([gene_column, variant_column]).size().unstack(fill_value=0)
        
        # Calcular el total para cada gen y ordenar
        gene_totals = gene_variant_counts.sum(axis=1).sort_values(ascending=False)
        top_genes = gene_totals.index[:count].tolist()
        
        # Seleccionar solo los top genes
        df_top = gene_variant_counts.loc[top_genes]
        
        # Ordenarlos manualmente de menor a mayor para visualización
        ordered_genes = sorted(top_genes, key=lambda g: gene_totals[g])
        df_plot = df_top.loc[ordered_genes]
        
        # Asignar colores para cada tipo de variante
        if color_map is None:
            cmap = plt.colormaps['tab20']  # En lugar de cm.get_cmap
            colors = [cmap(i % 20) for i in range(len(df_plot.columns))]
        else:
            colors = [color_map.get(variant, plt.colormaps['tab20'](i % 20)) 
                     for i, variant in enumerate(df_plot.columns)]
        
        # Crear gráfico de barras horizontales
        df_plot.plot(kind='barh', stacked=True, ax=ax, color=colors, width=0.65)
        
        # Añadir etiquetas con el recuento total a la derecha de cada barra
        for i, gene in enumerate(df_plot.index):
            total = gene_totals[gene]
            # Ajustar el offset para los números
            offset_variants = max(1, 0.01 * ax.get_xlim()[1]) if ax.get_xlim()[1] > 0 else 1
            ax.text(total + offset_variants, i, f'{int(total)}', va='center', fontsize=10)
        
        title_text_variants = f"Top {count} Mutated Genes (Total Variants)" # Usar count variable
        # Configuración común del título y ejes para el modo variants
        ax.set_title(title_text_variants, fontsize=14, fontweight='bold') 
        ax.set_ylabel("")  # Eliminar la etiqueta del eje Y
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False) # Ocultar línea del eje izquierdo
        ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=True) # Sin ticks en eje Y, pero con etiquetas

        # Mejorar la leyenda para el modo variants
        handles_v, labels_v = ax.get_legend_handles_labels()
        by_label_v = dict(zip(labels_v, handles_v))
        
        if by_label_v:
            num_legend_items_v = len(by_label_v)
            ncol_legend_v = min(num_legend_items_v, 4)
            base_offset_v = -0.20
            row_offset_factor_v = 0.06
            num_legend_rows_v = (num_legend_items_v + ncol_legend_v - 1) // ncol_legend_v
            vertical_offset_v = base_offset_v - (row_offset_factor_v * num_legend_rows_v)
            ax.legend(by_label_v.values(), by_label_v.keys(),
                      title="Variant Classification",
                      loc="lower center",
                      bbox_to_anchor=(0.5, vertical_offset_v),
                      ncol=ncol_legend_v)
        elif ax.get_legend() is not None:
            ax.get_legend().remove()
            
        return ax
    
    else:  # mode == "samples"
        # MODO "samples" - Contar el número de muestras afectadas por gen y tipo de variante
        
        # gene_variant_sample_counts: Dict[str (gene), Dict[str (variant_type), Set[str (sample_id)]]]
        gene_variant_sample_counts = {}

        # Iterar sobre cada gen único presente en los datos filtrados
        for gene_iter_val in data_filtered[gene_column].unique():
            if pd.isna(gene_iter_val) or gene_iter_val == "Unknown":
                continue
            
            gene_variant_sample_counts[gene_iter_val] = {}
            # Filtrar el DataFrame para obtener solo las filas correspondientes al gen actual
            gene_data_for_current_gene = data_filtered[data_filtered[gene_column] == gene_iter_val]

            if samples_as_columns:  # Formato ancho: las muestras son columnas TCGA-*
                for sample_col_name in sample_cols:
                    # Para cada fila (variante específica) del gen actual
                    for _, row_series in gene_data_for_current_gene.iterrows():
                        sample_genotype_value = str(row_series[sample_col_name]).strip().upper()
                        actual_variant_type = row_series[variant_column]
                        
                        if pd.isna(actual_variant_type) or actual_variant_type == "Unknown":
                            continue

                        is_mutation_present_in_sample = False
                        if '|' in sample_genotype_value:
                            alleles = sample_genotype_value.split('|')
                            if len(alleles) >= 2 and alleles[0] != alleles[1]: 
                                is_mutation_present_in_sample = True
                        elif '/' in sample_genotype_value: 
                            alleles = sample_genotype_value.split('/')
                            if len(alleles) >= 2 and alleles[0] != alleles[1]:
                                is_mutation_present_in_sample = True
                        elif sample_genotype_value not in ["", ".", "0", "0/0", "0|0"] and not pd.isna(row_series[sample_col_name]):
                            is_mutation_present_in_sample = True
                        
                        if is_mutation_present_in_sample:
                            if actual_variant_type not in gene_variant_sample_counts[gene_iter_val]:
                                gene_variant_sample_counts[gene_iter_val][actual_variant_type] = set()
                            gene_variant_sample_counts[gene_iter_val][actual_variant_type].add(sample_col_name)
            
            else:  # Formato largo: hay una columna 'sample_column'
                for _, row_series in gene_data_for_current_gene.iterrows():
                    actual_sample_id = row_series[sample_column]
                    actual_variant_type = row_series[variant_column]

                    if pd.isna(actual_variant_type) or actual_variant_type == "Unknown" or pd.isna(actual_sample_id):
                        continue
                    
                    if actual_variant_type not in gene_variant_sample_counts[gene_iter_val]:
                        gene_variant_sample_counts[gene_iter_val][actual_variant_type] = set()
                    gene_variant_sample_counts[gene_iter_val][actual_variant_type].add(actual_sample_id)

        plot_data_list = []
        for gene, variant_dict in gene_variant_sample_counts.items():
            row_for_df = {gene_column: gene}
            has_data_for_gene = False
            for variant_type, sample_set in variant_dict.items():
                if sample_set: 
                    row_for_df[variant_type] = len(sample_set)
                    has_data_for_gene = True
            if has_data_for_gene: 
                plot_data_list.append(row_for_df)
        
        if not plot_data_list:
            ax.text(0.5, 0.5, "No hay datos disponibles para analizar (modo samples)", 
                      ha='center', va='center', fontsize=12)
            ax.set_title(f"Top {count} Mutated Genes (Sample Prevalence)", fontsize=14, fontweight='bold')
            ax.axis('off')
            return ax

        gene_variant_counts_df = pd.DataFrame(plot_data_list).set_index(gene_column).fillna(0)
        
        gene_total_affected_samples = {}
        for gene, variant_dict in gene_variant_sample_counts.items():
            all_samples_for_gene = set()
            for _, sample_set in variant_dict.items(): 
                all_samples_for_gene.update(sample_set) 
            if all_samples_for_gene: 
                 gene_total_affected_samples[gene] = len(all_samples_for_gene)

        if not gene_total_affected_samples: 
            ax.text(0.5, 0.5, "No hay genes con muestras afectadas para mostrar (modo samples)", 
                      ha='center', va='center', fontsize=12)
            ax.set_title(f"Top {count} Mutated Genes (Sample Prevalence)", fontsize=14, fontweight='bold')
            ax.axis('off')
            return ax

        gene_totals_series = pd.Series(gene_total_affected_samples).sort_values(ascending=False)
        
        total_samples_in_dataset = len(sample_cols) if samples_as_columns else data_filtered[sample_column].nunique()

        top_genes_list = gene_totals_series.index[:count].tolist()
        
        df_top_plot = gene_variant_counts_df.loc[gene_variant_counts_df.index.isin(top_genes_list)]
        df_top_plot = df_top_plot.loc[:, (df_top_plot != 0).any(axis=0)] 

        valid_top_genes_for_plot = [g for g in top_genes_list if g in df_top_plot.index]
        ordered_genes_for_plot = sorted(valid_top_genes_for_plot, key=lambda g: gene_totals_series[g])
        
        if not ordered_genes_for_plot: 
            ax.text(0.5, 0.5, "No hay genes en top seleccionados para mostrar (modo samples)", 
                      ha='center', va='center', fontsize=12)
            ax.set_title(f"Top {count} Mutated Genes (Sample Prevalence)", fontsize=14, fontweight='bold')
            ax.axis('off')
            return ax

        df_plot_final = df_top_plot.loc[ordered_genes_for_plot]

        # NORMALIZACIÓN: Modificar df_plot_final para que la suma total de cada fila (gen) sea igual 
        # al número total de muestras únicas afectadas por ese gen
        normalized_df_plot = pd.DataFrame(index=df_plot_final.index, columns=df_plot_final.columns)
        
        for gene in df_plot_final.index:
            # Número total de muestras únicas afectadas por este gen (para el 100% de la barra)
            total_unique_samples = gene_totals_series[gene]
            
            # Valores actuales (sin normalizar) por tipo de variante para este gen
            current_values = df_plot_final.loc[gene]
            
            # Suma actual de los valores por tipo de variante
            current_sum = current_values.sum()
            
            if current_sum > 0:  # Evitar división por cero
                # Factor de normalización: cuánto representa cada unidad actual respecto al total de muestras únicas
                normalization_factor = total_unique_samples / current_sum
                
                # Calcular nuevos valores normalizados: los valores actuales * factor de normalización
                normalized_values = current_values * normalization_factor
                
                # Asignar los valores normalizados a la fila correspondiente al gen actual
                normalized_df_plot.loc[gene] = normalized_values
        
        # Reemplazar df_plot_final con la versión normalizada
        df_plot_final = normalized_df_plot.fillna(0)

        variant_types_in_plot = df_plot_final.columns.tolist()
        if color_map is None:
            cmap_instance = plt.colormaps.get_cmap('tab20')
            colors_for_plot = [cmap_instance(i % cmap_instance.N) for i in range(len(variant_types_in_plot))]
        else:
            cmap_instance = plt.colormaps.get_cmap('tab20')
            colors_for_plot = [color_map.get(vt, cmap_instance(i % cmap_instance.N)) 
                               for i, vt in enumerate(variant_types_in_plot)]
        
        df_plot_final.plot(kind='barh', stacked=True, ax=ax, color=colors_for_plot, width=0.65)
        
        for i, gene_name_in_plot in enumerate(df_plot_final.index):
            num_unique_samples_affected = gene_totals_series[gene_name_in_plot]
            percentage = (num_unique_samples_affected / total_samples_in_dataset) * 100 if total_samples_in_dataset > 0 else 0
            bar_length = df_plot_final.loc[gene_name_in_plot].sum()
            offset = 0.01 * ax.get_xlim()[1] if ax.get_xlim()[1] > 0 else 0.1 
            ax.text(bar_length + offset , i, f'{percentage:.1f}%', va='center', fontsize=10)
        
        title_text = f"Top {count} Mutated Genes (Sample Prevalence)"
        ax.set_title(title_text, fontsize=14, fontweight='bold') 
        ax.set_ylabel("")  # Eliminar la etiqueta del eje Y "Genes"
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False) # Ocultar la línea del eje izquierdo
        ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=True) # Sin ticks en eje Y, pero con etiquetas

        handles, labels = ax.get_legend_handles_labels()
        legend_elements = {label: handle for label, handle in zip(labels, handles)}
        # Asegurar que valid_handles y valid_labels se basen en las columnas de df_plot_final que realmente están en la leyenda
        valid_handles = [legend_elements[label] for label in df_plot_final.columns if label in legend_elements]
        valid_labels = [label for label in df_plot_final.columns if label in legend_elements]


        if valid_labels: 
            num_legend_items = len(valid_labels)
            ncol_legend = min(num_legend_items, 4) 
            # Ajuste más robusto para vertical_offset
            base_offset = -0.20 
            row_offset_factor = 0.06
            num_legend_rows = (num_legend_items + ncol_legend - 1) // ncol_legend
            vertical_offset = base_offset - (row_offset_factor * num_legend_rows)

            ax.legend(valid_handles, valid_labels, 
                    title="Variant Classification", 
                    loc="lower center", 
                    bbox_to_anchor=(0.5, vertical_offset), 
                    ncol=ncol_legend) 
        elif ax.get_legend() is not None: 
            ax.get_legend().remove()
            
        return ax