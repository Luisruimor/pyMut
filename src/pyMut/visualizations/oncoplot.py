"""
Funciones para crear oncoplots (también conocidos como waterfall plots).

Este módulo contiene las funciones necesarias para generar oncoplots, que son
visualizaciones heatmap que muestran patrones de mutación a través de muestras
y genes. Los oncoplots son fundamentales en la genómica del cáncer para
visualizar paisajes mutacionales.

La función principal es `create_oncoplot_plot()` que puede manejar:
- Detección automática de columnas de muestra (formato TCGA y .GT)
- Múltiples formatos de genotipo (A|G, A/G, etc.)
- Detección de Multi_Hit para muestras con múltiples tipos de mutación
- Esquemas de colores personalizables para clasificaciones de variantes
- Parámetros configurables para genes principales y muestras máximas

Funciones principales:
- is_mutated(): Determina si un genotipo representa una mutación
- detect_sample_columns(): Detecta automáticamente columnas de muestra en el DataFrame
- create_variant_color_mapping(): Crea mapeo de colores para tipos de variantes
- process_mutation_matrix(): Procesa datos de mutación en formato de matriz
- create_oncoplot_plot(): Función principal para crear el oncoplot
"""

import warnings
from typing import Dict, List, Optional, Set, Tuple, Union

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import seaborn as sns

from ..utils.constants import (
    GENE_COLUMN, VARIANT_CLASSIFICATION_COLUMN, REF_COLUMN, ALT_COLUMN,
    DEFAULT_ONCOPLOT_FIGSIZE, DEFAULT_ONCOPLOT_TOP_GENES, DEFAULT_ONCOPLOT_MAX_SAMPLES
)

def is_mutated(genotype: str, ref: str, alt: str) -> bool:
    """
    Determina si un genotipo representa una mutación comparando con REF/ALT.
    
    Maneja múltiples formatos de genotipo incluyendo:
    - Formato separado por pipe: "A|G", "A|A" 
    - Formato separado por barra: "A/G", "A/A"
    - Otros separadores comunes
    
    Args:
        genotype: Cadena de genotipo de la muestra
        ref: Alelo de referencia  
        alt: Alelo alternativo
        
    Returns:
        bool: True si el genotipo contiene el alelo alternativo
        
    Examples:
        >>> is_mutated("A|G", "A", "G")
        True
        >>> is_mutated("A|A", "A", "G") 
        False
        >>> is_mutated("G/G", "A", "G")
        True
    """
    if pd.isna(genotype) or pd.isna(ref) or pd.isna(alt):
        return False
        
    genotype = str(genotype).strip()
    ref = str(ref).strip()
    alt = str(alt).strip()
    
    # Manejo de casos especiales
    if genotype in ['', '.', 'nan', 'NaN']:
        return False
    
    # Separadores comunes para genotipos
    separators = ['|', '/', ':', ';', ',']
    alleles = [genotype]  # Por defecto, tratar como un solo alelo
    
    # Intentar dividir por separadores comunes
    for sep in separators:
        if sep in genotype:
            alleles = genotype.split(sep)
            break
    
    # Verificar si algún alelo coincide con el alternativo
    return alt in alleles

def detect_sample_columns(data: pd.DataFrame) -> List[str]:
    """
    Detecta automáticamente columnas de muestra en el DataFrame.
    
    Busca patrones comunes:
    - Formato TCGA: columnas que comienzan con "TCGA-"
    - Formato GT: columnas que terminan con ".GT"
    - Otros patrones de ID de muestra
    
    Args:
        data: DataFrame con datos de mutación
        
    Returns:
        List[str]: Lista de nombres de columnas identificadas como muestras
        
    Raises:
        ValueError: Si no se encuentran columnas de muestra
        
    Examples:
        >>> df = pd.DataFrame({'TCGA-AB-1234': [1], 'TCGA-CD-5678': [0]})
        >>> detect_sample_columns(df)
        ['TCGA-AB-1234', 'TCGA-CD-5678']
    """
    sample_columns = []
    
    for col in data.columns:
        # Formato TCGA (ej: TCGA-AB-1234)
        if col.startswith('TCGA-'):
            sample_columns.append(col)
        # Formato GT (ej: sample1.GT, sample2.GT)
        elif col.endswith('.GT'):
            sample_columns.append(col)
        # Otros patrones comunes de ID de muestra
        elif any(pattern in col.upper() for pattern in ['SAMPLE', 'PATIENT', 'SUBJECT']):
            if not any(keyword in col.lower() for keyword in ['barcode', 'id', 'name', 'type']):
                sample_columns.append(col)
    
    if not sample_columns:
        raise ValueError(
            "No se pudieron detectar columnas de muestra automáticamente. "
            "Asegúrese de que las columnas sigan formatos estándar como 'TCGA-*' o '*.GT', "
            "o especifique las columnas manualmente."
        )
    
    return sample_columns

def create_variant_color_mapping(variants: Set[str]) -> Dict[str, np.ndarray]:
    """
    Crea un mapeo de colores para clasificaciones de variantes.
    
    Genera colores únicos automáticamente para todos los tipos de variantes
    usando colormaps de matplotlib.
    
    Args:
        variants: Conjunto de clasificaciones de variantes únicas
        
    Returns:
        Dict[str, np.ndarray]: Mapeo de variante a color RGB
        
    Examples:
        >>> variants = {'Missense_Mutation', 'Nonsense_Mutation'}
        >>> colors = create_variant_color_mapping(variants)
        >>> 'Missense_Mutation' in colors
        True
    """
    color_mapping = {}
    
    # Ordenar variantes para consistencia
    sorted_variants = sorted(variants)
    
    # Usar colormap Set3 para generar colores únicos
    n_variants = len(sorted_variants)
    if n_variants > 0:
        # Set3 es bueno para hasta 12 colores, después usamos tab20
        if n_variants <= 12:
            colors = plt.cm.Set3(np.linspace(0, 1, max(n_variants, 3)))
        else:
            colors = plt.cm.tab20(np.linspace(0, 1, n_variants))
        
        for i, variant in enumerate(sorted_variants):
            if variant == 'None':
                # Color gris claro para 'None'
                color_mapping[variant] = np.array([0.9, 0.9, 0.9])
            else:
                color_mapping[variant] = colors[i][:3]  # Solo RGB, sin alpha
    
    return color_mapping

def process_mutation_matrix(data: pd.DataFrame,
                           gene_column: str = GENE_COLUMN,
                           variant_column: str = VARIANT_CLASSIFICATION_COLUMN,
                           ref_column: str = REF_COLUMN,
                           alt_column: str = ALT_COLUMN,
                           sample_columns: Optional[List[str]] = None) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """
    Procesa datos de mutación en formato de matriz para el oncoplot.
    
    Convierte datos de mutación de formato largo a formato de matriz donde:
    - Filas = genes
    - Columnas = muestras  
    - Valores = tipos de mutación o 'None' si no hay mutación
    
    Maneja detección de Multi_Hit cuando una muestra tiene múltiples
    tipos de mutación en el mismo gen.
    
    Args:
        data: DataFrame con datos de mutación
        gene_column: Nombre de la columna de genes
        variant_column: Nombre de la columna de clasificación de variantes
        ref_column: Nombre de la columna de alelo de referencia  
        alt_column: Nombre de la columna de alelo alternativo
        sample_columns: Lista de columnas de muestra (detecta automáticamente si es None)
        
    Returns:
        Tuple[pd.DataFrame, Dict[str, int]]: 
            - DataFrame de matriz de mutación con genes como índice y muestras como columnas
            - Diccionario con recuentos de mutación por gen
            
    Raises:
        ValueError: Si faltan columnas requeridas o no se encuentran muestras
        
    Examples:
        >>> df = pd.DataFrame({
        ...     'Hugo_Symbol': ['TP53', 'TP53'],
        ...     'Variant_Classification': ['Missense_Mutation', 'Nonsense_Mutation'],
        ...     'REF': ['A', 'C'],
        ...     'ALT': ['G', 'T'],
        ...     'TCGA-AB-1234': ['A|G', 'C|C'],
        ...     'TCGA-CD-5678': ['A|A', 'C|T']
        ... })
        >>> matrix, counts = process_mutation_matrix(df)
        >>> matrix.shape
        (1, 2)
    """
    # Validar columnas requeridas
    required_columns = [gene_column, variant_column, ref_column, alt_column]
    missing_columns = [col for col in required_columns if col not in data.columns]
    if missing_columns:
        raise ValueError(f"Faltan columnas requeridas: {missing_columns}")
    
    # Detectar columnas de muestra automáticamente si no se proporcionan
    if sample_columns is None:
        sample_columns = detect_sample_columns(data)
    
    # Validar que las columnas de muestra existen
    missing_samples = [col for col in sample_columns if col not in data.columns]
    if missing_samples:
        raise ValueError(f"Columnas de muestra no encontradas: {missing_samples}")
    
    # Preparar datos
    working_data = data[[gene_column, variant_column, ref_column, alt_column] + sample_columns].copy()
    
    # Obtener genes únicos
    genes = working_data[gene_column].unique()
    
    # Inicializar matriz de mutación
    mutation_matrix = pd.DataFrame(index=genes, columns=sample_columns, dtype=object)
    mutation_matrix = mutation_matrix.fillna('None')
    
    # Procesar cada fila de datos
    for _, row in working_data.iterrows():
        gene = row[gene_column]
        variant_type = row[variant_column]
        ref_allele = row[ref_column] 
        alt_allele = row[alt_column]
        
        # Verificar mutaciones en cada muestra
        for sample in sample_columns:
            genotype = row[sample]
            
            if is_mutated(genotype, ref_allele, alt_allele):
                current_value = mutation_matrix.loc[gene, sample]
                
                if current_value == 'None':
                    # Primera mutación en este gen/muestra
                    mutation_matrix.loc[gene, sample] = variant_type
                elif current_value != variant_type:
                    # Múltiples tipos de mutación = Multi_Hit
                    mutation_matrix.loc[gene, sample] = 'Multi_Hit'
    
    # Calcular recuentos de mutación por gen
    mutation_counts = {}
    for gene in genes:
        gene_row = mutation_matrix.loc[gene]
        mutation_counts[gene] = (gene_row != 'None').sum()
    
    return mutation_matrix, mutation_counts

def create_oncoplot_plot(data: pd.DataFrame,
                         gene_column: str = GENE_COLUMN,
                         variant_column: str = VARIANT_CLASSIFICATION_COLUMN,
                         ref_column: str = REF_COLUMN,
                         alt_column: str = ALT_COLUMN,
                         top_genes_count: int = 30,
                         max_samples: int = 180,
                         figsize: Optional[Tuple[int, int]] = None,
                         title: str = "Oncoplot",
                         ax: Optional[plt.Axes] = None) -> plt.Figure:
    """
    Crea un oncoplot (heatmap de mutaciones) mostrando patrones de mutación.
    
    Genera una visualización heatmap donde:
    - Eje Y: Genes (ordenados por frecuencia de mutación)
    - Eje X: Muestras (ordenadas por carga mutacional para efecto cascada)  
    - Colores: Tipos de clasificación de variantes
    
    Args:
        data: DataFrame con datos de mutación
        gene_column: Nombre de la columna de genes (por defecto 'Hugo_Symbol')
        variant_column: Nombre de la columna de clasificación de variantes
        ref_column: Nombre de la columna de alelo de referencia
        alt_column: Nombre de la columna de alelo alternativo  
        top_genes_count: Número de genes más mutados a mostrar (por defecto 30)
        max_samples: Número máximo de muestras a mostrar (por defecto 180)
        figsize: Tamaño de la figura (ancho, alto) en pulgadas
        title: Título para el gráfico
        ax: Ejes matplotlib existentes (opcional)
        
    Returns:
        plt.Figure: Objeto Figure de matplotlib con el oncoplot
        
    Raises:
        ValueError: Si faltan columnas requeridas o no hay datos para visualizar
        
    Examples:
        >>> fig = create_oncoplot_plot(mutation_data, top_genes_count=20)
        >>> fig.savefig('oncoplot.png', dpi=300, bbox_inches='tight')
    """
    try:
        # Procesar datos en matriz de mutación
        mutation_matrix, mutation_counts = process_mutation_matrix(
            data, gene_column, variant_column, ref_column, alt_column
        )
        
        if mutation_matrix.empty:
            raise ValueError("No hay datos de mutación para visualizar")
        
        # Seleccionar top genes más mutados
        top_genes = sorted(mutation_counts.items(), key=lambda x: x[1], reverse=True)
        top_genes = [gene for gene, count in top_genes[:top_genes_count] if count > 0]
        
        if not top_genes:
            raise ValueError("No se encontraron genes con mutaciones")
        
        # Filtrar matriz a top genes
        plot_matrix = mutation_matrix.loc[top_genes].copy()
        
        # ALGORITMO DE ORDENAMIENTO MEJORADO PARA EFECTO CASCADA
        # 1. Calcular TMB (Tumor Mutation Burden) por muestra
        sample_tmb = (plot_matrix != 'None').sum(axis=0)
        
        # 2. Crear score combinado para ordenamiento cascada
        # Priorizar muestras con alta TMB y mutaciones en genes importantes
        sample_scores = {}
        
        for sample in plot_matrix.columns:
            # Score base: TMB (carga mutacional total)
            tmb_score = sample_tmb[sample] * 100
            
            # Score adicional: mutaciones en genes top (más peso a genes más frecuentes)
            priority_score = 0
            for i, gene in enumerate(top_genes[:10]):  # Usar top 10 genes
                if plot_matrix.loc[gene, sample] != 'None':
                    # Dar más peso a genes más importantes (posición en ranking)
                    gene_weight = 10 - i
                    priority_score += gene_weight
            
            # Score final = TMB + bonus por mutaciones en genes importantes
            sample_scores[sample] = tmb_score + priority_score
        
        # 3. Ordenar muestras por score (descendente = efecto cascada)
        sorted_samples = sorted(sample_scores.keys(), 
                              key=lambda x: sample_scores[x], 
                              reverse=True)
        
        # 4. Limitar número de muestras si es necesario
        if len(sorted_samples) > max_samples:
            sorted_samples = sorted_samples[:max_samples]
            
        # 5. Reordenar matriz con muestras ordenadas
        plot_matrix = plot_matrix[sorted_samples]
        
        # Crear mapeo de colores
        unique_variants = set()
        for row in plot_matrix.values.flatten():
            unique_variants.add(row)
        
        color_mapping = create_variant_color_mapping(unique_variants)
        
        # Convertir valores categóricos a numéricos para heatmap
        unique_values = sorted(unique_variants)
        value_to_num = {value: i for i, value in enumerate(unique_values)}
        
        numeric_matrix = plot_matrix.applymap(lambda x: value_to_num[x])
        
        # Configurar figura
        if figsize is None:
            figsize = DEFAULT_ONCOPLOT_FIGSIZE
            
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()
        
        # Crear colormap personalizado
        colors = [color_mapping[value] for value in unique_values]
        custom_cmap = mcolors.ListedColormap(colors)
        
        # Crear heatmap principal
        im = ax.imshow(
            numeric_matrix.values,
            cmap=custom_cmap,
            aspect='auto',
            interpolation='nearest'
        )
        
        # Configurar ejes
        ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
        ax.set_xlabel(f'Muestras (n={len(sorted_samples)})', fontsize=12)
        ax.set_ylabel('Genes', fontsize=12)
        
        # Configurar ticks y labels
        ax.set_xticks(range(len(sorted_samples)))
        ax.set_yticks(range(len(top_genes)))
        ax.set_yticklabels(top_genes, fontsize=10)
        
        # Rotar etiquetas de muestras si hay muchas
        if len(sorted_samples) > 20:
            ax.set_xticklabels(sorted_samples, rotation=90, fontsize=8)
        else:
            ax.set_xticklabels(sorted_samples, rotation=45, fontsize=9)
        
        # Añadir grid sutil
        ax.set_xticks(np.arange(-0.5, len(sorted_samples), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(top_genes), 1), minor=True)
        ax.grid(which='minor', color='white', linestyle='-', linewidth=0.5)
        
        # Añadir leyenda
        legend_elements = []
        for variant in sorted(unique_variants):
            if variant != 'None':  # No mostrar 'None' en la leyenda
                color = color_mapping[variant]
                label = variant.replace('_', ' ')
                legend_elements.append(
                    plt.Rectangle((0, 0), 1, 1, facecolor=color, label=label)
                )
        
        if legend_elements:
            ax.legend(handles=legend_elements, 
                     bbox_to_anchor=(1.05, 1), 
                     loc='upper left',
                     fontsize=10,
                     title='Variant Classification',
                     title_fontsize=11,
                     frameon=True)
        
        # Ajustar layout
        plt.tight_layout()
        
        return fig
        
    except Exception as e:
        raise ValueError(f"Error al crear oncoplot: {str(e)}") 