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

# Importar la función de top mutated genes del módulo summary
from .summary import create_top_mutated_genes_plot

def is_mutated(genotype: str, ref: str, alt: str) -> bool:
    """
    Determina si un genotipo representa una mutación comparando con REF/ALT.
    
    Soporta múltiples formatos de genotipo:
    - Formato pipe: "A|G", "C|T" 
    - Formato slash: "A/G", "C/T"
    - Otros separadores: "A:G", "A;G"
    - Casos especiales: maneja valores faltantes y formatos no estándar
    
    Args:
        genotype: Valor del genotipo de la muestra
        ref: Alelo de referencia  
        alt: Alelo alternativo
        
    Returns:
        bool: True si representa una mutación, False en caso contrario
        
    Examples:
        >>> is_mutated("A|G", "A", "G")  # True
        >>> is_mutated("A|A", "A", "G")  # False  
        >>> is_mutated("G|G", "A", "G")  # True (homocigoto alternativo)
        >>> is_mutated("./.", "A", "G")  # False
    """
    # Verificar valores nulos o vacíos
    if pd.isna(genotype) or pd.isna(ref) or pd.isna(alt):
        return False
        
    genotype_str = str(genotype).strip()
    ref_str = str(ref).strip()
    alt_str = str(alt).strip()
    
    # Casos que claramente NO son mutaciones
    no_mutation_values = {"", ".", "0", "0/0", "0|0", "./.", ".|.", "NA", "NaN"}
    if genotype_str.upper() in no_mutation_values:
        return False
    
    # Separadores comunes para genotipos
    separators = ['|', '/', ':', ';', ',']
    alleles = [genotype_str]  # Por defecto, tratar como un solo alelo
    
    # Intentar dividir por separadores comunes
    for sep in separators:
        if sep in genotype_str:
            alleles = genotype_str.split(sep)
            break
    
    # Verificar si algún alelo coincide con el alternativo
    # Esto detecta tanto heterocigotos (A|G) como homocigotos alternativos (G|G)
    return alt_str in alleles

def detect_sample_columns(data: pd.DataFrame) -> List[str]:
    """
    Detecta automáticamente las columnas que representan muestras en el DataFrame.
    
    Busca patrones comunes de nomenclatura de muestras:
    - Formato TCGA: columnas que empiezan con "TCGA-"
    - Formato GT: columnas que terminan con ".GT"
    - Otros patrones de muestras según convenciones estándar
    
    Args:
        data: DataFrame con datos de mutación
        
    Returns:
        List[str]: Lista de nombres de columnas identificadas como muestras
        
    Raises:
        ValueError: Si no se detectan columnas de muestras
        
    Examples:
        >>> columns = ["Hugo_Symbol", "TCGA-AB-2988", "TCGA-AB-2869", "Variant_Classification"]
        >>> detect_sample_columns(pd.DataFrame(columns=columns))
        ['TCGA-AB-2988', 'TCGA-AB-2869']
    """
    sample_columns = []
    
    # Detectar formato TCGA (más común)
    tcga_cols = [col for col in data.columns if str(col).startswith("TCGA-")]
    sample_columns.extend(tcga_cols)
    
    # Detectar formato .GT (genotype)
    gt_cols = [col for col in data.columns if str(col).endswith(".GT")]
    sample_columns.extend(gt_cols)
    
    # Eliminar duplicados manteniendo orden
    sample_columns = list(dict.fromkeys(sample_columns))
    
    if not sample_columns:
        raise ValueError(
            "No se pudieron detectar columnas de muestra automáticamente. "
            "Asegúrese de que las columnas sigan formatos estándar como 'TCGA-*' o '*.GT', "
            "o especifique las columnas manualmente."
        )
        
    return sample_columns

def create_variant_color_mapping(variants: Set[str]) -> Dict[str, np.ndarray]:
    """
    Crea un mapeo de colores consistente para tipos de variantes.
    
    Asigna colores específicos a tipos de variantes conocidos y colores 
    automáticos para variantes no reconocidas.
    
    Args:
        variants: Conjunto de tipos de variantes únicos
        
    Returns:
        Dict[str, np.ndarray]: Diccionario mapeando variante -> color RGB
        
    Examples:
        >>> variants = {"Missense_Mutation", "Nonsense_Mutation", "None"}
        >>> mapping = create_variant_color_mapping(variants)
        >>> len(mapping) == 3
        True
    """
    # Definir colores específicos para tipos de variantes conocidos
    predefined_colors = {
        'Missense_Mutation': np.array([0.0, 0.8, 0.0]),        # Verde
        'Nonsense_Mutation': np.array([1.0, 0.0, 0.0]),        # Rojo
        'Frame_Shift_Del': np.array([0.8, 0.0, 0.8]),          # Magenta
        'Frame_Shift_Ins': np.array([0.6, 0.0, 0.6]),          # Magenta oscuro
        'In_Frame_Del': np.array([0.0, 0.8, 0.8]),             # Cian
        'In_Frame_Ins': np.array([0.0, 0.6, 0.6]),             # Cian oscuro
        'Splice_Site': np.array([1.0, 0.5, 0.0]),              # Naranja
        'Translation_Start_Site': np.array([0.8, 0.8, 0.0]),   # Amarillo
        'Nonstop_Mutation': np.array([0.5, 0.0, 1.0]),         # Azul violeta
        'Silent': np.array([0.7, 0.7, 0.7]),                   # Gris
        'Multi_Hit': np.array([0.0, 0.0, 0.0]),                # Negro
        'None': np.array([0.95, 0.95, 0.95])                   # Gris muy claro
    }
    
    color_mapping = {}
    
    # Usar colores predefinidos cuando estén disponibles
    for variant in variants:
        if variant in predefined_colors:
            color_mapping[variant] = predefined_colors[variant]
    
    # Asignar colores automáticos para variantes no reconocidas
    unassigned_variants = [v for v in variants if v not in color_mapping]
    
    if unassigned_variants:
        # Generar colores automáticos usando una paleta
        n_colors = len(unassigned_variants)
        auto_colors = plt.cm.tab20(np.linspace(0, 1, n_colors))
        
        for i, variant in enumerate(unassigned_variants):
            color_mapping[variant] = auto_colors[i][:3]  # Solo RGB, sin alpha
    
    return color_mapping

def process_mutation_matrix(data: pd.DataFrame,
                           gene_column: str = GENE_COLUMN,
                           variant_column: str = VARIANT_CLASSIFICATION_COLUMN,
                           ref_column: str = REF_COLUMN,
                           alt_column: str = ALT_COLUMN,
                           sample_columns: Optional[List[str]] = None) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """
    Procesa datos de mutación para crear una matriz genes x muestras.
    
    Transforma datos de mutación en formato largo o ancho a una matriz optimizada
    para visualización oncoplot. Detecta automáticamente múltiples mutaciones 
    (Multi_Hit) y maneja diferentes formatos de genotipo.
    
    Args:
        data: DataFrame con datos de mutación
        gene_column: Nombre de la columna de genes  
        variant_column: Nombre de la columna de clasificación de variantes
        ref_column: Nombre de la columna de alelo de referencia
        alt_column: Nombre de la columna de alelo alternativo
        sample_columns: Lista opcional de columnas de muestras (autodetecta si None)
        
    Returns:
        Tuple[pd.DataFrame, Dict[str, int]]: 
            - Matriz de mutaciones (genes x muestras)
            - Diccionario con conteos de mutaciones por gen
            
    Raises:
        ValueError: Si faltan columnas requeridas o no hay datos válidos
        
    Examples:
        >>> matrix, counts = process_mutation_matrix(mutation_data)
        >>> matrix.shape  # (genes, samples)
        (100, 50)
        >>> "Multi_Hit" in matrix.values
        True
    """
    # Validar columnas requeridas
    required_columns = [gene_column, variant_column]
    missing_columns = [col for col in required_columns if col not in data.columns]
    
    if missing_columns:
        raise ValueError(f"Faltan columnas requeridas: {missing_columns}")
    
    # Detectar columnas de muestras si no se proporcionaron
    if sample_columns is None:
        try:
            sample_columns = detect_sample_columns(data)
        except ValueError as e:
            # Si no se pueden detectar automáticamente, mostrar warning pero continuar
            print(f"Advertencia: {e}")
            print(f"Columnas disponibles: {list(data.columns)}")
            # Intentar identificar columnas que podrían ser muestras
            excluded_cols = {gene_column, variant_column, ref_column, alt_column, 'Tumor_Sample_Barcode'}
            potential_samples = [col for col in data.columns if col not in excluded_cols]
            if potential_samples:
                print(f"Usando columnas potenciales como muestras: {potential_samples}")
                sample_columns = potential_samples
            else:
                raise ValueError("No se detectaron columnas de muestras válidas")
        
    if not sample_columns:
        raise ValueError("No se detectaron columnas de muestras válidas")
    
    print(f"Procesando {len(sample_columns)} muestras y datos de {len(data)} variantes...")
    
    # Crear matriz base: genes x muestras, inicializada con 'None'
    unique_genes = data[gene_column].dropna().unique()
    unique_genes = [gene for gene in unique_genes if str(gene) != 'nan' and str(gene).strip() != '']
    
    if not unique_genes:
        raise ValueError("No se encontraron genes válidos en los datos")
    
    # Crear DataFrame con 'None' como valor por defecto
    mutation_matrix = pd.DataFrame(
        'None',
        index=unique_genes,
        columns=sample_columns
    )
    
    # Verificar que existen las columnas REF y ALT si están especificadas
    has_ref_alt = (ref_column in data.columns and alt_column in data.columns)
    
    # Procesar cada fila de datos de mutación
    mutation_counts = {gene: 0 for gene in unique_genes}
    
    for _, row in data.iterrows():
        gene = row[gene_column]
        variant_type = row[variant_column]
        
        # Saltar filas con valores faltantes
        if pd.isna(gene) or pd.isna(variant_type):
            continue
            
        gene_str = str(gene).strip()
        variant_str = str(variant_type).strip()
        
        if gene_str == '' or variant_str == '' or gene_str not in unique_genes:
            continue
        
        # Obtener alelos de referencia y alternativo si están disponibles
        ref_allele = row[ref_column] if has_ref_alt else None
        alt_allele = row[alt_column] if has_ref_alt else None
        
        # Procesar cada muestra para esta variante
        for sample_col in sample_columns:
            if sample_col not in data.columns:
                continue
                
            genotype = row[sample_col]
            
            # Determinar si hay mutación en esta muestra
            has_mutation = False
            
            if has_ref_alt and not pd.isna(ref_allele) and not pd.isna(alt_allele):
                # Usar información REF/ALT para determinar mutación
                has_mutation = is_mutated(genotype, ref_allele, alt_allele)
            else:
                # Usar heurística basada en genotype solamente
                # Para casos sin REF/ALT, considerar que valores no-vacíos/no-missing son mutaciones
                if not pd.isna(genotype):
                    genotype_str = str(genotype).strip()
                    no_mutation_values = {"", ".", "0", "0/0", "0|0", "./.", ".|.", "NA", "NaN", "None"}
                    has_mutation = genotype_str not in no_mutation_values
            
            if has_mutation:
                current_value = mutation_matrix.loc[gene_str, sample_col]
                
                if current_value == 'None':
                    # Primera mutación en esta celda
                    mutation_matrix.loc[gene_str, sample_col] = variant_str
                    mutation_counts[gene_str] += 1
                elif current_value != variant_str:
                    # Multiple mutaciones diferentes -> Multi_Hit
                    mutation_matrix.loc[gene_str, sample_col] = 'Multi_Hit'
                # Si es la misma variante, no hacer nada (evitar duplicados)
    
    print(f"Matriz procesada: {mutation_matrix.shape[0]} genes x {mutation_matrix.shape[1]} muestras")
    print(f"Genes con mutaciones: {sum(1 for count in mutation_counts.values() if count > 0)}")
    
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
    Crea un oncoplot completo con panel TMB superior y panel lateral de genes.
    
    Genera una visualización comprehensiva tipo oncoplot que incluye:
    - Panel superior: TMB (Tumor Mutation Burden) por muestra
    - Panel principal: Heatmap de mutaciones (genes x muestras)
    - Panel lateral derecho: Top genes mutados en modo "samples"
    - Leyenda: Tipos de clasificación de variantes
    - Algoritmo de ordenamiento waterfall/cascada para efecto visual
    
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
        ax: Ejes matplotlib existentes (opcional, para modo simple sin paneles)
        
    Returns:
        plt.Figure: Objeto Figure de matplotlib con el oncoplot completo
        
    Raises:
        ValueError: Si faltan columnas requeridas o no hay datos para visualizar
        
    Examples:
        >>> fig = create_oncoplot_plot(mutation_data, top_genes_count=20)
        >>> fig.savefig('oncoplot_completo.png', dpi=300, bbox_inches='tight')
    """
    try:
        # Procesar datos en matriz de mutación
        mutation_matrix, mutation_counts = process_mutation_matrix(
            data, gene_column, variant_column, ref_column, alt_column
        )
        
        if mutation_matrix.empty:
            raise ValueError("No hay datos de mutación para visualizar")
        
        # CALCULAR TMB USANDO TODOS LOS GENES (ANTES DEL FILTRADO)
        # Esto debe coincidir con el cálculo de variants_per_sample
        all_genes_tmb = (mutation_matrix != 'None').sum(axis=0)
        median_tmb_all = np.median(all_genes_tmb)
        
        # Seleccionar top genes más mutados
        top_genes = sorted(mutation_counts.items(), key=lambda x: x[1], reverse=True)
        top_genes = [gene for gene, count in top_genes[:top_genes_count] if count > 0]
        
        if not top_genes:
            raise ValueError("No se encontraron genes con mutaciones")
        
        # Filtrar matriz a top genes
        plot_matrix = mutation_matrix.loc[top_genes].copy()
        
        # ALGORITMO CLÁSICO DE WATERFALL/CASCADA PARA ONCOPLOTS
        # Inspirado en maftools y cBioPortal - crear efecto cascada real
        
        # 1. Calcular score de cascada para cada muestra
        # El efecto cascada clásico usa TMB + patrones específicos de mutación
        def waterfall_score(sample):
            # TMB total (factor principal)
            tmb = all_genes_tmb[sample]
            
            # Score adicional basado en mutaciones en top genes (con peso descendente)
            mutation_score = 0
            for i, gene in enumerate(top_genes):
                if plot_matrix.loc[gene, sample] != 'None':
                    # Peso descendente: primer gen más peso, último gen menos peso
                    weight = len(top_genes) - i
                    mutation_score += weight
            
            # Score final: TMB * 100 + patrón de mutaciones
            # Esto asegura que muestras con TMB similar se ordenen por patrones
            return tmb * 100 + mutation_score
        
        # 2. Ordenar muestras por score cascada (descendente)
        sorted_samples = sorted(plot_matrix.columns, 
                              key=waterfall_score, 
                              reverse=True)
        
        # 3. Limitar número de muestras si es necesario
        if len(sorted_samples) > max_samples:
            sorted_samples = sorted_samples[:max_samples]
            
        # 4. Reordenar matriz con muestras ordenadas
        plot_matrix = plot_matrix[sorted_samples]
        
        # 5. Obtener TMB para visualización
        # TMB para heatmap: ordenado según el efecto cascada  
        sorted_tmb_heatmap = [all_genes_tmb[sample] for sample in sorted_samples]
        median_tmb = np.median(sorted_tmb_heatmap)
        
        # TMB para panel superior: usar orden natural/alfabético de las muestras seleccionadas
        # Esto hace que el TMB no esté ordenado por carga mutacional
        tmb_display_order = sorted(sorted_samples)  # Orden alfabético
        sorted_tmb_display = [all_genes_tmb[sample] for sample in tmb_display_order]
        
        # Crear mapeo de colores
        unique_variants = set()
        for row in plot_matrix.values.flatten():
            unique_variants.add(row)
        
        color_mapping = create_variant_color_mapping(unique_variants)
        
        # Convertir valores categóricos a numéricos para heatmap
        unique_values = sorted(unique_variants)
        value_to_num = {value: i for i, value in enumerate(unique_values)}
        
        numeric_matrix = plot_matrix.applymap(lambda x: value_to_num[x])
        
        # Configurar figura y subplots
        if figsize is None:
            figsize = DEFAULT_ONCOPLOT_FIGSIZE
            
        if ax is None:
            # Crear figura completa con grid: TMB arriba, oncoplot principal, genes lateral derecho
            fig = plt.figure(figsize=figsize)
            
            # Definir grid layout: 3 filas x 2 columnas
            # height_ratios: [2.5, 10, 1.5] = TMB más alto arriba, heatmap grande, leyenda con espacio abajo
            # width_ratios: [10, 1] = heatmap ancho, panel genes lateral estrecho
            gs = fig.add_gridspec(3, 2, 
                                height_ratios=[2.5, 10, 1.5], 
                                width_ratios=[10, 1],
                                hspace=0.05, 
                                wspace=0.02)
            
            # Crear subplots según el grid
            ax_tmb = fig.add_subplot(gs[0, 0])      # Panel TMB (arriba izquierda)
            ax_main = fig.add_subplot(gs[1, 0])     # Panel principal heatmap (centro izquierda)
            ax_genes = fig.add_subplot(gs[1, 1])    # Panel genes (centro derecha)
            ax_legend = fig.add_subplot(gs[2, :])   # Panel leyenda (abajo, span completo)
            
            # === PANEL TMB (SUPERIOR) ===
            # Calcular TMB usando la misma lógica que variants_per_sample_plot
            # Usar las mismas columnas de muestra que se detectaron para el heatmap
            detected_sample_columns = mutation_matrix.columns.tolist()
            
            # Usar la misma lógica que variants_per_sample_plot
            variant_counts_tmb = {}
            
            # Agrupar por Variant_Classification (igual que en summary.py)
            for variant_type in data[variant_column].unique():
                variant_subset = data[data[variant_column] == variant_type]
                
                for sample_col in detected_sample_columns:
                    # Contar mutaciones usando la misma lógica que summary.py
                    sample_variants = variant_subset[sample_col].apply(
                        lambda x: 1 if '|' in str(x) and str(x).split('|')[0] != str(x).split('|')[1] else 0
                    ).sum()
                    
                    if sample_col not in variant_counts_tmb:
                        variant_counts_tmb[sample_col] = {}
                    
                    if sample_variants > 0:
                        variant_counts_tmb[sample_col][variant_type] = sample_variants
            
            # Crear DataFrame para TMB (igual que en summary.py)
            samples_df_tmb = []
            for sample, variants in variant_counts_tmb.items():
                for var_type, count in variants.items():
                    samples_df_tmb.append({'Sample': sample, 'Variant_Classification': var_type, 'Count': count})
            
            if samples_df_tmb:
                processed_df_tmb = pd.DataFrame(samples_df_tmb)
                tmb_df = processed_df_tmb.pivot(index='Sample', columns='Variant_Classification', values='Count').fillna(0)
                
                # Filtrar a las muestras seleccionadas y ordenarlas como en el heatmap
                tmb_df = tmb_df.loc[tmb_df.index.intersection(sorted_samples)]
                tmb_df = tmb_df.reindex(sorted_samples, fill_value=0)
                
                # Preparar colores para TMB (mismo mapeo que el heatmap)
                tmb_colors = [color_mapping.get(vt, [0.7, 0.7, 0.7]) for vt in tmb_df.columns]
                
                # Crear barras apiladas para TMB
                tmb_df.plot(kind='bar', stacked=True, ax=ax_tmb, color=tmb_colors, width=0.8, legend=False)
                
                # Configurar panel TMB
                ax_tmb.set_ylabel('TMB', fontsize=10)
                ax_tmb.set_xlim(-0.5, len(sorted_samples) - 0.5)
                
                # Los límites Y se ajustan automáticamente según los datos
                max_tmb = tmb_df.sum(axis=1).max() if not tmb_df.empty else 10
                ax_tmb.set_ylim(0, max_tmb * 1.1)
                
                ax_tmb.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
                ax_tmb.spines['top'].set_visible(False)
                ax_tmb.spines['right'].set_visible(False)
                ax_tmb.spines['bottom'].set_visible(False)
                ax_tmb.grid(axis='y', alpha=0.3)
            else:
                # Si no hay datos TMB, mostrar mensaje
                ax_tmb.text(0.5, 0.5, 'TMB no disponible', ha='center', va='center', fontsize=10)
                ax_tmb.set_xlim(0, 1)
                ax_tmb.set_ylim(0, 1)
                ax_tmb.axis('off')
            
            # === PANEL GENES LATERAL (DERECHO) ===
            # Crear panel lateral personalizado con barras apiladas y colores sincronizados
            try:
                # Calcular conteos de variantes por gen usando la matriz procesada
                gene_variant_counts = {}
                
                # Para cada gen, contar cuántas muestras tienen cada tipo de variante
                for gene in top_genes:
                    gene_variant_counts[gene] = {}
                    gene_row = plot_matrix.loc[gene]
                    
                    # Contar cada tipo de variante (excluyendo 'None')
                    for variant_type in gene_row.value_counts().index:
                        if variant_type != 'None':
                            count = gene_row.value_counts()[variant_type]
                            gene_variant_counts[gene][variant_type] = count
                
                # Crear DataFrame para barras apiladas
                variant_types = set()
                for gene_variants in gene_variant_counts.values():
                    variant_types.update(gene_variants.keys())
                
                variant_types = sorted(list(variant_types))
                
                # Construir DataFrame con genes como índice y tipos de variantes como columnas
                stacked_data = []
                for gene in top_genes:
                    row_data = {'gene': gene}
                    for variant_type in variant_types:
                        row_data[variant_type] = gene_variant_counts[gene].get(variant_type, 0)
                    stacked_data.append(row_data)
                
                df_stacked = pd.DataFrame(stacked_data).set_index('gene')
                
                # INVERTIR el orden para que coincida con el heatmap (arriba hacia abajo)
                # El heatmap muestra genes en el orden top_genes (de arriba hacia abajo)
                # Las barras horizontales se plotean de abajo hacia arriba por defecto
                # Por lo tanto, invertimos el DataFrame para que coincidan visualmente
                df_stacked = df_stacked.iloc[::-1]  # Invertir el orden de las filas
                
                # Calcular totales y porcentajes para cada gen
                gene_totals = df_stacked.sum(axis=1)
                total_samples = len(plot_matrix.columns)
                
                # Normalizar para mostrar prevalencia en muestras (como porcentaje del total de muestras)
                df_percentage = df_stacked.copy()
                for gene in df_percentage.index:
                    total_affected = gene_totals[gene]
                    percentage = (total_affected / total_samples) * 100
                    # Escalar cada variante proporcionalmente para que la suma sea el porcentaje total
                    if total_affected > 0:
                        scaling_factor = percentage / total_affected
                        df_percentage.loc[gene] = df_percentage.loc[gene] * scaling_factor
                
                # Crear mapeo de colores sincronizado con el heatmap principal
                stacked_colors = []
                for variant_type in variant_types:
                    if variant_type in color_mapping:
                        stacked_colors.append(color_mapping[variant_type])
                    else:
                        # Color por defecto si no está en el mapeo
                        stacked_colors.append([0.7, 0.7, 0.7])
                
                # Crear barras apiladas horizontales
                df_percentage.plot(kind='barh', stacked=True, ax=ax_genes, 
                                 color=stacked_colors, width=0.65)
                
                # Configurar panel genes - SIN etiquetas de genes y SIN etiqueta en X
                ax_genes.set_xlabel('')  # ELIMINAR la etiqueta "Muestras (%)" del eje X  
                ax_genes.set_ylabel('')  # Sin etiqueta Y
                
                # Calcular límite X basado en los datos
                max_percentage = df_percentage.sum(axis=1).max() if not df_percentage.empty else 50
                ax_genes.set_xlim(0, max_percentage * 1.1)
                ax_genes.set_ylim(-0.5, len(top_genes) - 0.5)
                
                # ELIMINAR etiquetas de genes del eje Y  
                ax_genes.set_yticklabels([''] * len(top_genes))  # Etiquetas vacías
                
                # CORREGIR PORCENTAJES: usar exactamente la misma lógica que summary.py
                # Detectar columnas de muestra igual que en summary.py
                sample_cols_for_percentage = [col for col in data.columns if str(col).startswith("TCGA-")]
                total_samples_in_dataset = len(sample_cols_for_percentage)  # Igual que en summary.py
                
                for i, gene in enumerate(df_percentage.index):
                    # Revertir el orden invertido para obtener el gen correcto
                    gene_in_original_order = top_genes[len(top_genes) - 1 - i]  
                    
                    # USAR LA MISMA LÓGICA QUE SUMMARY.PY:
                    # Contar muestras únicas afectadas por este gen usando los datos originales
                    gene_data_filtered = data[data[gene_column] == gene_in_original_order]
                    
                    # Contar muestras únicas afectadas (igual que en summary.py)
                    affected_samples_set = set()
                    
                    for sample_col_name in sample_cols_for_percentage:
                        for _, row_series in gene_data_filtered.iterrows():
                            sample_genotype_value = str(row_series[sample_col_name]).strip().upper()
                            
                            is_mutation_present = False
                            if '|' in sample_genotype_value:
                                alleles = sample_genotype_value.split('|')
                                if len(alleles) >= 2 and alleles[0] != alleles[1]: 
                                    is_mutation_present = True
                            elif '/' in sample_genotype_value: 
                                alleles = sample_genotype_value.split('/')
                                if len(alleles) >= 2 and alleles[0] != alleles[1]:
                                    is_mutation_present = True
                            elif sample_genotype_value not in ["", ".", "0", "0/0", "0|0"] and not pd.isna(row_series[sample_col_name]):
                                is_mutation_present = True
                            
                            if is_mutation_present:
                                affected_samples_set.add(sample_col_name)
                    
                    # Calcular porcentaje exactamente igual que summary.py
                    num_unique_samples_affected = len(affected_samples_set)
                    real_percentage = (num_unique_samples_affected / total_samples_in_dataset) * 100 if total_samples_in_dataset > 0 else 0
                    
                    if real_percentage > 0:
                        # Calcular offset para el texto
                        bar_length = df_percentage.loc[gene].sum()
                        offset = max_percentage * 0.02
                        ax_genes.text(bar_length + offset, i, f'{real_percentage:.1f}%', 
                                    va='center', fontsize=9)
                
                # Remover leyenda del panel lateral (se muestra abajo)
                legend = ax_genes.get_legend()
                if legend is not None:
                    legend.remove()
                
                # Configurar estilo del panel
                ax_genes.spines['top'].set_visible(False)
                ax_genes.spines['right'].set_visible(False)
                ax_genes.spines['left'].set_visible(False)
                ax_genes.tick_params(axis='y', which='both', left=False, right=False)
                ax_genes.tick_params(axis='x', labelsize=9)
                ax_genes.grid(axis='x', alpha=0.3)
                ax_genes.margins(x=0.05, y=0.01)
                
            except Exception as e:
                print(f"Advertencia: Error al crear panel de genes personalizado: {e}")
                import traceback
                traceback.print_exc()
                # Si falla, crear panel vacío
                ax_genes.text(0.5, 0.5, 'Panel de genes\nno disponible', 
                             ha='center', va='center', fontsize=10)
                ax_genes.set_xlim(0, 1)
                ax_genes.set_ylim(0, 1)
                ax_genes.axis('off')
            
            # === PANEL LEYENDA (INFERIOR) ===
            # Crear elementos de leyenda para tipos de variantes
            legend_elements = []
            for variant in sorted(unique_variants):
                if variant != 'None':  # No mostrar 'None' en la leyenda
                    color = color_mapping[variant]
                    label = variant.replace('_', ' ')
                    legend_elements.append(
                        plt.Rectangle((0, 0), 1, 1, facecolor=color, label=label)
                    )
            
            if legend_elements:
                ax_legend.legend(
                    handles=legend_elements,
                    title='Clasificación de Variantes',
                    loc='center',
                    ncol=min(len(legend_elements), 6),
                    fontsize=10,
                    title_fontsize=11,
                    frameon=False
                )
            
            ax_legend.axis('off')
            
            # Título principal de la figura
            fig.suptitle(title, fontsize=16, fontweight='bold', y=0.95)
            
        else:
            # Modo simple: solo heatmap sin paneles adicionales
            fig = ax.get_figure()
            ax_main = ax
        
        # === HEATMAP PRINCIPAL ===
        # Crear colormap personalizado
        colors = [color_mapping[value] for value in unique_values]
        custom_cmap = mcolors.ListedColormap(colors)
        
        # Crear heatmap principal
        im = ax_main.imshow(
            numeric_matrix.values,
            cmap=custom_cmap,
            aspect='auto',
            interpolation='nearest'
        )
        
        # Configurar ejes del heatmap principal
        if ax is not None:  # Solo poner título si no hay panel TMB
            ax_main.set_title(title, fontsize=16, fontweight='bold', pad=20)
        
        # ELIMINAR la etiqueta del eje X para evitar solapamiento con el TMB
        # ax_main.set_xlabel(f'Muestras (n={len(sorted_samples)})', fontsize=12)  # Comentado
        ax_main.set_ylabel('Genes', fontsize=12)
        
        # Configurar ticks y labels
        ax_main.set_xticks(range(len(sorted_samples)))
        ax_main.set_yticks(range(len(top_genes)))
        ax_main.set_yticklabels(top_genes, fontsize=10)
        
        # ELIMINAR completamente las marcas y etiquetas del eje X del waterfall plot
        ax_main.set_xticklabels([''] * len(sorted_samples))  # Etiquetas vacías
        ax_main.tick_params(axis='x', which='both', bottom=False, top=False)  # Sin marcas en X
        
        # Añadir grid sutil
        ax_main.set_xticks(np.arange(-0.5, len(sorted_samples), 1), minor=True)
        ax_main.set_yticks(np.arange(-0.5, len(top_genes), 1), minor=True)
        ax_main.grid(which='minor', color='white', linestyle='-', linewidth=0.5)
        
        # Si es modo simple, añadir leyenda al heatmap
        if ax is not None:
            legend_elements = []
            for variant in sorted(unique_variants):
                if variant != 'None':
                    color = color_mapping[variant]
                    label = variant.replace('_', ' ')
                    legend_elements.append(
                        plt.Rectangle((0, 0), 1, 1, facecolor=color, label=label)
                    )
            
            if legend_elements:
                ax_main.legend(
                    handles=legend_elements,
                    title='Clasificación de Variantes',
                    bbox_to_anchor=(1.05, 1),
                    loc='upper left',
                    fontsize=9,
                    title_fontsize=10
                )
        
        print(f"Oncoplot creado exitosamente:")
        print(f"  - {len(top_genes)} genes")
        print(f"  - {len(sorted_samples)} muestras")
        print(f"  - {len(unique_variants)} tipos de variantes")
        
        return fig
        
    except Exception as e:
        print(f"Error al crear oncoplot: {e}")
        raise 