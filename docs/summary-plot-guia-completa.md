# Summary Plot - Guía Completa y Técnica

Esta guía cubre en detalle todos los aspectos del Summary Plot de pyMut, desde la implementación técnica hasta casos de uso avanzados.

## Arquitectura del Summary Plot

### Estructura Interna

El Summary Plot de pyMut está construido como una composición de 6 visualizaciones independientes organizadas en una grilla 2x3:

```python
# Estructura de la grilla (2 filas x 3 columnas)
[Variant Classification] [Variant Type]       [SNV Class]
[Variants per Sample]    [Classification Box] [Top Genes]
```

### Implementación Técnica

#### Clase Principal: PyMutation.summary_plot()

```python
def summary_plot(self, 
               figsize: Tuple[int, int] = (16, 12),
               title: str = "Mutation Summary",
               max_samples: Optional[int] = 200,
               top_genes_count: int = 10,
               show_interactive: bool = False) -> plt.Figure
```

**Parámetros detallados:**

- `figsize`: Control preciso del tamaño de figura
  - Default: (16, 12) - Balanceado para pantalla y papel
  - Recomendado para publicaciones: (20, 16) o (18, 14)
  - Mínimo recomendado: (12, 8)

- `title`: Título principal de la figura completa
  - Se aplica como `fig.suptitle()` con fontsize=16
  - Posicionamiento automático en la parte superior

- `max_samples`: Control de rendimiento para TMB plot
  - None: Muestra todas las muestras (puede ser lento)
  - Valor recomendado: 100-200 para datasets grandes
  - Toma las primeras N muestras según orden original

- `top_genes_count`: Número de genes en el plot inferior derecho
  - Default: 10 genes
  - Rango recomendado: 5-25 genes
  - Si hay menos genes que el especificado, muestra todos

- `show_interactive`: Control de visualización
  - False: Solo genera la figura (recomendado para scripts)
  - True: Abre ventana interactiva usando manejo event-driven

## Subvisualizaciones Detalladas

### 1. Variant Classification Plot (Superior Izquierda)

**Función:** `create_variant_classification_plot()`

**Propósito:** Muestra la distribución de tipos de mutaciones en el dataset completo.

**Implementación:**
- Gráfico de barras horizontales
- Ordenamiento por frecuencia (más frecuente arriba)
- Colores consistentes con el resto del summary plot

**Interpretación biológica:**
- `Missense_Mutation`: Mutaciones que cambian aminoácidos (más comunes)
- `Silent`: Mutaciones sinónimas (no cambian proteína)
- `Nonsense_Mutation`: Mutaciones que crean codones stop prematuros
- `Frame_Shift_Del/Ins`: Mutaciones que alteran el marco de lectura

### 2. Variant Type Plot (Superior Centro)

**Función:** `create_variant_type_plot()`

**Propósito:** Clasifica las mutaciones por mecanismo molecular.

**Categorías típicas:**
- `SNP`: Single Nucleotide Polymorphisms (más comunes)
- `DEL`: Deleciones (pérdida de nucleótidos)
- `INS`: Inserciones (ganancia de nucleótidos)
- `DNP`: Dinucleotide Polymorphisms
- `TNP`: Trinucleotide Polymorphisms

### 3. SNV Class Plot (Superior Derecha)

**Función:** `create_snv_class_plot()`

**Propósito:** Analiza los patrones de cambios nucleotídicos específicos.

**Parámetros técnicos:**
```python
ref_column="REF"    # Columna con alelo de referencia
alt_column="ALT"    # Columna con alelo alternativo
```

**Categorías generadas:**
- Transiciones: A>G, G>A, C>T, T>C (más comunes, ~66%)
- Transversiones: A>C, A>T, G>C, G>T, C>A, C>G, T>A, T>G

**Análisis de signatures mutacionales:** Patrones pueden indicar:
- C>T: Deaminación espontánea de citosina
- G>A: Replicación errónea
- A>T: Daño oxidativo

### 4. Variants per Sample Plot / TMB (Inferior Izquierda)

**Función:** `create_variants_per_sample_plot()`

**Propósito:** Tumor Mutation Burden - distribución de carga mutacional por muestra.

**Características técnicas:**
- Barras apiladas mostrando composición de mutaciones
- Línea roja horizontal: mediana del TMB
- Ordenamiento por TMB total (mayor a menor)
- Limitación configurable con `max_samples`

**Interpretación clínica:**
- TMB alto: Posible beneficio de inmunoterapia
- TMB bajo: Tumor con pocas mutaciones
- Patrones de TMB: Identificación de subgrupos moleculares

### 5. Variant Classification Summary / Boxplot (Inferior Centro)

**Función:** `create_variant_classification_summary_plot()`

**Propósito:** Distribución estadística de cada tipo de mutación entre muestras.

**Elementos del boxplot:**
- Caja: Rango intercuartílico (IQR)
- Línea central: Mediana
- Whiskers: 1.5 × IQR desde Q1 y Q3
- Puntos: Outliers (muestras atípicas)

**Aplicaciones:**
- Identificar tipos de mutación con alta variabilidad
- Detectar muestras con patrones mutacionales atípicos
- Comparar consistencia entre tipos de mutación

### 6. Top Mutated Genes Plot (Inferior Derecha)

**Función:** `create_top_mutated_genes_plot()`

**Propósito:** Identificar genes más frecuentemente alterados.

**Modos de conteo:**
```python
mode="variants"  # Cuenta total de mutaciones por gen
mode="samples"   # Cuenta muestras afectadas por gen (prevalencia)
```

**En Summary Plot:** Siempre usa `mode="variants"`

**Interpretación:**
- Genes con muchas mutaciones: Posibles hotspots
- Genes driver: Frecuentemente mutados en cáncer
- Genes passenger: Mutaciones no funcionales

## Sistema de Colores y Leyenda

### Mapeo de Colores Consistente

El Summary Plot genera un mapa de colores único que se aplica a todas las subvisualizaciones:

```python
# Generación automática de colores
unique_variants = data['Variant_Classification'].unique()
cmap = plt.colormaps['tab20']
variant_color_map = {variant: cmap(i % 20) for i, variant in enumerate(unique_variants)}
```

### Leyenda Unificada

- Posición: Parte inferior de la figura
- Ordenamiento: Por frecuencia de mutaciones (descendente)
- Elementos: Solo variantes presentes en los datos
- Formato: Máximo 5 columnas para legibilidad

## Manejo de Datos y Preprocesamiento

### Detección Automática de Formato

pyMut detecta automáticamente el formato de los datos:

**Formato Largo (Long Format):**
```
Hugo_Symbol | Variant_Classification | Tumor_Sample_Barcode | REF | ALT
GENE1      | Missense_Mutation      | SAMPLE_001          | A   | G
```

**Formato Ancho (Wide Format - TCGA):**
```
Hugo_Symbol | Variant_Classification | TCGA-AB-2001 | TCGA-AB-2002
GENE1      | Missense_Mutation      | A|G          | C|T
```

### Procesamiento de Datos

#### 1. Extracción de Variant Classification

```python
# Si no existe la columna, extrae de FUNCOTATION
processed_data = extract_variant_classifications(
    self.data, 
    variant_column="Variant_Classification",
    funcotation_column="FUNCOTATION"
)
```

#### 2. Extracción de Variant Types

```python
processed_data = extract_variant_types(
    processed_data,
    variant_column="Variant_Type",
    funcotation_column="FUNCOTATION"
)
```

#### 3. Normalización de Columnas

- Búsqueda case-insensitive de nombres de columna
- Mapeo automático de variantes de nombres comunes
- Validación de existencia de datos requeridos

## Optimización de Rendimiento

### Estrategias para Datasets Grandes

#### 1. Limitación de Muestras
```python
# Limita el TMB plot para mejorar rendimiento
fig = py_mut.summary_plot(max_samples=100)
```

#### 2. Filtrado Previo de Datos
```python
# Filtrar antes de crear el objeto PyMutation
relevant_variants = ['Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del']
filtered_data = data[data['Variant_Classification'].isin(relevant_variants)]
py_mut = PyMutation(filtered_data)
```

#### 3. Reducción de Genes Mostrados
```python
# Menos genes = procesamiento más rápido
fig = py_mut.summary_plot(top_genes_count=10)
```

### Benchmarks de Rendimiento

| Tamaño Dataset | Muestras | Tiempo (seg) | Memoria (MB) |
|---------------|----------|--------------|--------------|
| Pequeño       | < 100    | < 5         | < 100        |
| Mediano       | 100-500  | 5-15        | 100-300      |
| Grande        | 500+     | 15-60       | 300-800      |

## Gestión de Calidad de Imagen

### Configuración Automática de Alta Calidad

```python
# Configurar una vez al inicio del script
PyMutation.configure_high_quality_plots()

# Configura automáticamente:
# - DPI = 300
# - bbox_inches = 'tight'
# - facecolor = 'white'
# - edgecolor = 'none'
```

### Métodos de Guardado

#### 1. Automático (Recomendado)
```python
PyMutation.configure_high_quality_plots()
fig = py_mut.summary_plot()
fig.savefig("summary.png")  # Automáticamente alta calidad
```

#### 2. Manual
```python
fig.savefig("summary.png", dpi=300, bbox_inches='tight')
```

#### 3. Método Centralizado de pyMut
```python
py_mut.save_figure(fig, "summary.png", dpi=600)  # Ultra alta calidad
```

### Formatos Recomendados por Uso

| Formato | Uso Principal | DPI | Tamaño |
|---------|---------------|-----|--------|
| PNG     | Presentaciones, web | 300 | Medio |
| PDF     | Publicaciones científicas | 300 | Grande |
| SVG     | Edición, vectorial | N/A | Pequeño |
| TIFF    | Impresión profesional | 600 | Muy grande |

## Manejo Interactivo

### Event-Driven Approach

pyMut implementa un sistema de ventanas interactivas que evita los problemas comunes de matplotlib:

```python
def _show_figure_interactive(self, figure):
    """Muestra figura usando enfoque event-driven"""
    # Configura callback para cierre de ventana
    # Maneja eventos sin bloquear la aplicación
    # Permite cerrar ventanas sin cuelgues
```

### Ventajas del Sistema

1. **No Global State**: No afecta otras figuras
2. **Responsive Windows**: Ventanas que responden correctamente
3. **Proper Cleanup**: Limpieza automática de recursos
4. **No Hanging**: Sin bloqueos de aplicación

## Casos de Uso Avanzados

### 1. Análisis Comparativo de Cohorts

```python
# Generar summary plots para múltiples cohorts
cohorts = ['TCGA-LAML', 'TCGA-GBM', 'TCGA-BRCA']

for cohort in cohorts:
    data = load_cohort_data(cohort)
    py_mut = PyMutation(data)
    
    fig = py_mut.summary_plot(
        title=f"{cohort} Mutation Landscape",
        figsize=(18, 14),
        max_samples=200,
        top_genes_count=20
    )
    
    fig.savefig(f"{cohort}_summary.png")
```

### 2. Análisis de Subgrupos Moleculares

```python
# Analizar subgrupos específicos
subgroups = data.groupby('Molecular_Subtype')

for subtype, subtype_data in subgroups:
    py_mut = PyMutation(subtype_data)
    
    fig = py_mut.summary_plot(
        title=f"Mutations in {subtype} Subtype",
        top_genes_count=15
    )
    
    fig.savefig(f"subtype_{subtype}_summary.png")
```

### 3. Análisis Temporal o Longitudinal

```python
# Comparar timepoints o tratamientos
timepoints = ['Baseline', 'Post_Treatment']

for tp in timepoints:
    tp_data = data[data['Timepoint'] == tp]
    py_mut = PyMutation(tp_data)
    
    fig = py_mut.summary_plot(
        title=f"Mutations at {tp}",
        figsize=(16, 12)
    )
    
    fig.savefig(f"mutations_{tp.lower()}.png")
```

## Solución de Problemas Técnicos

### Problemas de Memoria

**Síntoma:** Error de memoria con datasets grandes
**Solución:**
```python
# Reducir tamaño de datos procesados
fig = py_mut.summary_plot(
    max_samples=50,      # Limitar muestras
    top_genes_count=10   # Limitar genes
)

# O filtrar datos previamente
filtered_data = data.sample(n=1000)  # Muestra aleatoria
```

### Problemas de Visualización

**Síntoma:** Plots vacíos o errores de columnas
**Diagnóstico:**
```python
print(f"Columnas disponibles: {data.columns.tolist()}")
print(f"Clasificaciones: {data['Variant_Classification'].value_counts()}")
print(f"Tipos de variante: {data.get('Variant_Type', 'No disponible')}")
```

### Problemas de Rendimiento

**Síntoma:** Generación muy lenta
**Optimizaciones:**
```python
# 1. Reducir complejidad visual
fig = py_mut.summary_plot(figsize=(12, 8))  # Figura más pequeña

# 2. Procesar menos datos
data_subset = data.head(5000)  # Primeras 5000 mutaciones

# 3. Usar formato más eficiente
data = data.astype({'Hugo_Symbol': 'category'})  # Categorías para strings
```

## Integración con Pipelines de Análisis

### Script de Análisis Completo

```python
#!/usr/bin/env python3
"""
Pipeline completo de análisis mutacional con pyMut
"""
import pandas as pd
from pyMut import PyMutation
import os

def analyze_mutations(input_file, output_dir):
    """Análisis completo de mutaciones"""
    
    # Configurar calidad
    PyMutation.configure_high_quality_plots()
    
    # Cargar y validar datos
    print(f"Cargando datos de {input_file}...")
    data = pd.read_csv(input_file, sep='\t')
    
    print(f"Dataset: {data.shape[0]} mutaciones, {data.shape[1]} columnas")
    print(f"Muestras únicas: {data['Tumor_Sample_Barcode'].nunique()}")
    print(f"Genes únicos: {data['Hugo_Symbol'].nunique()}")
    
    # Crear análisis
    py_mut = PyMutation(data)
    
    # Generar summary plot principal
    summary_fig = py_mut.summary_plot(
        title="Comprehensive Mutation Analysis",
        figsize=(20, 16),
        max_samples=200,
        top_genes_count=25
    )
    
    # Guardar en múltiples formatos
    os.makedirs(output_dir, exist_ok=True)
    
    summary_fig.savefig(f"{output_dir}/mutation_summary.png")
    summary_fig.savefig(f"{output_dir}/mutation_summary.pdf")
    summary_fig.savefig(f"{output_dir}/mutation_summary.svg")
    
    # Generar plots individuales para análisis detallado
    individual_plots = [
        ('variant_classification', py_mut.variant_classification_plot),
        ('tmb', py_mut.variants_per_sample_plot),
        ('top_genes_variants', lambda: py_mut.top_mutated_genes_plot(mode='variants')),
        ('top_genes_samples', lambda: py_mut.top_mutated_genes_plot(mode='samples'))
    ]
    
    for name, plot_func in individual_plots:
        fig = plot_func()
        fig.savefig(f"{output_dir}/{name}.png")
    
    print(f"✅ Análisis completado. Resultados en {output_dir}/")

if __name__ == "__main__":
    analyze_mutations("mutations.tsv", "analysis_results")
```

### Integración con Jupyter Notebooks

```python
# Configuración para notebooks
%matplotlib inline
from pyMut import PyMutation

# Configurar una vez por sesión
PyMutation.configure_high_quality_plots()

# Análisis interactivo
py_mut = PyMutation(data)

# Ver interactivamente
fig = py_mut.summary_plot(show_interactive=True)

# Guardar para reporte
fig.savefig("notebook_summary.png")
```

## Extensiones y Personalización

### Modificación de Colores

```python
# Personalizar mapa de colores
custom_colors = {
    'Missense_Mutation': '#FF6B6B',
    'Silent': '#4ECDC4', 
    'Nonsense_Mutation': '#45B7D1',
    'Frame_Shift_Del': '#96CEB4'
}

# Aplicar en plots individuales
fig, ax = plt.subplots()
create_variant_classification_plot(data, ax=ax, color_map=custom_colors)
```

### Métricas Personalizadas

```python
# Calcular métricas adicionales
def calculate_mutation_metrics(data):
    """Calcula métricas personalizadas"""
    metrics = {
        'total_mutations': len(data),
        'unique_genes': data['Hugo_Symbol'].nunique(),
        'unique_samples': data['Tumor_Sample_Barcode'].nunique(),
        'avg_mutations_per_sample': len(data) / data['Tumor_Sample_Barcode'].nunique(),
        'most_mutated_gene': data['Hugo_Symbol'].value_counts().index[0]
    }
    return metrics

# Usar con summary plot
metrics = calculate_mutation_metrics(data)
fig = py_mut.summary_plot(
    title=f"Mutations Analysis - {metrics['total_mutations']} total mutations"
)
```

Esta guía completa cubre todos los aspectos técnicos y prácticos del Summary Plot, proporcionando una referencia exhaustiva para usuarios desde básicos hasta avanzados. 