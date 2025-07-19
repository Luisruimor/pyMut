# Summary Plot - Complete Mutation Analysis

El **Summary Plot** es la visualización principal de pyMut que combina múltiples análisis de mutaciones en una sola figura comprehensiva.

## ¿Qué es el Summary Plot?

Es un panel que incluye 6 visualizaciones diferentes en una sola figura:

1. **Variant Classification** - Distribución de tipos de mutación (Missense, Nonsense, etc.)
2. **Variant Type** - Distribución de tipos de variante (SNP, INS, DEL, etc.)  
3. **SNV Class** - Distribución de cambios nucleotídicos (A>G, C>T, etc.)
4. **Variants per Sample (TMB)** - Número de mutaciones por muestra con línea de mediana
5. **Variant Classification Summary** - Boxplot de clasificaciones por muestra
6. **Top Mutated Genes** - Genes más frecuentemente mutados

## Uso Básico

```python
from pyMut import PyMutation
import pandas as pd

# Paso 1: Configurar alta calidad (recomendado)
PyMutation.configure_high_quality_plots()

# Paso 2: Cargar datos y crear objeto
data = pd.read_csv("mutations.tsv", sep='\t')
py_mut = PyMutation(data)

# Paso 3: Generar summary plot
fig = py_mut.summary_plot()

# Paso 4: Guardar figura (automáticamente alta calidad)
fig.savefig("summary_plot.png")
```

## Personalización Avanzada

### Parámetros principales

```python
fig = py_mut.summary_plot(
    title="Mi Análisis de Mutaciones",    # Título personalizado
    figsize=(18, 14),                     # Tamaño de figura (más grande)
    max_samples=100,                      # Máximo muestras en TMB plot
    top_genes_count=10                    # Número de genes top a mostrar
)
```

### Control de rendimiento para datasets grandes

```python
# Para datasets con muchas muestras
fig = py_mut.summary_plot(
    max_samples=50,        # Limita el TMB plot a 50 muestras
    top_genes_count=10     # Muestra solo los 10 genes más mutados
)

# Para figuras más grandes y detalladas
fig = py_mut.summary_plot(
    figsize=(20, 16),      # Figura muy grande para publicaciones
    top_genes_count=20     # Más genes para análisis detallado
)
```

### Guardar en diferentes formatos

```python
# PNG para presentaciones (por defecto)
fig.savefig("summary.png")

# PDF para publicaciones científicas
fig.savefig("summary.pdf")

# SVG para edición vectorial
fig.savefig("summary.svg")

# Ultra alta calidad para impresión
py_mut.save_figure(fig, "summary_ultra.png", dpi=600)
```

## Formato de Datos

pyMut funciona con datos en formato largo (recomendado):

```
Hugo_Symbol | Variant_Classification | Tumor_Sample_Barcode | REF | ALT
GENE1      | Missense_Mutation      | SAMPLE_001          | A   | G
GENE2      | Nonsense_Mutation      | SAMPLE_001          | C   | T
```

**Columnas requeridas:**
- `Hugo_Symbol` - Símbolo del gen
- `Variant_Classification` - Tipo de mutación

**Columnas opcionales pero recomendadas:**
- `Tumor_Sample_Barcode` - Identificador de muestra
- `REF` / `ALT` - Alelos de referencia y alternativo

### Soporte para formato TCGA (ancho)
```python
# pyMut detecta automáticamente columnas de muestras TCGA
# Ejemplo: TCGA-AB-2001, TCGA-AB-2002, etc.
# Se procesan automáticamente sin configuración adicional
```

## Interpretación de Resultados

### Variant Classification
- **Missense_Mutation**: Cambia aminoácido (potencialmente funcional)
- **Silent**: No cambia aminoácido (generalmente neutral)
- **Nonsense_Mutation**: Crea codón de parada (muy impactante)
- **Frame_Shift_Del/Ins**: Cambia marco de lectura (muy impactante)
- **Splice_Site**: Afecta empalme de RNA (potencialmente impactante)

### TMB Plot (Tumor Mutation Burden)
- **Altura de barras**: Número total de mutaciones por muestra
- **Colores en barras**: Composición por tipos de mutación
- **Línea roja**: Mediana de todas las muestras
- **Utilidad**: Identificar muestras con alta carga mutacional

### Variant Classification Summary (Boxplot)
- **Cajas**: Distribución de cada tipo de mutación entre muestras
- **Línea central**: Mediana del tipo de mutación
- **Whiskers**: Rango de valores típicos
- **Puntos**: Muestras atípicas (outliers)

### Top Genes
- **Modo "variants"**: Genes con más mutaciones totales
- **Modo "samples"**: Genes mutados en más muestras (prevalencia)
- **Utilidad**: Identificar posibles genes driver o hotspots

## Características Técnicas

### Consistencia de colores
- Mapeo de colores coherente entre todas las subgráficas
- Leyenda unificada en la parte inferior
- Ordenamiento por frecuencia de mutaciones

### Gestión de memoria
- Procesamiento optimizado para datasets grandes
- Parámetro `max_samples` para limitar uso de memoria
- Detección automática de formato de datos

### Calidad de imagen
- DPI 300 por defecto con `configure_high_quality_plots()`
- Márgenes optimizados automáticamente
- Soporte para múltiples formatos vectoriales

## Ejemplo Completo con Datos TCGA

```python
import pandas as pd
from pyMut import PyMutation

# Cargar datos de ejemplo TCGA-LAML
data = pd.read_csv("src/pyMut/data/examples/tcga_laml_converted.tsv", sep='\t')

# Crear análisis
py_mut = PyMutation(data)

# Configurar alta calidad
PyMutation.configure_high_quality_plots()

# Generar summary plot personalizado para publicación
summary_fig = py_mut.summary_plot(
    title="TCGA-LAML Mutation Landscape Analysis",
    figsize=(20, 16),           # Figura grande para publicación
    max_samples=150,            # Mostrar hasta 150 muestras
    top_genes_count=20,         # Top 20 genes más mutados
    show_interactive=False      # No mostrar ventana interactiva
)

# Guardar en múltiples formatos
summary_fig.savefig("tcga_laml_summary.png")      # PNG alta calidad
summary_fig.savefig("tcga_laml_summary.pdf")      # PDF para publicación
summary_fig.savefig("tcga_laml_summary.svg")      # SVG para edición

print("✅ Summary plot generado exitosamente!")
```

## Solución de Problemas

### Datos vacíos o gráficos inesperados
```python
# Verificar estructura de datos
print(f"Forma de datos: {data.shape}")
print(f"Columnas: {data.columns.tolist()}")
print(f"Clasificaciones únicas: {data['Variant_Classification'].value_counts()}")
print(f"Genes únicos: {data['Hugo_Symbol'].nunique()}")
```

### Demasiadas muestras en TMB plot
```python
# Limitar número de muestras mostradas
fig = py_mut.summary_plot(max_samples=50)

# O verificar número de muestras primero
n_samples = data['Tumor_Sample_Barcode'].nunique()
print(f"Total de muestras: {n_samples}")
```

### Imágenes de baja calidad
```python
# Siempre configurar alta calidad al inicio
PyMutation.configure_high_quality_plots()

# O manualmente para un caso específico
fig.savefig("plot.png", dpi=300, bbox_inches='tight')
```

### Rendimiento con datasets grandes
```python
# Filtrar datos antes de crear el plot
filtered_data = data[data['Variant_Classification'].isin([
    'Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del'
])]
py_mut = PyMutation(filtered_data)

# O usar parámetros de rendimiento
fig = py_mut.summary_plot(
    max_samples=100,       # Limitar muestras
    top_genes_count=10     # Limitar genes
)
``` 