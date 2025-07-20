# Guía Completa del Oncoplot

El **oncoplot** (también conocido como waterfall plot en algunos contextos) es una visualización fundamental en genómica del cáncer que muestra los patrones de mutación a través de muestras y genes en formato heatmap. Esta guía cubre todo lo que necesitas saber para crear y personalizar oncoplots usando pyMut.

## Tabla de Contenidos

1. [Introducción](#introducción)
2. [Instalación y Configuración](#instalación-y-configuración)
3. [Uso Básico](#uso-básico)
4. [Parámetros Detallados](#parámetros-detallados)
5. [Formatos de Datos Soportados](#formatos-de-datos-soportados)
6. [Ejemplos Avanzados](#ejemplos-avanzados)
7. [Personalización Visual](#personalización-visual)
8. [Interpretación del Oncoplot](#interpretación-del-oncoplot)
9. [Solución de Problemas](#solución-de-problemas)
10. [Optimización de Rendimiento](#optimización-de-rendimiento)

## Introducción

El oncoplot es una herramienta esencial para visualizar paisajes mutacionales en datos genómicos. Proporciona una vista panorámica de:

- **Genes más frecuentemente mutados**
- **Patrones de co-ocurrencia de mutaciones**
- **Carga mutacional por muestra**
- **Tipos de mutaciones por color**
- **Distribución de mutaciones across samples**

### Características Principales

- **Detección automática** de columnas de muestra (TCGA y formato .GT)
- **Soporte multiformato** para genotipos (A|G, A/G, etc.)
- **Detección de Multi_Hit** para muestras con múltiples mutaciones
- **Esquemas de colores estándar** para tipos de mutación
- **Ordenamiento inteligente** por frecuencia de mutación
- **Alta calidad de exportación** (PNG, PDF, SVG)

## Instalación y Configuración

### Requisitos Previos

```bash
pip install pandas matplotlib seaborn numpy
```

### Configuración de Alta Calidad

```python
from pyMut import PyMutation

# Configurar alta calidad automáticamente (recomendado)
PyMutation.configure_high_quality_plots()
```

Esta configuración aplica automáticamente:
- DPI = 300 para imágenes nítidas
- bbox_inches = 'tight' para layouts optimizados
- Configuraciones de fuente mejoradas

## Uso Básico

### Ejemplo Mínimo

```python
import pandas as pd
from pyMut import PyMutation

# Cargar datos
data = pd.read_csv('mutations.tsv', sep='\t')

# Crear objeto PyMutation
py_mut = PyMutation(data)

# Generar oncoplot básico
fig = py_mut.oncoplot()

# Guardar
fig.savefig('oncoplot.png')
```

### Ejemplo con Personalización

```python
# Oncoplot personalizado
fig = py_mut.oncoplot(
    title="Paisaje Mutacional - TCGA LAML",
    top_genes_count=20,
    max_samples=100,
    figsize=(16, 10)
)
```

## Parámetros Detallados

### Parámetros del Método `oncoplot()`

| Parámetro | Tipo | Valor por Defecto | Descripción |
|-----------|------|-------------------|-------------|
| `figsize` | `Tuple[int, int]` | `(16, 10)` | Tamaño de figura (ancho, alto) en pulgadas |
| `title` | `str` | `"Oncoplot"` | Título del gráfico |
| `gene_column` | `str` | `"Hugo_Symbol"` | Columna con símbolos de genes |
| `variant_column` | `str` | `"Variant_Classification"` | Columna con tipos de mutación |
| `ref_column` | `str` | `"REF"` | Columna con alelos de referencia |
| `alt_column` | `str` | `"ALT"` | Columna con alelos alternativos |
| `top_genes_count` | `int` | `10` | Número de genes más mutados a mostrar |
| `max_samples` | `int` | `180` | Número máximo de muestras a mostrar |

### Validación de Parámetros

El método incluye validación robusta:

```python
# ✓ Válido
py_mut.oncoplot(top_genes_count=15, max_samples=50)

# ✗ Inválido - lanzará ValueError
py_mut.oncoplot(top_genes_count=0)         # Debe ser > 0
py_mut.oncoplot(max_samples=-1)            # Debe ser > 0
py_mut.oncoplot(figsize=(0, 10))           # Dimensiones deben ser > 0
```

## Formatos de Datos Soportados

### Estructura de Datos Requerida

El DataFrame debe contener al menos estas columnas:

#### Columnas Obligatorias
- `Hugo_Symbol`: Símbolos de genes (ej: "TP53", "KRAS")
- `Variant_Classification`: Tipos de mutación (ej: "Missense_Mutation")
- `REF`: Alelos de referencia (ej: "A", "T")
- `ALT`: Alelos alternativos (ej: "G", "C")

#### Columnas de Muestra
Las columnas de muestra se detectan automáticamente usando patrones:
- **Formato TCGA**: `TCGA-AB-1234`, `TCGA-CD-5678`
- **Formato GT**: `sample1.GT`, `sample2.GT`

### Formatos de Genotipo Soportados

El oncoplot soporta múltiples formatos de genotipo:

```python
# Formatos válidos:
"A|G"    # Formato estándar pipe-separated
"A/G"    # Formato slash-separated
"A:G"    # Colon-separated
"A;G"    # Semicolon-separated
"A,G"    # Comma-separated

# Casos especiales:
"A|A"    # Homocigoto referencia (no mutación)
"G|G"    # Homocigoto alternativo (mutación)
"."      # Datos faltantes
""       # Datos vacíos
```

### Estructura de Ejemplo

```python
import pandas as pd

data = pd.DataFrame({
    'Hugo_Symbol': ['TP53', 'KRAS', 'PIK3CA'],
    'Variant_Classification': ['Missense_Mutation', 'Nonsense_Mutation', 'In_Frame_Del'],
    'REF': ['A', 'C', 'G'],
    'ALT': ['G', 'T', 'A'],
    'TCGA-AB-1234': ['A|G', 'C|C', 'G|A'],    # Mutado en TP53 y PIK3CA
    'TCGA-CD-5678': ['A|A', 'C|T', 'G|G'],    # Mutado en KRAS y PIK3CA
    'TCGA-EF-9012': ['A|G', 'C|T', 'G|G']     # Mutado en todos los genes
})
```

## Ejemplos Avanzados

### Enfoque en Genes Específicos

```python
# Mostrar solo genes altamente mutados
fig = py_mut.oncoplot(
    title="Genes Driver Principales",
    top_genes_count=5,
    max_samples=200,
    figsize=(20, 6)
)
```

### Vista Panorámica

```python
# Vista amplia del paisaje mutacional
fig = py_mut.oncoplot(
    title="Paisaje Mutacional Completo",
    top_genes_count=30,
    max_samples=50,
    figsize=(18, 12)
)
```

### Exportación Multi-formato

```python
fig = py_mut.oncoplot(title="Análisis Principal")

# Exportar en múltiples formatos
for fmt in ['png', 'pdf', 'svg']:
    output_file = f'oncoplot.{fmt}'
    if fmt == 'svg':
        fig.savefig(output_file, format=fmt)
    else:
        fig.savefig(output_file, format=fmt, dpi=300)
```

### Análisis Interactivo

```python
# Para exploración interactiva
fig = py_mut.oncoplot(
    title="Exploración Interactiva",
    top_genes_count=25
)
```

## Personalización Visual

### Esquema de Colores

El oncoplot usa colores estándar para tipos de mutación comunes:

| Tipo de Mutación | Color | Descripción |
|------------------|-------|-------------|
| Missense_Mutation | Azul | Mutaciones que cambian aminoácido |
| Nonsense_Mutation | Rojo | Mutaciones que introducen stop codon |
| Frame_Shift_Del | Verde | Deleciones que alteran marco de lectura |
| Frame_Shift_Ins | Púrpura | Inserciones que alteran marco de lectura |
| In_Frame_Del | Cian claro | Deleciones que mantienen marco |
| In_Frame_Ins | Naranja | Inserciones que mantienen marco |
| Splice_Site | Rosa | Mutaciones en sitios de splicing |
| Multi_Hit | Negro | Múltiples tipos en misma muestra/gen |
| None | Gris claro | Sin mutación detectada |

### Tamaños de Figura Recomendados

```python
# Para pocos genes, muchas muestras
figsize=(20, 6)   # Formato panorámico

# Para muchos genes, pocas muestras  
figsize=(12, 16)  # Formato vertical

# Balanceado general
figsize=(16, 10)  # Formato estándar

# Para presentaciones
figsize=(24, 12)  # Extra grande
```

## Interpretación del Oncoplot

### Lectura del Gráfico

1. **Eje Y (Genes)**: Ordenados por frecuencia de mutación (más mutados arriba)
2. **Eje X (Muestras)**: Ordenadas por carga mutacional total
3. **Colores**: Indican tipo de mutación según leyenda
4. **Celdas grises**: No hay mutación detectada
5. **Celdas negras**: Múltiples tipos de mutación (Multi_Hit)

### Patrones Importantes

- **Genes frecuentemente mutados**: Aparecen en la parte superior
- **Muestras con alta carga mutacional**: Hacia la izquierda
- **Co-ocurrencias**: Genes mutados en las mismas muestras
- **Exclusividad mutua**: Genes que raramente se mutan juntos

### Análisis de Multi_Hit

Cuando una muestra tiene múltiples tipos de mutación en el mismo gen:

```python
# Ejemplo: TP53 con Missense y Nonsense en la misma muestra
# Se mostrará como Multi_Hit (negro)
data_multihit = pd.DataFrame({
    'Hugo_Symbol': ['TP53', 'TP53'],
    'Variant_Classification': ['Missense_Mutation', 'Nonsense_Mutation'],
    'REF': ['A', 'C'],
    'ALT': ['G', 'T'],
    'TCGA-AB-1234': ['A|G', 'C|T']  # Ambas mutaciones presentes
})
```

## Solución de Problemas

### Errores Comunes

#### 1. Columnas Faltantes

```python
# Error: ValueError: Faltan las siguientes columnas requeridas: ['REF']
# Solución: Verificar que todas las columnas requeridas estén presentes
required_columns = ['Hugo_Symbol', 'Variant_Classification', 'REF', 'ALT']
missing = [col for col in required_columns if col not in data.columns]
print(f"Columnas faltantes: {missing}")
```

#### 2. No se Detectan Muestras

```python
# Error: No se pudieron detectar columnas de muestra automáticamente
# Solución: Verificar formato de nombres de columnas
sample_cols = [col for col in data.columns if col.startswith('TCGA-')]
print(f"Columnas TCGA encontradas: {sample_cols}")
```

#### 3. No Hay Datos de Mutación

```python
# Error: No hay datos de mutación para visualizar
# Solución: Verificar que hay mutaciones detectables
mutations_detected = 0
for col in sample_cols[:5]:  # Revisar primeras 5 muestras
    mutations = data[col].apply(lambda x: '|' in str(x) and 'A|A' not in str(x)).sum()
    mutations_detected += mutations
print(f"Mutaciones detectadas en primeras 5 muestras: {mutations_detected}")
```

#### 4. Parámetros Inválidos

```python
# Error: ValueError: top_genes_count debe ser un entero positivo
# Solución: Usar valores válidos
assert top_genes_count > 0, "top_genes_count debe ser > 0"
assert max_samples > 0, "max_samples debe ser > 0"
assert len(figsize) == 2 and all(x > 0 for x in figsize), "figsize debe ser (ancho, alto) con valores > 0"
```

### Debugging Paso a Paso

```python
# 1. Verificar estructura de datos
print("Columnas disponibles:", data.columns.tolist())
print("Forma de datos:", data.shape)

# 2. Verificar detección de muestras
from pyMut.visualizations.oncoplot import detect_sample_columns
try:
    samples = detect_sample_columns(data)
    print(f"Muestras detectadas: {len(samples)}")
except ValueError as e:
    print(f"Error en detección: {e}")

# 3. Verificar datos de mutación
from pyMut.visualizations.oncoplot import process_mutation_matrix
try:
    matrix, counts = process_mutation_matrix(data)
    print(f"Matriz generada: {matrix.shape}")
    print(f"Genes con mutaciones: {len([g for g, c in counts.items() if c > 0])}")
except Exception as e:
    print(f"Error en procesamiento: {e}")

# 4. Generar oncoplot
try:
    py_mut = PyMutation(data)
    fig = py_mut.oncoplot(top_genes_count=5)
    print("✓ Oncoplot generado exitosamente")
except Exception as e:
    print(f"Error en oncoplot: {e}")
```

### Problemas de Rendimiento

#### Datasets Grandes

```python
# Para datasets con >1000 genes o >500 muestras
fig = py_mut.oncoplot(
    top_genes_count=20,     # Limitar genes
    max_samples=100,        # Limitar muestras  
    figsize=(16, 8)         # Tamaño moderado
)
```

#### Memoria Insuficiente

```python
# Procesar en chunks para datasets muy grandes
chunk_size = 50
for i in range(0, len(data), chunk_size):
    chunk = data.iloc[i:i+chunk_size]
    py_mut_chunk = PyMutation(chunk)
    fig = py_mut_chunk.oncoplot(
        title=f"Chunk {i//chunk_size + 1}",
        max_samples=chunk_size
    )
    fig.savefig(f'oncoplot_chunk_{i//chunk_size + 1}.png')
```

## Optimización de Rendimiento

### Mejores Prácticas

1. **Filtrar datos antes**: Remover genes/muestras irrelevantes
2. **Usar parámetros apropiados**: Limitar `top_genes_count` y `max_samples`
3. **Configurar calidad una vez**: Llamar `configure_high_quality_plots()` al inicio
4. **Cerrar figuras**: Usar `plt.close(fig)` para liberar memoria

### Parámetros de Rendimiento

```python
# Para análisis exploratorio rápido
fig = py_mut.oncoplot(
    top_genes_count=10,
    max_samples=50,
    figsize=(12, 6)
)

# Para análisis final de alta calidad
fig = py_mut.oncoplot(
    top_genes_count=30,
    max_samples=200,
    figsize=(20, 12)
)
```

### Monitoreo de Recursos

```python
import time
import psutil

# Monitorear tiempo y memoria
start_time = time.time()
memory_before = psutil.Process().memory_info().rss / 1024 / 1024  # MB

fig = py_mut.oncoplot()

end_time = time.time()
memory_after = psutil.Process().memory_info().rss / 1024 / 1024  # MB

print(f"Tiempo: {end_time - start_time:.2f} segundos")
print(f"Memoria usada: {memory_after - memory_before:.1f} MB")
```

---

## Más Recursos

- **Guía de Instalación**: [docs/index.md](index.md)
- **Código Fuente**: [src/pyMut/visualizations/oncoplot.py](../src/pyMut/visualizations/oncoplot.py)
- **Ejemplos Completos**: [examples/example_oncoplot.py](../examples/example_oncoplot.py)
- **Tests**: [tests/test_oncoplot.py](../tests/test_oncoplot.py)

Para preguntas específicas o reportar problemas, consulta la documentación del proyecto o contacta al equipo de desarrollo. 