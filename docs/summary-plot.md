# Summary Plot - Complete Mutation Analysis

El **Summary Plot** es la visualización principal de pyMut que combina múltiples análisis de mutaciones en una sola figura comprehensiva.

## ¿Qué es el Summary Plot?

Es un panel que incluye 6 visualizaciones diferentes en una sola figura:

1. **Variant Classification** - Distribución de tipos de mutación (Missense, Nonsense, etc.)
2. **Variant Type** - Distribución de tipos de variante (SNP, INS, DEL, etc.)  
3. **SNV Class** - Distribución de cambios nucleotídicos (A>G, C>T, etc.)
4. **Variants per Sample (TMB)** - Número de mutaciones por muestra
5. **Variant Classification Summary** - Boxplot de clasificaciones por muestra
6. **Top Mutated Genes** - Genes más frecuentemente mutados

## Uso Básico

```python
from pyMut import PyMutation
import pandas as pd

# Cargar datos
data = pd.read_csv("mutations.tsv", sep='\t')

# Crear objeto PyMutation
py_mut = PyMutation(data)

# Configurar alta calidad (recomendado)
PyMutation.configure_high_quality_plots()

# Generar summary plot
fig = py_mut.summary_plot()

# Guardar figura
fig.savefig("summary_plot.png")
```

## Personalización

### Parámetros principales

```python
fig = py_mut.summary_plot(
    title="Mi Análisis de Mutaciones",    # Título personalizado
    figsize=(16, 12),                     # Tamaño de figura
    max_samples=100,                      # Máximo muestras en TMB plot
    top_genes_count=10,                   # Número de genes top a mostrar
    show_interactive=True                 # Mostrar interactivamente
)
```

### Guardar en diferentes formatos

```python
# PNG para presentaciones
fig.savefig("summary.png")

# PDF para publicaciones
fig.savefig("summary.pdf", dpi=300)

# SVG para edición
fig.savefig("summary.svg")
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

## Interpretación

### Variant Classification
- **Missense_Mutation**: Cambia aminoácido
- **Silent**: No cambia aminoácido  
- **Nonsense_Mutation**: Crea codón de parada
- **Frame_Shift**: Cambia marco de lectura

### TMB Plot
- **Altura**: Número de mutaciones por muestra
- **Colores**: Diferentes tipos de mutación
- **Línea roja**: Mediana de todas las muestras

### Top Genes
- Muestra los genes con más mutaciones en el dataset
- Útil para identificar posibles genes driver

## Ejemplo Completo

```python
import pandas as pd
from pyMut import PyMutation

# Cargar datos de ejemplo
data = pd.read_csv("src/pyMut/data/examples/tcga_laml_converted.tsv", sep='\t')

# Crear análisis
py_mut = PyMutation(data)

# Configurar alta calidad
PyMutation.configure_high_quality_plots()

# Generar summary plot personalizado
summary_fig = py_mut.summary_plot(
    title="TCGA-LAML Mutation Analysis",
    figsize=(18, 14),
    max_samples=150,
    top_genes_count=15
)

# Guardar en múltiples formatos
summary_fig.savefig("tcga_laml_summary.png")
summary_fig.savefig("tcga_laml_summary.pdf", dpi=300)

print("✅ Summary plot generado exitosamente!")
```

## Solución de Problemas

### Datos vacíos o gráficos inesperados
```python
# Verificar estructura de datos
print(data.head())
print(data.columns.tolist())
print(data['Variant_Classification'].value_counts())
```

### Demasiadas muestras en TMB plot
```python
# Limitar número de muestras mostradas
fig = py_mut.summary_plot(max_samples=50)
```

### Imágenes de baja calidad
```python
# Siempre configurar alta calidad al inicio
PyMutation.configure_high_quality_plots()
``` 