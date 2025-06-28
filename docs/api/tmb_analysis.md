# TMB Analysis - Análisis de Carga Mutacional Tumoral

El método **calculate_tmb_analysis** permite calcular la Carga Mutacional Tumoral (TMB - Tumor Mutation Burden) para cada muestra en un objeto PyMutation.

## ¿Qué es el Análisis TMB?

El TMB es una métrica que cuantifica el número total de mutaciones en una muestra tumoral, normalizado por el tamaño del genoma interrogado. Es un biomarcador importante en oncología para predecir la respuesta a inmunoterapia.

## Características Principales

- **Cálculo automático por muestra**: Analiza cada muestra individualmente
- **Detección automática de columnas**: Identifica automáticamente las columnas de clasificación de variantes
- **Mutaciones no sinónimas**: Distingue entre mutaciones totales y no sinónimas
- **Normalización por tamaño de genoma**: Calcula TMB por millón de bases
- **Estadísticas globales**: Genera estadísticas descriptivas del conjunto de datos
- **Exportación automática**: Guarda resultados en archivos TSV

## Uso Básico

```python
from pyMut.io import read_maf

# Cargar datos
py_mut = read_maf("mutations.maf")

# Análisis TMB básico
tmb_results = py_mut.calculate_tmb_analysis()

# Acceder a resultados
analysis_df = tmb_results['analysis']
statistics_df = tmb_results['statistics']

print(f"Muestras analizadas: {len(analysis_df)}")
print(f"TMB promedio: {analysis_df['TMB_Total_Normalized'].mean():.2f} mut/Mb")
```

## Parámetros

### variant_classification_column (str, opcional)
- **Descripción**: Nombre de la columna con clasificación de variantes
- **Detección automática**: Si es None, detecta automáticamente columnas como 'Variant_Classification'
- **Ejemplo**: `"Variant_Classification"` o `"gencode_19_variant_classification"`

### genome_size_bp (int, default=60456963)
- **Descripción**: Tamaño de la región interrogada en pares de bases
- **WES (default)**: 60,456,963 bp (exoma completo)
- **WGS**: ~3,000,000,000 bp (genoma completo)
- **Paneles**: Tamaño específico del panel utilizado

### output_dir (str, default=".")
- **Descripción**: Directorio donde guardar los archivos de resultados
- **Archivos generados**: `TMB_analysis.tsv` y `TMB_statistics.tsv`

### save_files (bool, default=True)
- **Descripción**: Si guardar los resultados en archivos TSV
- **False**: Solo retorna los DataFrames sin guardar archivos

## Valor de Retorno

Retorna un diccionario con dos DataFrames:

```python
{
    'analysis': DataFrame,      # Análisis por muestra
    'statistics': DataFrame     # Estadísticas globales
}
```

### DataFrame 'analysis'
```
Sample | Total_Mutations | Non_Synonymous_Mutations | TMB_Total_Normalized | TMB_Non_Synonymous_Normalized
SAMPLE_001 | 45 | 32 | 0.744 | 0.529
SAMPLE_002 | 123 | 89 | 2.034 | 1.472
```

### DataFrame 'statistics'
```
Metric | Count | Mean | Median | Min | Max | Q1 | Q3 | Std
Total_Mutations | 100 | 67.5 | 45.0 | 12 | 234 | 32.0 | 89.0 | 45.2
TMB_Total_Normalized | 100 | 1.117 | 0.744 | 0.198 | 3.871 | 0.529 | 1.472 | 0.748
```

## Tipos de Mutaciones No Sinónimas

El análisis identifica automáticamente mutaciones con impacto biológico:

```python
# Tipos considerados no sinónimos
non_synonymous_types = {
    'MISSENSE_MUTATION',        # Cambio de aminoácido
    'NONSENSE_MUTATION',        # Codón de parada prematuro
    'FRAME_SHIFT_DEL',          # Deleción que cambia marco
    'FRAME_SHIFT_INS',          # Inserción que cambia marco
    'NONSTOP_MUTATION',         # Pérdida de codón de parada
    'TRANSLATION_START_SITE',   # Afecta sitio de inicio
    'SPLICE_SITE',              # Afecta sitio de splicing
    'IN_FRAME_DEL',             # Deleción en marco
    'IN_FRAME_INS',             # Inserción en marco
    'START_CODON_SNP',          # SNP en codón de inicio
    'START_CODON_DEL',          # Deleción en codón de inicio
    'START_CODON_INS',          # Inserción en codón de inicio
    'STOP_CODON_DEL',           # Deleción en codón de parada
    'STOP_CODON_INS'            # Inserción en codón de parada
}
```

## Ejemplo Completo

```python
from pyMut.io import read_maf
import pandas as pd

# Cargar datos MAF
py_mut = read_maf("src/pyMut/data/examples/tcga_laml.maf.gz")

# Análisis TMB personalizado
tmb_results = py_mut.calculate_tmb_analysis(
    variant_classification_column="Variant_Classification",
    genome_size_bp=60456963,  # WES estándar
    output_dir="results/tmb_analysis",
    save_files=True
)

# Explorar resultados
analysis_df = tmb_results['analysis']
statistics_df = tmb_results['statistics']

# Estadísticas básicas
print("=== RESUMEN TMB ===")
print(f"Muestras analizadas: {len(analysis_df)}")
print(f"TMB promedio: {analysis_df['TMB_Total_Normalized'].mean():.3f} mut/Mb")
print(f"TMB mediano: {analysis_df['TMB_Total_Normalized'].median():.3f} mut/Mb")

# Identificar muestras con TMB alto (>10 mut/Mb)
high_tmb = analysis_df[analysis_df['TMB_Total_Normalized'] > 10]
print(f"Muestras con TMB alto (>10 mut/Mb): {len(high_tmb)}")

# Top 5 muestras con mayor TMB
top_tmb = analysis_df.nlargest(5, 'TMB_Total_Normalized')
print("\nTop 5 muestras con mayor TMB:")
for _, row in top_tmb.iterrows():
    print(f"  {row['Sample']}: {row['TMB_Total_Normalized']:.3f} mut/Mb")

# Visualización opcional
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
plt.hist(analysis_df['TMB_Total_Normalized'], bins=30, alpha=0.7, edgecolor='black')
plt.xlabel('TMB (mutaciones/Mb)')
plt.ylabel('Número de muestras')
plt.title('Distribución de TMB en el conjunto de datos')
plt.axvline(analysis_df['TMB_Total_Normalized'].median(), color='red', 
           linestyle='--', label=f'Mediana: {analysis_df["TMB_Total_Normalized"].median():.3f}')
plt.legend()
plt.savefig('results/tmb_distribution.png', dpi=300, bbox_inches='tight')
plt.show()

print("✅ Análisis TMB completado exitosamente!")
```

## Configuración por Tipo de Secuenciación

### Whole Exome Sequencing (WES)
```python
tmb_results = py_mut.calculate_tmb_analysis(
    genome_size_bp=60456963  # ~60.5 Mb (default)
)
```

### Whole Genome Sequencing (WGS)
```python
tmb_results = py_mut.calculate_tmb_analysis(
    genome_size_bp=3000000000  # ~3 Gb
)
```

### Panel de Genes Específico
```python
# Ejemplo: Panel de 500 genes (~1.5 Mb)
tmb_results = py_mut.calculate_tmb_analysis(
    genome_size_bp=1500000  # 1.5 Mb
)
```

## Interpretación de Resultados

### Valores TMB Típicos
- **TMB bajo**: < 6 mutaciones/Mb
- **TMB intermedio**: 6-20 mutaciones/Mb  
- **TMB alto**: > 20 mutaciones/Mb

### Contexto Clínico
- **TMB alto**: Asociado con mejor respuesta a inmunoterapia
- **Tipos de cáncer**: Melanoma y cáncer de pulmón suelen tener TMB más alto
- **Microsatellite instability (MSI)**: Tumores MSI-H típicamente tienen TMB muy alto

## Archivos de Salida

### TMB_analysis.tsv
Contiene el análisis detallado por muestra:
```
Sample	Total_Mutations	Non_Synonymous_Mutations	TMB_Total_Normalized	TMB_Non_Synonymous_Normalized
TCGA-AB-2802	45	32	0.744063	0.529244
TCGA-AB-2803	123	89	2.034127	1.472063
```

### TMB_statistics.tsv
Contiene estadísticas descriptivas globales:
```
Metric	Count	Mean	Median	Min	Max	Q1	Q3	Std
Total_Mutations	191	67.534031	45.000000	12	234	32.000000	89.000000	45.234567
TMB_Total_Normalized	191	1.117234	0.744063	0.198456	3.871234	0.529244	1.472063	0.748123
```

## Validación y Control de Calidad

### Verificar detección de columnas
```python
# Verificar qué columna se detectó automáticamente
tmb_results = py_mut.calculate_tmb_analysis()
# Revisar logs para ver: "Auto-detected variant classification column: ..."
```

### Validar resultados
```python
analysis_df = tmb_results['analysis']

# Verificar que no hay valores negativos
assert (analysis_df['Total_Mutations'] >= 0).all()
assert (analysis_df['TMB_Total_Normalized'] >= 0).all()

# Verificar coherencia
assert (analysis_df['Non_Synonymous_Mutations'] <= analysis_df['Total_Mutations']).all()
```

## Manejo de Errores

### PyMutation inválido
```python
try:
    tmb_results = py_mut.calculate_tmb_analysis()
except ValueError as e:
    print(f"❌ Error en objeto PyMutation: {e}")
```

### Columna de clasificación no encontrada
```python
try:
    tmb_results = py_mut.calculate_tmb_analysis(
        variant_classification_column="columna_inexistente"
    )
except ValueError as e:
    print(f"❌ Columna no encontrada: {e}")
```

### Sin muestras válidas
```python
try:
    tmb_results = py_mut.calculate_tmb_analysis()
except ValueError as e:
    print(f"❌ Sin muestras válidas: {e}")
```

## Optimización para Datasets Grandes

```python
# Para datasets muy grandes, desactivar guardado de archivos
# y procesar en lotes si es necesario
tmb_results = py_mut.calculate_tmb_analysis(
    save_files=False  # Solo retorna DataFrames
)

# Procesar resultados en memoria
analysis_df = tmb_results['analysis']
# ... análisis personalizado ...

# Guardar manualmente si es necesario
analysis_df.to_csv("custom_tmb_analysis.tsv", sep='\t', index=False)
```