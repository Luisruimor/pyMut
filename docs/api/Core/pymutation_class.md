# PyMutation - Clase Principal para Análisis de Mutaciones

La clase **PyMutation** es el objeto central de pyMut que encapsula datos de mutaciones y proporciona métodos para análisis y visualización.

## ¿Qué es PyMutation?

PyMutation es la clase principal que representa un conjunto de datos de mutaciones. Contiene los datos en formato estructurado, metadatos sobre el origen de los datos, y métodos para realizar análisis y generar visualizaciones.

## Estructura del Objeto

### Atributos Principales

```python
class PyMutation:
    def __init__(self, data: pd.DataFrame, metadata: MutationMetadata, samples: List[str]):
        self.data = data           # DataFrame con mutaciones en formato VCF-like
        self.samples = samples     # Lista de muestras en el dataset
        self.metadata = metadata   # Metadatos sobre el origen y configuración
```

### Atributo `data` (pd.DataFrame)
Contiene las mutaciones en formato VCF-like con columnas estándar:

```
CHROM | POS | ID | REF | ALT | QUAL | FILTER | SAMPLE_001 | SAMPLE_002 | Hugo_Symbol | Variant_Classification
chr1  | 100 | .  | A   | G   | .    | .      | A|G        | A|A        | GENE1       | Missense_Mutation
chr2  | 200 | .  | C   | T   | .    | .      | C|C        | C|T        | GENE2       | Nonsense_Mutation
```

### Atributo `samples` (List[str])
Lista de identificadores de muestras:
```python
['SAMPLE_001', 'SAMPLE_002', 'SAMPLE_003', ...]
```

### Atributo `metadata` (MutationMetadata)
Información sobre el origen y configuración:
```python
metadata.source_format    # "MAF", "VCF", etc.
metadata.file_path       # Ruta del archivo original
metadata.loaded_at       # Timestamp de carga
metadata.filters         # Filtros aplicados
metadata.fasta          # Archivo FASTA de referencia
metadata.notes          # Comentarios del archivo original
```

## Creación de Objetos PyMutation

### Desde archivo MAF
```python
from pyMut.io import read_maf

# Método recomendado
py_mut = read_maf("mutations.maf")
```

### Creación manual (avanzado)
```python
import pandas as pd
from pyMut.core import PyMutation, MutationMetadata

# Crear DataFrame con formato requerido
data = pd.DataFrame({
    'CHROM': ['chr1', 'chr2'],
    'POS': [100, 200],
    'REF': ['A', 'C'],
    'ALT': ['G', 'T'],
    'SAMPLE_001': ['A|G', 'C|C'],
    'SAMPLE_002': ['A|A', 'C|T'],
    'Hugo_Symbol': ['GENE1', 'GENE2'],
    'Variant_Classification': ['Missense_Mutation', 'Nonsense_Mutation']
})

# Crear metadatos
metadata = MutationMetadata(
    source_format="Manual",
    file_path="manual_creation",
    filters=["."],
    fasta="",
    notes="Creado manualmente"
)

# Crear objeto PyMutation
samples = ['SAMPLE_001', 'SAMPLE_002']
py_mut = PyMutation(data, metadata, samples)
```

## Métodos de Visualización

### Summary Plot - Análisis Completo
```python
# Generar resumen completo con 6 visualizaciones
fig = py_mut.summary_plot(
    title="Mi Análisis de Mutaciones",
    figsize=(16, 12),
    max_samples=100,
    top_genes_count=15,
    show_interactive=True
)
fig.savefig("summary_analysis.png")
```

### Visualizaciones Individuales

#### Clasificación de Variantes
```python
fig = py_mut.variant_classification_plot(
    title="Distribución de Tipos de Mutación",
    figsize=(10, 6)
)
```

#### Tipos de Variante
```python
fig = py_mut.variant_type_plot(
    title="Distribución de Tipos de Variante",
    figsize=(10, 6)
)
```

#### Clases de SNV
```python
fig = py_mut.snv_class_plot(
    title="Distribución de Cambios Nucleotídicos",
    figsize=(10, 6)
)
```

#### Variantes por Muestra (TMB)
```python
fig = py_mut.variants_per_sample_plot(
    title="Carga Mutacional por Muestra",
    max_samples=50,
    figsize=(12, 6)
)
```

#### Resumen de Clasificación por Muestra
```python
fig = py_mut.variant_classification_summary_plot(
    title="Resumen de Clasificaciones",
    figsize=(10, 6)
)
```

#### Genes Más Mutados
```python
fig = py_mut.top_mutated_genes_plot(
    title="Top Genes Mutados",
    count=20,
    figsize=(10, 8)
)
```

## Métodos de Análisis

### Análisis TMB (Tumor Mutation Burden)
```python
# Análisis completo de carga mutacional
tmb_results = py_mut.calculate_tmb_analysis(
    variant_classification_column="Variant_Classification",
    genome_size_bp=60456963,  # WES
    output_dir="results",
    save_files=True
)

# Acceder a resultados
analysis_df = tmb_results['analysis']
statistics_df = tmb_results['statistics']
```

## Métodos de Utilidad

### Configuración de Alta Calidad
```python
# Configurar matplotlib para alta calidad (recomendado)
PyMutation.configure_high_quality_plots()

# Ahora todas las figuras se guardan automáticamente en alta calidad
fig = py_mut.summary_plot()
fig.savefig("high_quality_plot.png")  # Automáticamente DPI=300
```

### Guardado de Figuras
```python
# Método centralizado para guardar figuras
fig = py_mut.summary_plot()

# Guardado básico
py_mut.save_figure(fig, "analysis.png")

# Guardado personalizado
py_mut.save_figure(fig, "analysis.pdf", dpi=600, bbox_inches='tight')
```

## Ejemplo Completo de Uso

```python
from pyMut.io import read_maf
from pyMut import PyMutation
import matplotlib.pyplot as plt

# 1. CARGAR DATOS
print("📂 Cargando datos MAF...")
py_mut = read_maf("src/pyMut/data/examples/tcga_laml.maf.gz")

# 2. EXPLORAR ESTRUCTURA
print(f"✅ Datos cargados exitosamente:")
print(f"   • Muestras: {len(py_mut.samples)}")
print(f"   • Mutaciones: {len(py_mut.data)}")
print(f"   • Formato fuente: {py_mut.metadata.source_format}")
print(f"   • Archivo: {py_mut.metadata.file_path}")

# 3. CONFIGURAR ALTA CALIDAD
PyMutation.configure_high_quality_plots()

# 4. ANÁLISIS VISUAL COMPLETO
print("\n📊 Generando análisis visual...")
summary_fig = py_mut.summary_plot(
    title="TCGA-LAML: Análisis Completo de Mutaciones",
    figsize=(18, 14),
    max_samples=150,
    top_genes_count=20
)

# Guardar en múltiples formatos
py_mut.save_figure(summary_fig, "tcga_laml_summary.png")
py_mut.save_figure(summary_fig, "tcga_laml_summary.pdf", dpi=300)

# 5. ANÁLISIS TMB
print("\n🧬 Calculando TMB...")
tmb_results = py_mut.calculate_tmb_analysis(
    genome_size_bp=60456963,  # WES estándar
    output_dir="results/tmb_analysis",
    save_files=True
)

# Explorar resultados TMB
analysis_df = tmb_results['analysis']
statistics_df = tmb_results['statistics']

print(f"   • Muestras analizadas: {len(analysis_df)}")
print(f"   • TMB promedio: {analysis_df['TMB_Total_Normalized'].mean():.3f} mut/Mb")
print(f"   • TMB mediano: {analysis_df['TMB_Total_Normalized'].median():.3f} mut/Mb")

# 6. VISUALIZACIONES INDIVIDUALES
print("\n📈 Generando visualizaciones específicas...")

# TMB por muestra
tmb_fig = py_mut.variants_per_sample_plot(
    title="Carga Mutacional Tumoral (TMB)",
    max_samples=100
)
py_mut.save_figure(tmb_fig, "tmb_per_sample.png")

# Genes más mutados
genes_fig = py_mut.top_mutated_genes_plot(
    title="Top 25 Genes Más Mutados",
    count=25
)
py_mut.save_figure(genes_fig, "top_mutated_genes.png")

# Tipos de mutación
classification_fig = py_mut.variant_classification_plot(
    title="Distribución de Tipos de Mutación"
)
py_mut.save_figure(classification_fig, "variant_classification.png")

# 7. ANÁLISIS PERSONALIZADO
print("\n🔍 Análisis personalizado...")

# Identificar muestras con TMB alto
high_tmb_threshold = 10  # mut/Mb
high_tmb_samples = analysis_df[
    analysis_df['TMB_Total_Normalized'] > high_tmb_threshold
]

print(f"   • Muestras con TMB alto (>{high_tmb_threshold} mut/Mb): {len(high_tmb_samples)}")

if len(high_tmb_samples) > 0:
    print("   • Top 5 muestras con mayor TMB:")
    top_samples = high_tmb_samples.nlargest(5, 'TMB_Total_Normalized')
    for _, row in top_samples.iterrows():
        print(f"     - {row['Sample']}: {row['TMB_Total_Normalized']:.3f} mut/Mb")

# Genes más frecuentemente mutados
gene_counts = py_mut.data['Hugo_Symbol'].value_counts().head(10)
print(f"\n   • Top 10 genes más mutados:")
for gene, count in gene_counts.items():
    print(f"     - {gene}: {count} mutaciones")

print("\n✅ Análisis completo finalizado!")
print("📁 Archivos generados:")
print("   • tcga_laml_summary.png/pdf - Resumen completo")
print("   • tmb_per_sample.png - TMB por muestra")
print("   • top_mutated_genes.png - Genes más mutados")
print("   • variant_classification.png - Tipos de mutación")
print("   • results/tmb_analysis/ - Análisis TMB detallado")
```

## Acceso a Datos

### Explorar estructura de datos
```python
# Información básica
print(f"Dimensiones: {py_mut.data.shape}")
print(f"Columnas: {list(py_mut.data.columns)}")
print(f"Muestras: {py_mut.samples}")

# Primeras filas
print(py_mut.data.head())

# Información de columnas
print(py_mut.data.info())
```

### Filtrar datos
```python
# Filtrar por gen específico
gene_data = py_mut.data[py_mut.data['Hugo_Symbol'] == 'TP53']

# Filtrar por tipo de mutación
missense_data = py_mut.data[
    py_mut.data['Variant_Classification'] == 'Missense_Mutation'
]

# Filtrar por cromosoma
chr1_data = py_mut.data[py_mut.data['CHROM'] == 'chr1']
```

### Estadísticas descriptivas
```python
# Conteos por tipo de mutación
mutation_counts = py_mut.data['Variant_Classification'].value_counts()
print(mutation_counts)

# Genes únicos
unique_genes = py_mut.data['Hugo_Symbol'].nunique()
print(f"Genes únicos: {unique_genes}")

# Mutaciones por muestra
mutations_per_sample = {}
for sample in py_mut.samples:
    # Contar genotipos que no son REF|REF
    sample_mutations = 0
    for idx, row in py_mut.data.iterrows():
        genotype = row[sample]
        ref_allele = row['REF']
        if genotype != f"{ref_allele}|{ref_allele}":
            sample_mutations += 1
    mutations_per_sample[sample] = sample_mutations

print("Mutaciones por muestra:")
for sample, count in sorted(mutations_per_sample.items(), 
                           key=lambda x: x[1], reverse=True)[:10]:
    print(f"  {sample}: {count}")
```

## Integración con Otros Análisis

### Exportar datos para análisis externos
```python
# Exportar datos completos
py_mut.data.to_csv("mutations_export.tsv", sep='\t', index=False)

# Exportar solo mutaciones de una muestra específica
sample_name = py_mut.samples[0]
sample_mutations = py_mut.data[py_mut.data[sample_name] != f"{py_mut.data['REF']}|{py_mut.data['REF']}"]
sample_mutations.to_csv(f"{sample_name}_mutations.tsv", sep='\t', index=False)

# Exportar metadatos
metadata_info = {
    'source_format': py_mut.metadata.source_format,
    'file_path': py_mut.metadata.file_path,
    'loaded_at': str(py_mut.metadata.loaded_at),
    'total_samples': len(py_mut.samples),
    'total_mutations': len(py_mut.data)
}

import json
with open("metadata.json", "w") as f:
    json.dump(metadata_info, f, indent=2)
```

### Combinar con pandas para análisis avanzados
```python
import pandas as pd
import numpy as np

# Análisis de co-ocurrencia de mutaciones
def analyze_gene_cooccurrence(py_mut, genes_of_interest):
    """Analizar co-ocurrencia de mutaciones en genes específicos"""
    results = []
    
    for sample in py_mut.samples:
        sample_genes = []
        for gene in genes_of_interest:
            gene_data = py_mut.data[py_mut.data['Hugo_Symbol'] == gene]
            if not gene_data.empty:
                # Verificar si hay mutación en esta muestra
                has_mutation = False
                for idx, row in gene_data.iterrows():
                    genotype = row[sample]
                    ref_allele = row['REF']
                    if genotype != f"{ref_allele}|{ref_allele}":
                        has_mutation = True
                        break
                if has_mutation:
                    sample_genes.append(gene)
        
        results.append({
            'Sample': sample,
            'Mutated_Genes': sample_genes,
            'Gene_Count': len(sample_genes)
        })
    
    return pd.DataFrame(results)

# Ejemplo de uso
oncogenes = ['TP53', 'KRAS', 'PIK3CA', 'APC', 'EGFR']
cooccurrence_df = analyze_gene_cooccurrence(py_mut, oncogenes)
print(cooccurrence_df.head())
```

## Mejores Prácticas

### 1. Configuración inicial recomendada
```python
from pyMut.io import read_maf
from pyMut import PyMutation

# Siempre configurar alta calidad al inicio
PyMutation.configure_high_quality_plots()

# Cargar datos
py_mut = read_maf("data.maf")
```

### 2. Verificación de datos
```python
# Verificar carga exitosa
assert len(py_mut.samples) > 0, "No se encontraron muestras"
assert len(py_mut.data) > 0, "No se encontraron mutaciones"
assert 'Hugo_Symbol' in py_mut.data.columns, "Falta columna de genes"
```

### 3. Análisis sistemático
```python
# 1. Resumen visual
summary_fig = py_mut.summary_plot()

# 2. Análisis TMB
tmb_results = py_mut.calculate_tmb_analysis()

# 3. Visualizaciones específicas según necesidad
# 4. Análisis personalizado con los datos
```

### 4. Gestión de archivos
```python
import os

# Crear directorio de resultados
os.makedirs("results", exist_ok=True)

# Guardar todas las figuras en el directorio
py_mut.save_figure(summary_fig, "results/summary.png")
```