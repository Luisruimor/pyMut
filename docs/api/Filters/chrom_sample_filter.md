# filter_by_chrom_sample - Filtro por Cromosoma y Muestra

El método **filter_by_chrom_sample** permite filtrar datos de PyMutation por cromosoma y/o muestra, proporcionando un control granular sobre qué datos incluir en el análisis.

## ¿Qué es filter_by_chrom_sample?

Es un método versátil que permite filtrar datos por cromosoma, muestra, o ambos criterios simultáneamente. Además de filtrar filas, también maneja el filtrado de columnas cuando se especifican muestras, manteniendo la integridad del formato de datos.

## Características Principales

- **Filtrado dual**: Por cromosoma y/o muestra en una sola operación
- **Filtrado de filas y columnas**: Elimina tanto filas como columnas irrelevantes
- **Soporte para múltiples valores**: Acepta listas de cromosomas y muestras
- **Compatibilidad MAF/VCF**: Maneja ambos formatos automáticamente
- **Preservación de metadatos**: Registra todos los filtros aplicados
- **Validación automática**: Verifica la existencia de cromosomas y muestras
- **Logging detallado**: Proporciona información sobre el proceso de filtrado

## Uso Básico

```python
from pyMut.input import read_maf

# Cargar datos
py_mut = read_maf("mutations.maf")

# Filtrar por cromosoma únicamente
chr17_data = py_mut.filter_by_chrom_sample(chrom="chr17")

# Filtrar por muestra únicamente
sample_data = py_mut.filter_by_chrom_sample(sample="TCGA-AB-2802")

# Filtrar por ambos criterios
chr17_sample = py_mut.filter_by_chrom_sample(
    chrom="chr17", 
    sample="TCGA-AB-2802"
)

print(f"Mutaciones en chr17: {len(chr17_data.data)}")
print(f"Mutaciones en muestra: {len(sample_data.data)}")
print(f"Mutaciones chr17 + muestra: {len(chr17_sample.data)}")
```

## Parámetros

### chrom (str, list, opcional)
- **Descripción**: Cromosoma(s) a filtrar
- **Formatos aceptados**: `"chr17"`, `"17"`, `["chr1", "chr17"]`, `["X", "Y"]`
- **Normalización**: Se normaliza automáticamente al formato estándar
- **Ejemplo**: `"chr17"` o `["chr1", "chr2", "chrX"]`

### sample (str, list, opcional)
- **Descripción**: Muestra(s) a filtrar
- **Formato**: Identificadores de muestra como aparecen en los datos
- **Ejemplo**: `"TCGA-AB-2802"` o `["TCGA-AB-2802", "TCGA-AB-2803"]`

### sample_column (str, opcional)
- **Descripción**: Nombre de la columna que contiene información de muestras
- **Default**: `"Tumor_Sample_Barcode"` (estándar MAF)
- **Uso**: Para datos con nombres de columna no estándar

## Comportamiento del Filtrado

### Filtrado Solo por Cromosoma
```python
# Mantiene todas las columnas, filtra solo filas
chr_filtered = py_mut.filter_by_chrom_sample(chrom="chr17")

# Resultado: Solo mutaciones en chr17, todas las muestras preservadas
print(f"Cromosomas únicos: {chr_filtered.data['CHROM'].unique()}")
print(f"Muestras preservadas: {len(chr_filtered.samples)}")
```

### Filtrado Solo por Muestra
```python
# Filtra tanto filas como columnas
sample_filtered = py_mut.filter_by_chrom_sample(sample="TCGA-AB-2802")

# Resultado: Solo mutaciones de la muestra, solo columnas relevantes
print(f"Muestras en datos: {sample_filtered.samples}")
print(f"Columnas de muestra: {[col for col in sample_filtered.data.columns if 'TCGA' in col]}")
```

### Filtrado Combinado
```python
# Aplica ambos filtros simultáneamente
combined = py_mut.filter_by_chrom_sample(
    chrom=["chr17", "chrX"], 
    sample=["TCGA-AB-2802", "TCGA-AB-2803"]
)

print(f"Cromosomas: {combined.data['CHROM'].unique()}")
print(f"Muestras: {combined.samples}")
```

## Ejemplos Detallados

### Análisis de Cromosomas Específicos

```python
from pyMut.input import read_maf
import logging

# Configurar logging
logging.basicConfig(level=logging.INFO)

# Cargar datos TCGA
py_mut = read_maf("src/pyMut/data/examples/tcga_laml.maf.gz")
print(f"Datos originales: {len(py_mut.data)} mutaciones, {len(py_mut.samples)} muestras")

# Análisis de cromosomas de interés oncológico
cromosomas_cancer = ["chr17", "chr13", "chr3", "chr7", "chr12"]

for chrom in cromosomas_cancer:
    chr_data = py_mut.filter_by_chrom_sample(chrom=chrom)
    
    print(f"\n=== Análisis {chrom} ===")
    print(f"Mutaciones: {len(chr_data.data)}")
    
    if len(chr_data.data) > 0:
        # Genes más mutados en este cromosoma
        top_genes = chr_data.data['Hugo_Symbol'].value_counts().head(3)
        print("Genes más mutados:")
        for gene, count in top_genes.items():
            print(f"  - {gene}: {count}")
        
        # Tipos de mutación más comunes
        top_types = chr_data.data['Variant_Classification'].value_counts().head(3)
        print("Tipos de mutación:")
        for mut_type, count in top_types.items():
            print(f"  - {mut_type}: {count}")
```

### Análisis de Muestras Específicas

```python
# Seleccionar muestras de interés
muestras_interes = ["TCGA-AB-2802", "TCGA-AB-2803", "TCGA-AB-2869"]

# Filtrar por muestras específicas
samples_data = py_mut.filter_by_chrom_sample(sample=muestras_interes)

print(f"=== Análisis de Muestras Específicas ===")
print(f"Muestras analizadas: {len(samples_data.samples)}")
print(f"Mutaciones totales: {len(samples_data.data)}")

# Análisis por muestra individual
for sample in samples_data.samples:
    sample_mutations = samples_data.data[
        samples_data.data['Tumor_Sample_Barcode'] == sample
    ]
    
    print(f"\n{sample}:")
    print(f"  Mutaciones: {len(sample_mutations)}")
    
    if len(sample_mutations) > 0:
        # Cromosomas afectados
        chroms = sample_mutations['CHROM'].value_counts()
        print(f"  Cromosomas afectados: {len(chroms)}")
        print(f"  Cromosoma más mutado: {chroms.index[0]} ({chroms.iloc[0]} mut)")
```

### Filtrado por Múltiples Criterios

```python
# Análisis de genes específicos en cromosomas de interés
def analizar_genes_cromosomas(py_mut, cromosomas, genes_interes):
    """
    Analiza genes específicos en cromosomas determinados
    """
    # Filtrar por cromosomas
    chr_data = py_mut.filter_by_chrom_sample(chrom=cromosomas)
    
    print(f"=== Análisis: {cromosomas} ===")
    print(f"Mutaciones en cromosomas: {len(chr_data.data)}")
    
    # Analizar cada gen de interés
    for gene in genes_interes:
        gene_mutations = chr_data.data[
            chr_data.data['Hugo_Symbol'] == gene
        ]
        
        if len(gene_mutations) > 0:
            print(f"\n{gene}:")
            print(f"  Mutaciones: {len(gene_mutations)}")
            
            # Muestras afectadas
            affected_samples = gene_mutations['Tumor_Sample_Barcode'].nunique()
            print(f"  Muestras afectadas: {affected_samples}")
            
            # Tipos de mutación
            mut_types = gene_mutations['Variant_Classification'].value_counts()
            print(f"  Tipo más común: {mut_types.index[0]} ({mut_types.iloc[0]})")

# Ejemplo de uso
cromosomas_tp53_brca = ["chr17", "chr13"]
genes_supresores = ["TP53", "BRCA1", "BRCA2", "RB1"]

analizar_genes_cromosomas(py_mut, cromosomas_tp53_brca, genes_supresores)
```

## Manejo de Formatos de Datos

### Datos Estilo MAF
```python
# Para datos MAF (formato largo con columna de muestra)
maf_data = read_maf("data.maf")

# Filtrado estándar por muestra
filtered = maf_data.filter_by_chrom_sample(
    sample=["TCGA-AB-2802", "TCGA-AB-2803"]
)

# El método detecta automáticamente la columna 'Tumor_Sample_Barcode'
print(f"Columna de muestra detectada: Tumor_Sample_Barcode")
```

### Datos Estilo VCF
```python
# Para datos VCF (formato ancho con columnas por muestra)
vcf_data = read_vcf("data.vcf.gz")

# Filtrado por muestra en formato VCF
filtered = vcf_data.filter_by_chrom_sample(
    sample=["SAMPLE_001", "SAMPLE_002"]
)

# El método maneja automáticamente el formato VCF
print(f"Muestras VCF filtradas: {filtered.samples}")
```

## Casos de Uso Avanzados

### Pipeline de Análisis Secuencial

```python
# Pipeline de filtrado progresivo
def pipeline_analisis_progresivo(py_mut):
    """
    Aplica filtros progresivos para análisis específico
    """
    print("=== Pipeline de Análisis Progresivo ===")
    
    # Paso 1: Filtrar cromosomas de interés
    step1 = py_mut.filter_by_chrom_sample(
        chrom=["chr17", "chr13", "chr3"]
    )
    print(f"Paso 1 - Cromosomas cáncer: {len(step1.data)} mutaciones")
    
    # Paso 2: Filtrar muestras con alta carga mutacional
    # (Ejemplo: muestras con >50 mutaciones en estos cromosomas)
    sample_counts = step1.data['Tumor_Sample_Barcode'].value_counts()
    high_burden_samples = sample_counts[sample_counts > 50].index.tolist()
    
    step2 = step1.filter_by_chrom_sample(sample=high_burden_samples)
    print(f"Paso 2 - Alta carga mutacional: {len(step2.data)} mutaciones")
    print(f"Muestras seleccionadas: {len(step2.samples)}")
    
    # Paso 3: Análisis final
    if len(step2.data) > 0:
        genes_frecuentes = step2.data['Hugo_Symbol'].value_counts().head(5)
        print("\nGenes más frecuentemente mutados:")
        for gene, count in genes_frecuentes.items():
            print(f"  {gene}: {count} mutaciones")
    
    return step2

# Ejecutar pipeline
resultado_final = pipeline_analisis_progresivo(py_mut)
```

### Análisis Comparativo de Subgrupos

```python
def comparar_subgrupos_cromosomicos(py_mut, grupo1_samples, grupo2_samples):
    """
    Compara patrones mutacionales entre dos grupos de muestras
    """
    # Filtrar cada grupo
    grupo1 = py_mut.filter_by_chrom_sample(sample=grupo1_samples)
    grupo2 = py_mut.filter_by_chrom_sample(sample=grupo2_samples)
    
    print("=== Comparación de Subgrupos ===")
    print(f"Grupo 1: {len(grupo1.data)} mutaciones en {len(grupo1.samples)} muestras")
    print(f"Grupo 2: {len(grupo2.data)} mutaciones en {len(grupo2.samples)} muestras")
    
    # Análisis por cromosoma
    cromosomas = ["chr1", "chr17", "chrX", "chrY"]
    
    for chrom in cromosomas:
        g1_chr = grupo1.filter_by_chrom_sample(chrom=chrom)
        g2_chr = grupo2.filter_by_chrom_sample(chrom=chrom)
        
        print(f"\n{chrom}:")
        print(f"  Grupo 1: {len(g1_chr.data)} mutaciones")
        print(f"  Grupo 2: {len(g2_chr.data)} mutaciones")
        
        # Calcular densidad mutacional
        if len(grupo1.samples) > 0 and len(grupo2.samples) > 0:
            densidad_g1 = len(g1_chr.data) / len(grupo1.samples)
            densidad_g2 = len(g2_chr.data) / len(grupo2.samples)
            print(f"  Densidad G1: {densidad_g1:.2f} mut/muestra")
            print(f"  Densidad G2: {densidad_g2:.2f} mut/muestra")

# Ejemplo de uso
grupo_a = ["TCGA-AB-2802", "TCGA-AB-2803", "TCGA-AB-2869"]
grupo_b = ["TCGA-AB-2988", "TCGA-AB-2989", "TCGA-AB-2990"]

comparar_subgrupos_cromosomicos(py_mut, grupo_a, grupo_b)
```

## Validación y Manejo de Errores

### Validación de Parámetros

```python
# Error: No se proporciona ningún parámetro
try:
    invalid = py_mut.filter_by_chrom_sample()
except ValueError as e:
    print(f"❌ Error: {e}")
    # "At least one of 'chrom' or 'sample' must be provided"

# Error: Cromosoma inexistente
try:
    invalid = py_mut.filter_by_chrom_sample(chrom="chr99")
except KeyError as e:
    print(f"❌ Cromosoma no encontrado: {e}")

# Error: Muestra inexistente
try:
    invalid = py_mut.filter_by_chrom_sample(sample="SAMPLE_INEXISTENTE")
except ValueError as e:
    print(f"❌ Muestra no encontrada: {e}")
```

### Manejo de Datos Vacíos

```python
# Filtro que no encuentra datos
empty_result = py_mut.filter_by_chrom_sample(
    chrom="chr17", 
    sample="MUESTRA_INEXISTENTE"
)

if len(empty_result.data) == 0:
    print("⚠️ El filtro no encontró datos")
    print("Verificar:")
    print("- Existencia del cromosoma en los datos")
    print("- Existencia de la muestra en los datos")
    print("- Compatibilidad entre cromosoma y muestra")
```

## Optimización y Rendimiento

### Para Datasets Grandes

```python
import time

# Medir rendimiento de filtrado
start_time = time.time()

# Filtrado eficiente: cromosoma primero (reduce datos)
chr_subset = py_mut.filter_by_chrom_sample(chrom=["chr17", "chr13"])
final_result = chr_subset.filter_by_chrom_sample(
    sample=["TCGA-AB-2802", "TCGA-AB-2803"]
)

end_time = time.time()

print(f"Filtrado completado en {end_time - start_time:.2f} segundos")
print(f"Datos finales: {len(final_result.data)} mutaciones")
```

### Monitoreo de Memoria

```python
import psutil
import os

def filtrar_con_monitoreo(py_mut, chrom, sample):
    """
    Filtra datos monitoreando el uso de memoria
    """
    process = psutil.Process(os.getpid())
    memoria_inicial = process.memory_info().rss / 1024 / 1024  # MB
    
    print(f"Memoria inicial: {memoria_inicial:.1f} MB")
    
    # Aplicar filtro
    resultado = py_mut.filter_by_chrom_sample(chrom=chrom, sample=sample)
    
    memoria_final = process.memory_info().rss / 1024 / 1024  # MB
    print(f"Memoria final: {memoria_final:.1f} MB")
    print(f"Incremento: {memoria_final - memoria_inicial:.1f} MB")
    
    return resultado

# Ejemplo de uso
resultado = filtrar_con_monitoreo(
    py_mut, 
    chrom=["chr17", "chr13"], 
    sample=["TCGA-AB-2802", "TCGA-AB-2803"]
)
```

## Integración con Otros Métodos

### Combinación con Filtros Genómicos

```python
# Combinar filtros cromosómicos con filtros de región
chr17_data = py_mut.filter_by_chrom_sample(chrom="chr17")
tp53_region = chr17_data.region("chr17", 7571720, 7590868)

print(f"Mutaciones en TP53: {len(tp53_region.data)}")
```

### Preparación para Análisis TMB

```python
# Filtrar datos para análisis TMB específico
samples_tmb = py_mut.filter_by_chrom_sample(
    sample=["TCGA-AB-2802", "TCGA-AB-2803", "TCGA-AB-2869"]
)

# Realizar análisis TMB en subconjunto
tmb_results = samples_tmb.calculate_tmb_analysis(
    genome_size_bp=60456963,  # WES
    output_dir="results_subset"
)

print(f"TMB calculado para {len(samples_tmb.samples)} muestras")
```