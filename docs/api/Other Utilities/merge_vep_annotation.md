# merge_vep_annotation - Fusión de Anotaciones VEP

La función **merge_maf_with_vep_annotations** permite combinar archivos MAF con anotaciones VEP (Variant Effect Predictor), enriqueciendo los datos de mutaciones con información funcional detallada.

## ¿Qué es merge_maf_with_vep_annotations?

Es una función que fusiona datos de mutaciones en formato MAF con anotaciones funcionales generadas por VEP, utilizando pandas y DuckDB para optimización. Esto permite añadir información detallada sobre el impacto funcional de las mutaciones.

## Características Principales

- **Fusión inteligente**: Mapea mutaciones MAF con anotaciones VEP usando coordenadas genómicas
- **Optimización con DuckDB**: Utiliza DuckDB para operaciones de fusión eficientes
- **Parsing automático**: Extrae y estructura información del campo "Extra" de VEP
- **Filtrado de calidad**: Elimina anotaciones de bajo impacto automáticamente
- **Compresión opcional**: Puede generar archivos de salida comprimidos
- **Preservación de datos**: Mantiene toda la información original del MAF

## Uso Básico

```python
from pyMut.utils.merge_vep_annotation import merge_maf_with_vep_annotations

# Fusionar MAF con anotaciones VEP
merged_df, output_path = merge_maf_with_vep_annotations(
    maf_file="mutations.maf",
    vep_file="vep_annotations.txt",
    output_file="mutations_vep_annotated.maf",
    compress=False
)

print(f"Archivo fusionado guardado en: {output_path}")
print(f"Mutaciones anotadas: {len(merged_df)}")
```

## Parámetros

### maf_file (str | Path) [requerido]
- **Descripción**: Ruta al archivo MAF original
- **Formatos soportados**: `.maf`, `.maf.gz`
- **Ejemplo**: `"data/tcga_mutations.maf.gz"`

### vep_file (str | Path) [requerido]
- **Descripción**: Ruta al archivo de anotaciones VEP
- **Formato**: Archivo de texto con formato VEP estándar (`.txt`)
- **Ejemplo**: `"annotations/vep_output.txt"`

### output_file (str | Path, opcional)
- **Descripción**: Ruta del archivo de salida fusionado
- **Default**: Se genera automáticamente con sufijo "_VEP_annotated"
- **Ejemplo**: `"results/mutations_annotated.maf"`

### compress (bool, default=False)
- **Descripción**: Si comprimir el archivo de salida con gzip
- **True**: Genera archivo `.maf.gz`
- **False**: Genera archivo `.maf` sin comprimir

## Valor de Retorno

Retorna una **tupla** con:
- **DataFrame**: Datos fusionados con anotaciones VEP
- **Path**: Ruta del archivo de salida creado

## Proceso de Fusión

### 1. Lectura de Archivos

```python
# El proceso interno incluye:
# - Lectura del archivo MAF (maneja compresión automáticamente)
# - Lectura del archivo VEP (detecta header automáticamente)
# - Creación de claves de región para mapeo
```

### 2. Creación de Claves de Mapeo

```python
# Para MAF: chr:pos-pos:1/alt
# Ejemplo: "chr17:7577121-7577121:1/T"

# Para VEP: usa el campo "Uploaded_variation"
# Debe coincidir con el formato de clave MAF
```

### 3. Parsing de Anotaciones VEP

```python
# El campo "Extra" de VEP se parsea automáticamente:
# "IMPACT=HIGH;SIFT=deleterious;PolyPhen=probably_damaging"
# Se convierte en columnas separadas:
# - IMPACT: HIGH
# - SIFT: deleterious  
# - PolyPhen: probably_damaging
```

## Ejemplos Detallados

### Fusión Básica

```python
from pyMut.utils.merge_vep_annotation import merge_maf_with_vep_annotations
import logging

# Configurar logging para ver progreso
logging.basicConfig(level=logging.INFO)

# Fusión con configuración básica
merged_df, output_path = merge_maf_with_vep_annotations(
    maf_file="src/pyMut/data/examples/tcga_laml.maf.gz",
    vep_file="annotations/tcga_laml_vep.txt"
)

print(f"=== Resultados de Fusión ===")
print(f"Archivo de salida: {output_path}")
print(f"Mutaciones originales: {len(merged_df)}")

# Verificar nuevas columnas VEP
vep_columns = [col for col in merged_df.columns if any(
    term in col.upper() for term in ['VEP', 'SIFT', 'POLYPHEN', 'IMPACT']
)]
print(f"Columnas VEP añadidas: {len(vep_columns)}")
print(f"Ejemplos: {vep_columns[:5]}")
```

### Fusión con Compresión

```python
# Generar archivo comprimido para ahorrar espacio
merged_df, output_path = merge_maf_with_vep_annotations(
    maf_file="large_mutations.maf.gz",
    vep_file="large_vep_annotations.txt",
    output_file="results/annotated_mutations.maf",
    compress=True  # Genera .maf.gz
)

print(f"Archivo comprimido: {output_path}")

# Verificar tamaño del archivo
import os
file_size_mb = os.path.getsize(output_path) / (1024 * 1024)
print(f"Tamaño del archivo: {file_size_mb:.2f} MB")
```

### Pipeline de Anotación Completo

```python
def pipeline_anotacion_completo(maf_input, vep_input, output_dir):
    """
    Pipeline completo de anotación con VEP
    """
    from pathlib import Path
    import pandas as pd
    
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    print("=== Pipeline de Anotación VEP ===")
    
    # 1. Fusionar con VEP
    print("1. Fusionando MAF con anotaciones VEP...")
    merged_df, annotated_file = merge_maf_with_vep_annotations(
        maf_file=maf_input,
        vep_file=vep_input,
        output_file=output_dir / "annotated.maf",
        compress=True
    )
    
    # 2. Análisis de calidad de anotaciones
    print("2. Analizando calidad de anotaciones...")
    
    # Contar mutaciones con diferentes tipos de impacto
    if 'IMPACT' in merged_df.columns:
        impact_counts = merged_df['IMPACT'].value_counts()
        print("Distribución de impacto:")
        for impact, count in impact_counts.items():
            print(f"  {impact}: {count}")
    
    # 3. Filtrar mutaciones de alto impacto
    print("3. Filtrando mutaciones de alto impacto...")
    if 'IMPACT' in merged_df.columns:
        high_impact = merged_df[merged_df['IMPACT'].isin(['HIGH', 'MODERATE'])]
        high_impact_file = output_dir / "high_impact_mutations.maf"
        high_impact.to_csv(high_impact_file, sep='\t', index=False)
        print(f"Mutaciones de alto impacto: {len(high_impact)} -> {high_impact_file}")
    
    # 4. Generar resumen
    print("4. Generando resumen...")
    summary = {
        'total_mutations': len(merged_df),
        'annotated_mutations': merged_df['Gene'].notna().sum(),
        'high_impact': len(high_impact) if 'IMPACT' in merged_df.columns else 0,
        'output_file': str(annotated_file)
    }
    
    summary_file = output_dir / "annotation_summary.txt"
    with open(summary_file, 'w') as f:
        for key, value in summary.items():
            f.write(f"{key}: {value}\n")
    
    print(f"Resumen guardado en: {summary_file}")
    return merged_df, summary

# Ejemplo de uso
resultado, resumen = pipeline_anotacion_completo(
    "mutations.maf",
    "vep_annotations.txt",
    "annotated_results/"
)
```

## Información de Anotaciones VEP

### Campos VEP Comunes

```python
# Después de la fusión, campos VEP típicos incluyen:
vep_fields = {
    'Gene': 'Símbolo del gen afectado',
    'Feature': 'ID del transcrito/feature',
    'Feature_type': 'Tipo de feature (Transcript, etc.)',
    'Consequence': 'Consecuencia de la mutación',
    'IMPACT': 'Nivel de impacto (HIGH, MODERATE, LOW, MODIFIER)',
    'SYMBOL': 'Símbolo oficial del gen',
    'SIFT': 'Predicción SIFT (deleterious/tolerated)',
    'PolyPhen': 'Predicción PolyPhen',
    'CADD_PHRED': 'Score CADD',
    'gnomAD_AF': 'Frecuencia alélica en gnomAD',
    'CLIN_SIG': 'Significado clínico'
}

# Verificar campos disponibles después de la fusión
merged_df, _ = merge_maf_with_vep_annotations("data.maf", "vep.txt")
available_fields = [field for field in vep_fields.keys() if field in merged_df.columns]
print(f"Campos VEP disponibles: {available_fields}")
```

### Filtrado por Impacto

```python
def filtrar_por_impacto(merged_df, impactos=['HIGH', 'MODERATE']):
    """
    Filtra mutaciones por nivel de impacto VEP
    """
    if 'IMPACT' not in merged_df.columns:
        print("⚠️ Campo IMPACT no disponible")
        return merged_df
    
    filtered = merged_df[merged_df['IMPACT'].isin(impactos)]
    
    print(f"=== Filtrado por Impacto ===")
    print(f"Impactos incluidos: {impactos}")
    print(f"Mutaciones originales: {len(merged_df)}")
    print(f"Mutaciones filtradas: {len(filtered)}")
    print(f"Porcentaje retenido: {len(filtered)/len(merged_df)*100:.1f}%")
    
    return filtered

# Ejemplo de uso
high_impact_mutations = filtrar_por_impacto(merged_df, ['HIGH'])
moderate_high_impact = filtrar_por_impacto(merged_df, ['HIGH', 'MODERATE'])
```

## Casos de Uso Comunes

### Enriquecimiento de Datos TCGA

```python
# Anotar datos TCGA con información funcional
tcga_annotated, output_file = merge_maf_with_vep_annotations(
    maf_file="TCGA-LAML.maf.gz",
    vep_file="TCGA-LAML_VEP.txt",
    compress=True
)

# Analizar genes con mutaciones de alto impacto
if 'IMPACT' in tcga_annotated.columns:
    high_impact_genes = tcga_annotated[
        tcga_annotated['IMPACT'] == 'HIGH'
    ]['Hugo_Symbol'].value_counts()
    
    print("Genes con más mutaciones de alto impacto:")
    print(high_impact_genes.head(10))
```

### Preparación para Análisis Funcional

```python
# Preparar datos para análisis de pathways
def preparar_analisis_funcional(merged_df):
    """
    Prepara datos anotados para análisis de pathways
    """
    # Filtrar mutaciones con información funcional
    functional_data = merged_df[
        (merged_df['IMPACT'].isin(['HIGH', 'MODERATE'])) &
        (merged_df['Gene'].notna()) &
        (merged_df['Consequence'].notna())
    ].copy()
    
    # Crear tabla de genes-muestras
    gene_sample_matrix = functional_data.pivot_table(
        index='Hugo_Symbol',
        columns='Tumor_Sample_Barcode',
        values='IMPACT',
        aggfunc='count',
        fill_value=0
    )
    
    # Guardar para análisis downstream
    gene_sample_matrix.to_csv("gene_sample_matrix.csv")
    
    # Resumen de consecuencias
    consequence_summary = functional_data.groupby(['Hugo_Symbol', 'Consequence']).size().reset_index(name='Count')
    consequence_summary.to_csv("consequence_summary.csv", index=False)
    
    print(f"Genes con mutaciones funcionales: {len(gene_sample_matrix)}")
    print(f"Muestras analizadas: {len(gene_sample_matrix.columns)}")
    
    return gene_sample_matrix, consequence_summary

# Ejemplo de uso
gene_matrix, consequences = preparar_analisis_funcional(merged_df)
```

### Integración con pyMut

```python
from pyMut.utils.merge_vep_annotation import merge_maf_with_vep_annotations
from pyMut.input import read_maf

# 1. Fusionar con VEP
merged_df, annotated_file = merge_maf_with_vep_annotations(
    maf_file="original.maf",
    vep_file="vep_annotations.txt",
    compress=True
)

# 2. Cargar en pyMut
py_mut_annotated = read_maf(annotated_file)

# 3. Análisis con datos enriquecidos
# Filtrar por genes de alto impacto
if 'IMPACT' in py_mut_annotated.data.columns:
    high_impact_data = py_mut_annotated.data[
        py_mut_annotated.data['IMPACT'] == 'HIGH'
    ]
    print(f"Mutaciones de alto impacto: {len(high_impact_data)}")

# 4. Análisis TMB con datos anotados
tmb_results = py_mut_annotated.calculate_tmb_analysis()
print(f"TMB calculado con {len(py_mut_annotated.samples)} muestras anotadas")
```

## Manejo de Errores

### Archivos VEP Malformados

```python
try:
    merged_df, output_path = merge_maf_with_vep_annotations(
        "data.maf",
        "malformed_vep.txt"
    )
except ValueError as e:
    print(f"❌ Error en formato VEP: {e}")
    # Posibles causas:
    # - Falta header #Uploaded_variation
    # - Formato de columnas incorrecto
    # - Archivo corrupto
```

### Problemas de Mapeo

```python
# Verificar éxito del mapeo
merged_df, _ = merge_maf_with_vep_annotations("data.maf", "vep.txt")

# Contar mutaciones mapeadas
original_count = len(merged_df)
mapped_count = merged_df['Gene'].notna().sum()
mapping_rate = mapped_count / original_count * 100

print(f"Tasa de mapeo: {mapping_rate:.1f}%")

if mapping_rate < 50:
    print("⚠️ Baja tasa de mapeo - verificar:")
    print("- Compatibilidad de coordenadas genómicas")
    print("- Formato de claves de región")
    print("- Versión del genoma de referencia")
```

### Archivos de Salida Grandes

```python
import os

# Monitorear tamaño del archivo de salida
def fusionar_con_monitoreo(maf_file, vep_file, max_size_gb=5):
    """
    Fusiona archivos monitoreando el tamaño de salida
    """
    merged_df, output_path = merge_maf_with_vep_annotations(
        maf_file, vep_file, compress=True
    )
    
    file_size_gb = os.path.getsize(output_path) / (1024**3)
    
    if file_size_gb > max_size_gb:
        print(f"⚠️ Archivo grande: {file_size_gb:.2f} GB")
        print("Considerar:")
        print("- Filtrar por impacto antes de guardar")
        print("- Usar compresión adicional")
        print("- Dividir en archivos más pequeños")
    
    return merged_df, output_path

# Ejemplo de uso
resultado = fusionar_con_monitoreo("large.maf", "large_vep.txt", max_size_gb=2)
```

## Optimización y Rendimiento

### Para Datasets Grandes

```python
# Optimización de memoria para datasets grandes
def fusionar_optimizado(maf_file, vep_file, chunk_size=10000):
    """
    Fusión optimizada para datasets muy grandes
    """
    print("Procesando en modo optimizado...")
    
    # La función ya usa DuckDB internamente para optimización
    # Para datasets extremadamente grandes, considerar:
    # 1. Filtrar MAF antes de la fusión
    # 2. Usar compresión
    # 3. Procesar por cromosomas
    
    merged_df, output_path = merge_maf_with_vep_annotations(
        maf_file, vep_file, compress=True
    )
    
    return merged_df, output_path
```

### Monitoreo de Recursos

```python
import psutil
import time

def fusionar_con_estadisticas(maf_file, vep_file):
    """
    Fusiona archivos monitoreando recursos del sistema
    """
    start_time = time.time()
    start_memory = psutil.virtual_memory().used / 1024**3  # GB
    
    print(f"Memoria inicial: {start_memory:.2f} GB")
    
    merged_df, output_path = merge_maf_with_vep_annotations(
        maf_file, vep_file
    )
    
    end_time = time.time()
    end_memory = psutil.virtual_memory().used / 1024**3  # GB
    
    print(f"Tiempo total: {end_time - start_time:.2f} segundos")
    print(f"Memoria final: {end_memory:.2f} GB")
    print(f"Incremento de memoria: {end_memory - start_memory:.2f} GB")
    print(f"Mutaciones procesadas: {len(merged_df)}")
    
    return merged_df, output_path

# Ejemplo de uso
resultado = fusionar_con_estadisticas("data.maf", "vep.txt")
```

## Solución de Problemas

### Claves de Región No Coinciden

```python
# Verificar formato de claves de región
def diagnosticar_mapeo(maf_file, vep_file):
    """
    Diagnostica problemas de mapeo entre MAF y VEP
    """
    import pandas as pd
    
    # Leer archivos por separado
    maf_df = pd.read_csv(maf_file, sep='\t', comment='#')
    
    # Crear claves de ejemplo
    sample_maf_keys = []
    for _, row in maf_df.head(5).iterrows():
        chrom = f"chr{row['Chromosome']}" if not str(row['Chromosome']).startswith('chr') else str(row['Chromosome'])
        start = int(row['Start_Position'])
        alt = str(row['Tumor_Seq_Allele2'])
        key = f"{chrom}:{start}-{start}:1/{alt}"
        sample_maf_keys.append(key)
    
    print("Claves MAF de ejemplo:")
    for key in sample_maf_keys:
        print(f"  {key}")
    
    # Leer VEP y mostrar claves de ejemplo
    with open(vep_file, 'r') as f:
        vep_lines = [line.strip() for line in f if not line.startswith('#')]
    
    print("\nClaves VEP de ejemplo:")
    for line in vep_lines[:5]:
        if line:
            fields = line.split('\t')
            if len(fields) > 0:
                print(f"  {fields[0]}")

# Ejemplo de uso
diagnosticar_mapeo("data.maf", "vep.txt")
```