# Known Cancer Gene Annotation - Anotación de Genes de Cáncer Conocidos

La función **knownCancer** permite anotar archivos MAF con información específica de genes relacionados con cáncer utilizando datos de COSMIC Cancer Gene Census y opcionalmente OncoKB.

## ¿Qué es la Anotación de Genes de Cáncer Conocidos?

Esta función filtra y añade anotaciones específicas relacionadas con cáncer, enfocándose únicamente en campos relevantes para la investigación oncológica. Combina información de múltiples fuentes (COSMIC y OncoKB) y crea un campo unificado que indica si un gen es oncogénico según cualquier fuente.

## Características Principales

- **Anotación específica de cáncer**: Filtra solo campos relevantes para oncología
- **Múltiples fuentes**: Combina datos de COSMIC Cancer Gene Census y OncoKB
- **Campo unificado**: Crea `Is_Oncogene_any` que combina información de todas las fuentes
- **Detección automática**: Maneja automáticamente diferentes formatos de archivo
- **Optimización inteligente**: Usa DuckDB para archivos grandes (>2GB) y pandas para archivos pequeños
- **Manejo de sinónimos**: Resuelve automáticamente sinónimos de genes

## Uso Básico

```python
from pyMut.annotate.cosmic_cancer_annotate import knownCancer

# Anotación básica con COSMIC solamente
annotated_df, output_path = knownCancer(
    maf_file="mutations.maf",
    annotation_table="cosmic_cancer_gene_census.tsv"
)

print(f"Archivo anotado guardado en: {output_path}")
print(f"Genes anotados: {annotated_df['Is_Oncogene_any'].sum()}")
```

## Uso Avanzado con OncoKB

```python
# Anotación combinada con COSMIC y OncoKB
annotated_df, output_path = knownCancer(
    maf_file="mutations.maf.gz",
    annotation_table="cosmic_cancer_gene_census.tsv.gz",
    oncokb_table="oncokb_cancer_gene_list.tsv",
    output_path="annotated_mutations.maf.gz",
    compress_output=True
)

# Analizar resultados
cosmic_genes = annotated_df['COSMIC_ROLE_IN_CANCER'].notna().sum()
oncokb_genes = annotated_df['OncoKB_Is Oncogene'].eq('True').sum()
any_oncogene = annotated_df['Is_Oncogene_any'].sum()

print(f"Genes con anotación COSMIC: {cosmic_genes}")
print(f"Genes oncogénicos según OncoKB: {oncokb_genes}")
print(f"Genes oncogénicos según cualquier fuente: {any_oncogene}")
```

## Parámetros

### maf_file (str | Path, requerido)
- **Descripción**: Ruta al archivo MAF de entrada
- **Formatos soportados**: `.maf`, `.maf.gz`
- **Ejemplo**: `"mutations.maf"` o `"data/mutations.maf.gz"`

### annotation_table (str | Path, requerido)
- **Descripción**: Ruta a la tabla de anotación COSMIC Cancer Gene Census
- **Formatos soportados**: `.tsv`, `.tsv.gz`
- **Ejemplo**: `"cosmic_cancer_gene_census.tsv.gz"`

### output_path (str | Path, opcional)
- **Descripción**: Ruta del archivo de salida
- **Default**: Se genera automáticamente con sufijo `_knownCancer_annotated.maf`
- **Ejemplo**: `"annotated_mutations.maf.gz"`

### compress_output (bool, default=True)
- **Descripción**: Si comprimir el archivo de salida con gzip
- **True**: Genera archivo `.maf.gz`
- **False**: Genera archivo `.maf` sin comprimir

### join_column (str, default="Hugo_Symbol")
- **Descripción**: Columna para realizar el join entre archivos
- **Detección automática**: Utiliza aliases definidos en `fields.py`
- **Alternativas**: `"Gene_Symbol"`, `"GENE"`, etc.

### oncokb_table (str | Path, opcional)
- **Descripción**: Ruta a la tabla de OncoKB cancer gene list
- **Efecto**: Añade anotaciones OncoKB al resultado
- **Ejemplo**: `"oncokb_cancer_gene_list.tsv"`

## Valor de Retorno

Retorna una tupla con dos elementos:

```python
annotated_df, output_path = knownCancer(...)

# annotated_df: pandas.DataFrame con datos anotados
# output_path: pathlib.Path del archivo generado
```

### DataFrame Anotado

El DataFrame resultante contiene las columnas originales del MAF más las siguientes anotaciones específicas de cáncer:

#### Anotaciones COSMIC
- **COSMIC_ROLE_IN_CANCER**: Rol del gen en cáncer según COSMIC
- **COSMIC_TIER**: Nivel de evidencia en COSMIC

#### Anotaciones OncoKB (si se proporciona oncokb_table)
- **OncoKB_Is Oncogene**: Si el gen es oncogénico según OncoKB
- **OncoKB_Is Tumor Suppressor Gene**: Si es gen supresor tumoral
- **OncoKB_OncoKB Annotated**: Si está anotado en OncoKB
- **OncoKB_MSK-IMPACT**: Si está en el panel MSK-IMPACT
- **OncoKB_MSK-HEME**: Si está en el panel MSK-HEME
- **OncoKB_FOUNDATION ONE**: Si está en Foundation One
- **OncoKB_FOUNDATION ONE HEME**: Si está en Foundation One Heme
- **OncoKB_Vogelstein**: Si está en la lista de Vogelstein

#### Campo Unificado
- **Is_Oncogene_any**: `True` si el gen es oncogénico según cualquier fuente

## Ejemplos de Análisis

### Análisis de Genes Oncogénicos

```python
annotated_df, _ = knownCancer(
    maf_file="mutations.maf",
    annotation_table="cosmic.tsv",
    oncokb_table="oncokb.tsv"
)

# Genes oncogénicos por fuente
cosmic_oncogenes = annotated_df[
    annotated_df['COSMIC_ROLE_IN_CANCER'].notna()
]['Hugo_Symbol'].unique()

oncokb_oncogenes = annotated_df[
    annotated_df['OncoKB_Is Oncogene'] == 'True'
]['Hugo_Symbol'].unique()

any_oncogenes = annotated_df[
    annotated_df['Is_Oncogene_any'] == True
]['Hugo_Symbol'].unique()

print(f"Genes oncogénicos COSMIC: {len(cosmic_oncogenes)}")
print(f"Genes oncogénicos OncoKB: {len(oncokb_oncogenes)}")
print(f"Genes oncogénicos (cualquier fuente): {len(any_oncogenes)}")
```

### Filtrado por Tipo de Gen

```python
# Filtrar solo mutaciones en genes oncogénicos conocidos
oncogenic_mutations = annotated_df[
    annotated_df['Is_Oncogene_any'] == True
]

# Filtrar por genes supresores tumorales
tumor_suppressors = annotated_df[
    annotated_df['OncoKB_Is Tumor Suppressor Gene'] == 'True'
]

# Genes con evidencia de alto nivel (COSMIC Tier 1)
high_evidence_genes = annotated_df[
    annotated_df['COSMIC_TIER'] == '1'
]

print(f"Mutaciones en genes oncogénicos: {len(oncogenic_mutations)}")
print(f"Mutaciones en supresores tumorales: {len(tumor_suppressors)}")
print(f"Mutaciones en genes Tier 1: {len(high_evidence_genes)}")
```

### Análisis por Panel Clínico

```python
# Genes en paneles clínicos específicos
msk_impact_genes = annotated_df[
    annotated_df['OncoKB_MSK-IMPACT'] == 'True'
]['Hugo_Symbol'].unique()

foundation_one_genes = annotated_df[
    annotated_df['OncoKB_FOUNDATION ONE'] == 'True'
]['Hugo_Symbol'].unique()

print(f"Genes en MSK-IMPACT: {len(msk_impact_genes)}")
print(f"Genes en Foundation One: {len(foundation_one_genes)}")
```

## Archivos de Salida

La función genera automáticamente archivos con nombres descriptivos:

```python
# Ejemplos de nombres de archivo generados automáticamente
"mutations_knownCancer_annotated.maf.gz"           # Solo COSMIC
"mutations_knownCancer_annotated.maf"              # Sin compresión
"sample_data_knownCancer_annotated.maf.gz"         # Con OncoKB
```

## Optimización de Rendimiento

La función utiliza una estrategia híbrida para optimizar el rendimiento:

- **Archivos pequeños (<2GB)**: Utiliza pandas + pyarrow para máxima velocidad
- **Archivos grandes (≥2GB)**: Utiliza DuckDB para eficiencia de memoria
- **Detección automática**: El tamaño se evalúa automáticamente

## Manejo de Sinónimos

La función resuelve automáticamente sinónimos de genes:

```python
# Sinónimos manejados automáticamente
# COSMIC: columna "SYNONIMS"
# OncoKB: columna "Gene Aliases"

# Ejemplo: si el MAF contiene "ERBB2" y COSMIC tiene "HER2"
# La función los reconoce como el mismo gen
```

## Casos de Uso Comunes

1. **Priorización de variantes**: Enfocar análisis en genes oncogénicos conocidos
2. **Validación clínica**: Verificar si mutaciones están en genes de paneles clínicos
3. **Clasificación funcional**: Distinguir entre oncogenes y supresores tumorales
4. **Análisis de evidencia**: Filtrar por nivel de evidencia (COSMIC Tier)
5. **Investigación traslacional**: Conectar hallazgos con conocimiento clínico

## Consideraciones Técnicas

- **Memoria**: Optimizado para archivos grandes mediante DuckDB
- **Compatibilidad**: Soporta archivos comprimidos y no comprimidos
- **Robustez**: Manejo automático de diferentes formatos de columnas
- **Logging**: Información detallada del proceso de anotación
- **Limpieza**: Elimina automáticamente archivos temporales

## Manejo de Errores

```python
try:
    annotated_df, output_path = knownCancer(
        maf_file="nonexistent.maf",
        annotation_table="cosmic.tsv"
    )
except FileNotFoundError as e:
    print(f"Archivo no encontrado: {e}")

try:
    annotated_df, output_path = knownCancer(
        maf_file="mutations.maf",
        annotation_table="cosmic.tsv",
        join_column="NonexistentColumn"
    )
except ValueError as e:
    print(f"Error de columna: {e}")
```

## Notas Importantes

- Los datos COSMIC requieren registro y descarga desde el sitio oficial
- OncoKB requiere licencia académica o comercial según el uso
- La función preserva todas las columnas originales del MAF
- El campo `Is_Oncogene_any` facilita análisis posteriores
- Compatible con flujos de trabajo de análisis de mutaciones somáticas