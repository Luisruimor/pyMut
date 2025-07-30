# COSMIC Annotation - Anotación COSMIC

La función **_maf_COSMIC_OncoKB_annotation_aux** permite anotar archivos MAF con información completa del COSMIC Cancer Gene Census y opcionalmente OncoKB, proporcionando todas las anotaciones disponibles.

## ¿Qué es la Anotación COSMIC?

COSMIC (Catalogue of Somatic Mutations in Cancer) es la base de datos más completa del mundo sobre mutaciones somáticas en cáncer. Esta función añade información detallada sobre genes relacionados con cáncer, incluyendo su rol oncológico, nivel de evidencia, y clasificaciones funcionales.

## Características Principales

- **Anotación completa**: Incluye todas las columnas disponibles de COSMIC y OncoKB
- **Estrategia híbrida**: Optimización automática según el tamaño del archivo
- **Manejo de sinónimos**: Resolución automática de nombres alternativos de genes
- **Múltiples formatos**: Soporte para archivos comprimidos y no comprimidos
- **Logging detallado**: Información completa del proceso de anotación
- **Limpieza automática**: Gestión de archivos temporales

## Uso Básico

```python
from pyMut.annotate.cosmic_cancer_annotate import _maf_COSMIC_OncoKB_annotation_aux

# Anotación COSMIC básica
annotated_df, output_path = _maf_COSMIC_OncoKB_annotation_aux(
    maf_file="mutations.maf",
    annotation_table="cosmic_cancer_gene_census.tsv"
)

print(f"Archivo anotado: {output_path}")
print(f"Columnas añadidas: {[col for col in annotated_df.columns if col.startswith('COSMIC_')]}")
```

## Uso Avanzado con OncoKB

```python
# Anotación completa con COSMIC y OncoKB
annotated_df, output_path = _maf_COSMIC_OncoKB_annotation_aux(
    maf_file="large_dataset.maf.gz",
    annotation_table="cosmic_v102.tsv.gz",
    oncokb_table="oncokb_cancer_gene_list.tsv",
    output_path="fully_annotated.maf.gz",
    compress_output=True,
    join_column="Hugo_Symbol"
)

# Explorar anotaciones añadidas
cosmic_cols = [col for col in annotated_df.columns if col.startswith('COSMIC_')]
oncokb_cols = [col for col in annotated_df.columns if col.startswith('OncoKB_')]

print(f"Columnas COSMIC añadidas: {len(cosmic_cols)}")
print(f"Columnas OncoKB añadidas: {len(oncokb_cols)}")
```

## Parámetros

### maf_file (str | Path, requerido)
- **Descripción**: Ruta al archivo MAF de entrada
- **Formatos soportados**: `.maf`, `.maf.gz`
- **Detección automática**: Maneja compresión automáticamente

### annotation_table (str | Path, requerido)
- **Descripción**: Ruta a la tabla COSMIC Cancer Gene Census
- **Formatos soportados**: `.tsv`, `.tsv.gz`
- **Versión recomendada**: COSMIC v102 o superior

### output_path (str | Path, opcional)
- **Descripción**: Ruta del archivo de salida
- **Default**: Se genera automáticamente con sufijo apropiado
- **Ejemplos**: 
  - Solo COSMIC: `_COSMIC_annotated.maf`
  - Con OncoKB: `_COSMIC_OncoKB_annotated.maf`

### compress_output (bool, default=True)
- **Descripción**: Si comprimir el archivo de salida
- **Recomendación**: `True` para archivos grandes

### join_column (str, default="Hugo_Symbol")
- **Descripción**: Columna para el join entre archivos
- **Aliases soportados**: Definidos en `fields.py`
- **Detección automática**: Busca columnas equivalentes

### oncokb_table (str | Path, opcional)
- **Descripción**: Ruta a la tabla OncoKB cancer gene list
- **Efecto**: Añade anotaciones OncoKB completas
- **Formato**: `.tsv` con columnas estándar OncoKB

## Valor de Retorno

```python
annotated_df, output_path = _maf_COSMIC_OncoKB_annotation_aux(...)

# annotated_df: pandas.DataFrame con todas las anotaciones
# output_path: pathlib.Path del archivo generado
```

## Anotaciones COSMIC Incluidas

La función añade todas las columnas disponibles del COSMIC Cancer Gene Census:

### Información Básica del Gen
- **COSMIC_Gene Symbol**: Símbolo oficial del gen
- **COSMIC_Name**: Nombre completo del gen
- **COSMIC_Entrez GeneId**: ID de Entrez Gene
- **COSMIC_Genome Location**: Localización genómica

### Clasificación Funcional
- **COSMIC_Tier**: Nivel de evidencia (1, 2, 3)
- **COSMIC_Hallmark**: Características del cáncer asociadas
- **COSMIC_Chr Band**: Banda cromosómica
- **COSMIC_Somatic**: Si tiene mutaciones somáticas
- **COSMIC_Germline**: Si tiene mutaciones germinales

### Rol en Cáncer
- **COSMIC_Tumour Types(Somatic)**: Tipos de tumor con mutaciones somáticas
- **COSMIC_Tumour Types(Germline)**: Tipos de tumor con mutaciones germinales
- **COSMIC_Cancer Syndrome**: Síndromes de cáncer asociados
- **COSMIC_Tissue Type**: Tipos de tejido afectados
- **COSMIC_Molecular Genetics**: Genética molecular
- **COSMIC_Role in Cancer**: Rol específico en cáncer

### Información Adicional
- **COSMIC_Mutation Types**: Tipos de mutación observados
- **COSMIC_Translocation Partner**: Parejas de translocación
- **COSMIC_Other Germline Mut**: Otras mutaciones germinales
- **COSMIC_Other Syndrome**: Otros síndromes asociados

## Anotaciones OncoKB Incluidas

Si se proporciona `oncokb_table`, se añaden todas las columnas OncoKB:

### Clasificación Oncológica
- **OncoKB_Is Oncogene**: Si es oncogén
- **OncoKB_Is Tumor Suppressor Gene**: Si es supresor tumoral
- **OncoKB_OncoKB Annotated**: Si está anotado en OncoKB

### Paneles Clínicos
- **OncoKB_MSK-IMPACT**: Panel MSK-IMPACT
- **OncoKB_MSK-HEME**: Panel hematológico MSK
- **OncoKB_FOUNDATION ONE**: Panel Foundation One
- **OncoKB_FOUNDATION ONE HEME**: Panel Foundation One hematológico
- **OncoKB_Vogelstein**: Lista de Vogelstein
- **OncoKB_SANGER CGC(05/30/2017)**: Cancer Gene Census histórico

### Información Adicional
- **OncoKB_Gene Aliases**: Sinónimos del gen
- **OncoKB_Background**: Información de fondo

## Optimización de Rendimiento

La función utiliza una estrategia híbrida inteligente:

```python
# Detección automática del método óptimo
if file_size >= 2_000_000_000:  # 2GB
    # Usa DuckDB para archivos grandes
    method = "DuckDB"
else:
    # Usa pandas + pyarrow para archivos pequeños
    method = "pandas"

print(f"Método seleccionado: {method}")
```

### DuckDB (Archivos Grandes ≥2GB)
- **Ventajas**: Eficiencia de memoria, procesamiento paralelo
- **Uso**: Datasets grandes, servidores con memoria limitada
- **Rendimiento**: Óptimo para archivos multi-GB

### Pandas + PyArrow (Archivos Pequeños <2GB)
- **Ventajas**: Velocidad máxima, compatibilidad completa
- **Uso**: Análisis interactivo, datasets medianos
- **Rendimiento**: Óptimo para archivos hasta 2GB

## Manejo de Sinónimos

La función resuelve automáticamente sinónimos de genes:

```python
# Sinónimos COSMIC (columna "SYNONIMS")
cosmic_synonyms = {
    "ERBB2": ["HER2", "NEU"],
    "CDKN2A": ["P16", "INK4A"],
    "TP53": ["P53"]
}

# Sinónimos OncoKB (columna "Gene Aliases")
oncokb_synonyms = {
    "ERBB2": "HER2, NEU, MLN19",
    "BRAF": "B-RAF1, BRAF1"
}

# Resolución automática durante el join
```

## Ejemplos de Análisis

### Exploración de Anotaciones

```python
annotated_df, _ = _maf_COSMIC_OncoKB_annotation_aux(
    maf_file="mutations.maf",
    annotation_table="cosmic.tsv",
    oncokb_table="oncokb.tsv"
)

# Genes por tier COSMIC
tier_counts = annotated_df['COSMIC_Tier'].value_counts()
print("Distribución por Tier COSMIC:")
print(tier_counts)

# Roles en cáncer
roles = annotated_df['COSMIC_Role in Cancer'].value_counts()
print("\nRoles en cáncer:")
print(roles.head())

# Tipos de tumor más frecuentes
tumor_types = annotated_df['COSMIC_Tumour Types(Somatic)'].value_counts()
print("\nTipos de tumor (top 5):")
print(tumor_types.head())
```

### Análisis Comparativo de Fuentes

```python
# Comparar clasificaciones entre COSMIC y OncoKB
comparison = annotated_df.groupby([
    'COSMIC_Role in Cancer', 
    'OncoKB_Is Oncogene'
]).size().reset_index(name='count')

print("Comparación COSMIC vs OncoKB:")
print(comparison)

# Genes solo en COSMIC
cosmic_only = annotated_df[
    annotated_df['COSMIC_Role in Cancer'].notna() & 
    annotated_df['OncoKB_Is Oncogene'].isna()
]['Hugo_Symbol'].unique()

print(f"Genes solo en COSMIC: {len(cosmic_only)}")
```

### Filtrado por Evidencia

```python
# Genes con evidencia alta (Tier 1)
tier1_genes = annotated_df[
    annotated_df['COSMIC_Tier'] == '1'
]

# Genes con mutaciones somáticas y germinales
dual_mutations = annotated_df[
    (annotated_df['COSMIC_Somatic'] == 'yes') & 
    (annotated_df['COSMIC_Germline'] == 'yes')
]

# Genes en paneles clínicos múltiples
multi_panel_genes = annotated_df[
    (annotated_df['OncoKB_MSK-IMPACT'] == 'True') & 
    (annotated_df['OncoKB_FOUNDATION ONE'] == 'True')
]

print(f"Genes Tier 1: {len(tier1_genes)}")
print(f"Genes con mutaciones duales: {len(dual_mutations)}")
print(f"Genes en múltiples paneles: {len(multi_panel_genes)}")
```

## Archivos de Salida

La función genera archivos con nombres descriptivos automáticamente:

```python
# Ejemplos de nombres generados
"mutations_COSMIC_annotated.maf.gz"                    # Solo COSMIC
"large_dataset_COSMIC_OncoKB_annotated.maf.gz"        # COSMIC + OncoKB
"sample_data_COSMIC_annotated.maf"                     # Sin compresión
```

## Logging y Monitoreo

La función proporciona logging detallado:

```python
import logging
logging.basicConfig(level=logging.INFO)

# Ejemplo de salida de logging
"""
INFO - Starting COSMIC annotation process
INFO - MAF file: mutations.maf (1.2GB)
INFO - Annotation table: cosmic_v102.tsv.gz
INFO - OncoKB table: oncokb_cancer_gene_list.tsv
INFO - Using DuckDB for large file processing
INFO - Processing 45,231 mutations
INFO - Found 1,247 genes with COSMIC annotations
INFO - Found 892 genes with OncoKB annotations
INFO - Annotation completed successfully
INFO - Output file saved: mutations_COSMIC_OncoKB_annotated.maf.gz
"""
```

## Consideraciones de Memoria

### Para Archivos Grandes
```python
# Configuración recomendada para archivos >5GB
annotated_df, output_path = _maf_COSMIC_OncoKB_annotation_aux(
    maf_file="huge_dataset.maf.gz",
    annotation_table="cosmic.tsv.gz",
    compress_output=True,  # Importante para ahorrar espacio
    # DuckDB se selecciona automáticamente
)
```

### Para Análisis Interactivo
```python
# Configuración para análisis rápido
annotated_df, output_path = _maf_COSMIC_OncoKB_annotation_aux(
    maf_file="small_sample.maf",
    annotation_table="cosmic.tsv",
    compress_output=False,  # Acceso más rápido
    # pandas se selecciona automáticamente
)
```

## Casos de Uso Comunes

1. **Investigación básica**: Caracterización completa de genes mutados
2. **Análisis clínico**: Identificación de variantes accionables
3. **Estudios comparativos**: Análisis entre diferentes cohortes
4. **Desarrollo de paneles**: Selección de genes para paneles clínicos
5. **Validación de hallazgos**: Confirmación con bases de datos curadas

## Manejo de Errores

```python
import logging

try:
    annotated_df, output_path = _maf_COSMIC_OncoKB_annotation_aux(
        maf_file="mutations.maf",
        annotation_table="cosmic.tsv"
    )
except FileNotFoundError as e:
    logging.error(f"Archivo no encontrado: {e}")
except ValueError as e:
    logging.error(f"Error de validación: {e}")
except Exception as e:
    logging.error(f"Error inesperado: {e}")
```

## Notas Técnicas

- **Compatibilidad**: Funciona con todas las versiones de COSMIC v95+
- **Memoria**: Gestión automática según el tamaño del archivo
- **Paralelización**: DuckDB utiliza múltiples cores automáticamente
- **Limpieza**: Archivos temporales se eliminan automáticamente
- **Robustez**: Manejo de errores y validación de entrada

## Diferencias con knownCancer

| Aspecto | _maf_COSMIC_OncoKB_annotation_aux | knownCancer |
|---------|-----------------------------------|-------------|
| **Propósito** | Anotación completa | Anotación específica de cáncer |
| **Columnas** | Todas las disponibles | Solo relevantes para cáncer |
| **Tamaño resultado** | Mayor | Menor |
| **Uso recomendado** | Investigación exploratoria | Análisis clínico |
| **Rendimiento** | Más lento | Más rápido |

## Recomendaciones de Uso

- **Usar para investigación**: Cuando necesites todas las anotaciones disponibles
- **Combinar con filtrado**: Aplicar filtros posteriores según necesidades
- **Monitorear memoria**: Especialmente importante para archivos grandes
- **Validar resultados**: Verificar que las anotaciones son las esperadas
- **Documentar versiones**: Registrar versiones de COSMIC y OncoKB utilizadas