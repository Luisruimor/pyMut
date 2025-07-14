# VEP Annotation - Anotación VEP

La función **wrap_maf_vep_annotate** permite anotar archivos MAF utilizando VEP (Variant Effect Predictor) de Ensembl, proporcionando predicciones detalladas del efecto funcional de las variantes.

## ¿Qué es la Anotación VEP?

VEP (Variant Effect Predictor) es una herramienta de Ensembl que predice el efecto funcional de variantes genómicas. Esta función automatiza el proceso completo: convierte archivos MAF al formato requerido por VEP, ejecuta la anotación, y fusiona los resultados de vuelta al formato MAF original.

## Características Principales

- **Proceso automatizado**: Conversión MAF → región → VEP → fusión automática
- **Detección automática**: Extrae automáticamente assembly y versión del cache
- **Anotación completa**: Incluye proteínas, dominios, Uniprot y símbolos génicos
- **Gestión de sinónimos**: Manejo automático de sinónimos cromosómicos
- **Limpieza automática**: Elimina archivos temporales automáticamente
- **Logging detallado**: Información completa del proceso de anotación
- **Compresión opcional**: Soporte para archivos comprimidos de salida

## Uso Básico

```python
from pyMut.annotate.vep_annotate import wrap_maf_vep_annotate

# Anotación VEP básica
success, output_info = wrap_maf_vep_annotate(
    maf_file="mutations.maf",
    cache_dir="/path/to/vep_cache",
    fasta="/path/to/reference.fa"
)

if success:
    print(f"Anotación exitosa: {output_info}")
else:
    print(f"Error en anotación: {output_info}")
```

## Uso Avanzado

```python
# Anotación con parámetros específicos
success, output_info = wrap_maf_vep_annotate(
    maf_file="large_dataset.maf.gz",
    cache_dir="/data/vep_cache",
    fasta="/data/GRCh38.fa",
    output_file="/results/vep_annotations.txt",
    synonyms_file="/data/chr_synonyms.txt",
    assembly="GRCh38",
    version="110",
    compress=True
)

print(f"Estado: {'Éxito' if success else 'Error'}")
print(f"Información: {output_info}")
```

## Parámetros

### maf_file (str | Path, requerido)
- **Descripción**: Ruta al archivo MAF de entrada
- **Formatos soportados**: `.maf`, `.maf.gz`
- **Ejemplo**: `"mutations.maf"` o `"data/variants.maf.gz"`

### cache_dir (str | Path, requerido)
- **Descripción**: Ruta al directorio de cache de VEP
- **Estructura esperada**: `cache_dir/homo_sapiens/version_assembly/`
- **Ejemplo**: `"/opt/vep/cache"` o `"/data/vep_cache"`

### fasta (str | Path, requerido)
- **Descripción**: Ruta al archivo FASTA de referencia
- **Formato**: Archivo FASTA indexado (.fa, .fasta)
- **Ejemplo**: `"/data/GRCh38.fa"` o `"/ref/hg19.fasta"`

### output_file (str | Path, opcional)
- **Descripción**: Ruta del archivo de salida VEP
- **Default**: Se crea directorio automático `vep_annotation_HHMMDDMM`
- **Ejemplo**: `"/results/annotations.txt"`

### synonyms_file (str | Path, opcional)
- **Descripción**: Archivo de sinónimos cromosómicos
- **Default**: Se construye automáticamente desde cache
- **Ruta automática**: `cache_dir/homo_sapiens/version_assembly/chr_synonyms.txt`

### assembly (str, opcional)
- **Descripción**: Nombre del ensamblaje genómico
- **Default**: Se extrae automáticamente del nombre del cache
- **Ejemplos**: `"GRCh38"`, `"GRCh37"`, `"hg19"`

### version (str, opcional)
- **Descripción**: Versión del cache de VEP
- **Default**: Se extrae automáticamente del nombre del cache
- **Ejemplos**: `"110"`, `"109"`, `"108"`

### compress (bool, default=True)
- **Descripción**: Si comprimir el archivo MAF fusionado final
- **True**: Genera archivo `.maf.gz`
- **False**: Genera archivo `.maf` sin comprimir

## Valor de Retorno

Retorna una tupla con dos elementos:

```python
success, output_info = wrap_maf_vep_annotate(...)

# success: bool - True si la anotación fue exitosa
# output_info: str - Información sobre archivos generados o errores
```

### Casos de Retorno

#### Éxito Completo
```python
success = True
output_info = "VEP folder: /path/to/vep_output.txt, Merged file: /path/to/merged.maf.gz"
```

#### Éxito VEP, Error en Fusión
```python
success = True
output_info = "VEP folder: /path/to/vep_output.txt, Merge failed: error_message"
```

#### Error en VEP
```python
success = False
output_info = "/path/to/attempted_output.txt"
```

## Parámetros VEP Utilizados

La función utiliza parámetros VEP fijos optimizados para análisis de mutaciones:

```bash
vep \
  --input_file regions.txt \
  --format region \
  --offline --cache --cache_version 110 \
  --dir_cache /path/to/cache \
  --assembly GRCh38 \
  --synonyms /path/to/chr_synonyms.txt \
  --fasta /path/to/reference.fa \
  --protein --uniprot --domains --symbol \
  --pick \
  --no_stats \
  --output_file output.txt
```

### Explicación de Parámetros VEP

- **--offline --cache**: Modo offline usando cache local
- **--protein**: Incluye consecuencias en proteínas
- **--uniprot**: Añade IDs de Uniprot
- **--domains**: Incluye información de dominios proteicos
- **--symbol**: Añade símbolos génicos
- **--pick**: Selecciona una transcripción por variante
- **--no_stats**: No genera estadísticas (más rápido)

## Anotaciones VEP Añadidas

La función añade columnas con prefijo `VEP_` al MAF original:

### Información Básica
- **VEP_Allele**: Alelo alternativo
- **VEP_Consequence**: Consecuencia funcional principal
- **VEP_IMPACT**: Impacto de la variante (HIGH, MODERATE, LOW, MODIFIER)
- **VEP_SYMBOL**: Símbolo del gen
- **VEP_Gene**: ID del gen Ensembl

### Información de Transcripción
- **VEP_Feature_type**: Tipo de característica (Transcript, RegulatoryFeature)
- **VEP_Feature**: ID de la transcripción/característica
- **VEP_BIOTYPE**: Biotipo de la transcripción
- **VEP_EXON**: Número de exón
- **VEP_INTRON**: Número de intrón

### Cambios en Proteína
- **VEP_HGVSc**: Nomenclatura HGVS para DNA codificante
- **VEP_HGVSp**: Nomenclatura HGVS para proteína
- **VEP_cDNA_position**: Posición en cDNA
- **VEP_CDS_position**: Posición en secuencia codificante
- **VEP_Protein_position**: Posición en proteína
- **VEP_Amino_acids**: Cambio de aminoácidos
- **VEP_Codons**: Cambio de codones

### Información Adicional
- **VEP_STRAND**: Cadena del gen
- **VEP_SWISSPROT**: ID de SwissProt/Uniprot
- **VEP_DOMAINS**: Dominios proteicos afectados
- **VEP_DISTANCE**: Distancia a la característica más cercana

## Ejemplos de Análisis

### Análisis de Consecuencias

```python
# Después de la anotación VEP
success, output_info = wrap_maf_vep_annotate(
    maf_file="mutations.maf",
    cache_dir="/vep_cache",
    fasta="/reference.fa"
)

if success:
    # Cargar el archivo fusionado
    import pandas as pd
    merged_file = output_info.split("Merged file: ")[1]
    df = pd.read_csv(merged_file, sep='\t', compression='gzip' if merged_file.endswith('.gz') else None)
    
    # Análisis de consecuencias
    consequences = df['VEP_Consequence'].value_counts()
    print("Distribución de consecuencias:")
    print(consequences.head(10))
    
    # Variantes de alto impacto
    high_impact = df[df['VEP_IMPACT'] == 'HIGH']
    print(f"Variantes de alto impacto: {len(high_impact)}")
```

### Filtrado por Tipo de Variante

```python
# Filtrar variantes missense
missense_variants = df[df['VEP_Consequence'].str.contains('missense_variant', na=False)]

# Filtrar variantes nonsense
nonsense_variants = df[df['VEP_Consequence'].str.contains('stop_gained', na=False)]

# Filtrar variantes en sitios de splicing
splice_variants = df[df['VEP_Consequence'].str.contains('splice', na=False)]

print(f"Variantes missense: {len(missense_variants)}")
print(f"Variantes nonsense: {len(nonsense_variants)}")
print(f"Variantes de splicing: {len(splice_variants)}")
```

### Análisis de Dominios Proteicos

```python
# Variantes que afectan dominios conocidos
domain_variants = df[df['VEP_DOMAINS'].notna()]

# Dominios más frecuentemente afectados
domain_counts = df['VEP_DOMAINS'].str.split('&').explode().value_counts()
print("Dominios más afectados:")
print(domain_counts.head())

# Genes con más variantes en dominios
genes_with_domains = domain_variants.groupby('VEP_SYMBOL').size().sort_values(ascending=False)
print("Genes con más variantes en dominios:")
print(genes_with_domains.head())
```

## Configuración del Entorno

### Instalación de VEP

```bash
# Instalar VEP usando conda
conda install -c bioconda ensembl-vep

# O usando Docker
docker pull ensemblorg/ensembl-vep
```

### Descarga de Cache

```bash
# Descargar cache para GRCh38
vep_install -a cf -s homo_sapiens -y GRCh38 -c /path/to/cache --CONVERT

# Descargar cache para GRCh37
vep_install -a cf -s homo_sapiens -y GRCh37 -c /path/to/cache --CONVERT
```

### Estructura de Directorios

```
/path/to/vep_cache/
├── homo_sapiens/
│   ├── 110_GRCh38/
│   │   ├── chr_synonyms.txt
│   │   └── [otros archivos de cache]
│   └── 110_GRCh37/
│       ├── chr_synonyms.txt
│       └── [otros archivos de cache]
```

## Manejo de Errores

### Errores Comunes

```python
try:
    success, output_info = wrap_maf_vep_annotate(
        maf_file="mutations.maf",
        cache_dir="/vep_cache",
        fasta="/reference.fa"
    )
except FileNotFoundError as e:
    print(f"Archivo no encontrado: {e}")
except ValueError as e:
    print(f"Error de configuración: {e}")
except Exception as e:
    print(f"Error inesperado: {e}")
```

### Diagnóstico de Problemas

```python
import logging
logging.basicConfig(level=logging.INFO)

# El logging mostrará información detallada:
"""
INFO - Converting MAF file to region format: mutations.maf
INFO - Successfully converted MAF to region format: regions.txt
INFO - Extracted from cache: assembly=GRCh38, version=110
INFO - Auto-constructed chr synonyms path: /cache/homo_sapiens/110_GRCh38/chr_synonyms.txt
INFO - Running VEP annotation: vep --input_file regions.txt ...
INFO - VEP annotation completed successfully
INFO - Merging VEP annotations with original MAF file...
INFO - Successfully merged VEP annotations. Merged file: merged.maf.gz
"""
```

## Consideraciones de Rendimiento

### Para Archivos Grandes
- **Memoria**: VEP puede requerir varios GB de RAM
- **Tiempo**: El proceso puede tomar horas para datasets grandes
- **Almacenamiento**: Los archivos anotados son significativamente más grandes

### Optimizaciones
```python
# Para datasets grandes, considerar procesamiento por lotes
import pandas as pd

# Dividir MAF grande en chunks
chunk_size = 10000
maf_df = pd.read_csv("large_mutations.maf", sep='\t', chunksize=chunk_size)

for i, chunk in enumerate(maf_df):
    chunk_file = f"chunk_{i}.maf"
    chunk.to_csv(chunk_file, sep='\t', index=False)
    
    success, output_info = wrap_maf_vep_annotate(
        maf_file=chunk_file,
        cache_dir="/vep_cache",
        fasta="/reference.fa"
    )
```

## Casos de Uso Comunes

1. **Predicción de efectos**: Determinar el impacto funcional de variantes
2. **Priorización clínica**: Identificar variantes de alto impacto
3. **Análisis de dominios**: Estudiar variantes en dominios proteicos específicos
4. **Anotación funcional**: Añadir información detallada para análisis downstream
5. **Control de calidad**: Validar la clasificación de variantes

## Archivos Generados

La función genera varios archivos:

```
vep_annotation_HHMMDDMM/
├── mutations_vep.txt              # Anotaciones VEP raw
└── mutations_vep_merged.maf.gz    # MAF fusionado con anotaciones VEP
```

## Integración con Otros Análisis

```python
# Combinar VEP con otras anotaciones
from pyMut.annotate.cosmic_cancer_annotate import knownCancer

# Primero VEP
success, vep_output = wrap_maf_vep_annotate(
    maf_file="mutations.maf",
    cache_dir="/vep_cache",
    fasta="/reference.fa"
)

if success:
    # Extraer archivo fusionado
    merged_file = vep_output.split("Merged file: ")[1]
    
    # Luego COSMIC
    cosmic_df, cosmic_output = knownCancer(
        maf_file=merged_file,
        annotation_table="cosmic.tsv"
    )
    
    print(f"Archivo final con VEP + COSMIC: {cosmic_output}")
```

## Notas Técnicas

- **Dependencias**: Requiere VEP instalado y accesible en PATH
- **Cache**: Necesita cache de VEP descargado localmente
- **Memoria**: Uso de memoria proporcional al tamaño del archivo
- **Limpieza**: Archivos temporales se eliminan automáticamente
- **Compatibilidad**: Funciona con VEP v104+ y cache correspondiente

## Limitaciones

- **Dependencia externa**: Requiere instalación de VEP
- **Tamaño de cache**: Los caches de VEP son grandes (varios GB)
- **Tiempo de procesamiento**: Puede ser lento para datasets grandes
- **Formato fijo**: Utiliza parámetros VEP predefinidos
- **Una transcripción**: `--pick` selecciona solo una transcripción por variante