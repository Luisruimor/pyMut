# read_maf - Lectura de Archivos MAF

La función **read_maf** es la función principal para cargar archivos MAF (Mutation Annotation Format) en pyMut y convertirlos en objetos PyMutation.

## ¿Qué es read_maf?

Es una función que lee archivos MAF (formato estándar para datos de mutaciones) y los convierte en objetos PyMutation que pueden ser utilizados para análisis y visualización de mutaciones.

## Características Principales

- **Soporte para archivos comprimidos**: Lee archivos `.maf` y `.maf.gz`
- **Validación automática**: Verifica que el archivo tenga las columnas requeridas
- **Conversión a formato VCF-like**: Convierte los datos MAF al formato interno de pyMut
- **Manejo de comentarios**: Preserva los comentarios del archivo MAF en los metadatos
- **Expansión de muestras**: Convierte el formato largo de MAF a formato ancho con columnas por muestra

## Uso Básico

```python
from pyMut.io import read_maf

# Cargar archivo MAF simple
py_mut = read_maf("mutations.maf")

# Cargar archivo MAF comprimido
py_mut = read_maf("mutations.maf.gz")

# Cargar con archivo FASTA de referencia
py_mut = read_maf("mutations.maf", fasta="reference.fasta")
```

## Parámetros

### path (str | Path) [requerido]
- **Descripción**: Ruta al archivo MAF
- **Formatos soportados**: `.maf`, `.maf.gz`
- **Ejemplo**: `"data/mutations.maf"`

### fasta (str | Path, opcional)
- **Descripción**: Ruta al archivo FASTA de referencia
- **Uso**: Se incluye en los metadatos del objeto PyMutation resultante
- **Ejemplo**: `"reference/hg38.fasta"`

## Valor de Retorno

Retorna un objeto **PyMutation** que contiene:
- **data**: DataFrame con las mutaciones en formato VCF-like
- **samples**: Lista de muestras detectadas en el archivo
- **metadata**: Información sobre el archivo fuente, comentarios y configuración

## Formato de Datos MAF

### Columnas Requeridas
```
Hugo_Symbol              - Símbolo del gen
Variant_Classification   - Tipo de mutación
Chromosome              - Cromosoma
Start_Position          - Posición de inicio
Reference_Allele        - Alelo de referencia
Tumor_Seq_Allele2       - Alelo alternativo
Tumor_Sample_Barcode    - Identificador de muestra
```

### Columnas Opcionales
```
dbSNP_RS                - Identificador dbSNP
Tumor_Seq_Allele1       - Alelo alternativo secundario
End_Position            - Posición final
Variant_Type            - Tipo de variante (SNP, INS, DEL)
```

## Ejemplo Completo

```python
from pyMut.io import read_maf
from pyMut import PyMutation

# Cargar datos MAF
py_mut = read_maf("src/pyMut/data/examples/tcga_laml.maf.gz")

# Verificar carga exitosa
print(f"Muestras cargadas: {len(py_mut.samples)}")
print(f"Mutaciones totales: {len(py_mut.data)}")
print(f"Columnas disponibles: {list(py_mut.data.columns)}")

# Configurar alta calidad para visualizaciones
PyMutation.configure_high_quality_plots()

# Generar análisis
summary_fig = py_mut.summary_plot(title="Análisis TCGA-LAML")
summary_fig.savefig("tcga_laml_summary.png")

# Análisis TMB
tmb_results = py_mut.calculate_tmb_analysis(
    genome_size_bp=60456963,  # WES
    output_dir="results"
)

print("✅ Análisis completado exitosamente!")
```

## Conversión de Formato

La función convierte automáticamente el formato MAF al formato interno de pyMut:

### Formato MAF (entrada)
```
Hugo_Symbol | Variant_Classification | Tumor_Sample_Barcode | Reference_Allele | Tumor_Seq_Allele2
GENE1      | Missense_Mutation      | SAMPLE_001          | A               | G
GENE2      | Nonsense_Mutation      | SAMPLE_002          | C               | T
```

### Formato pyMut (salida)
```
CHROM | POS | REF | ALT | SAMPLE_001 | SAMPLE_002 | Hugo_Symbol | Variant_Classification
chr1  | 100 | A   | G   | A|G        | A|A        | GENE1       | Missense_Mutation
chr2  | 200 | C   | T   | C|C        | C|T        | GENE2       | Nonsense_Mutation
```

## Manejo de Errores

### FileNotFoundError
```python
try:
    py_mut = read_maf("archivo_inexistente.maf")
except FileNotFoundError:
    print("❌ Archivo no encontrado")
```

### ValueError (columnas faltantes)
```python
try:
    py_mut = read_maf("archivo_incompleto.maf")
except ValueError as e:
    print(f"❌ Error en formato MAF: {e}")
```

### ImportError (pyarrow no disponible)
```python
# La función automáticamente usa el motor 'c' como alternativa
# No requiere manejo especial del usuario
```

## Metadatos Incluidos

El objeto PyMutation resultante incluye metadatos completos:

```python
py_mut = read_maf("mutations.maf")

# Acceder a metadatos
print(f"Formato fuente: {py_mut.metadata.source_format}")
print(f"Archivo: {py_mut.metadata.file_path}")
print(f"Cargado en: {py_mut.metadata.loaded_at}")
print(f"Comentarios: {py_mut.metadata.notes}")
```

## Optimización y Rendimiento

- **Motor pyarrow**: Utiliza pyarrow por defecto para mejor rendimiento
- **Fallback automático**: Si pyarrow falla, usa el motor 'c' automáticamente
- **Manejo de memoria**: Optimizado para archivos grandes con `low_memory=False`
- **Compresión**: Soporte nativo para archivos `.gz` sin descompresión manual

## Solución de Problemas

### Archivo muy grande
```python
# Para archivos muy grandes, monitorear el progreso
import logging
logging.basicConfig(level=logging.INFO)

py_mut = read_maf("archivo_grande.maf.gz")
```

### Columnas con nombres no estándar
```python
# La función normaliza automáticamente nombres de columnas comunes
# No requiere acción del usuario
```

### Archivos con muchos comentarios
```python
# Los comentarios se preservan automáticamente en metadata.notes
py_mut = read_maf("archivo_con_comentarios.maf")
print(py_mut.metadata.notes)  # Ver comentarios preservados
```