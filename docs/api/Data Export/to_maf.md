# to_maf - Exportación a Formato MAF

El método **to_maf** permite exportar datos de un objeto PyMutation de vuelta al formato MAF (Mutation Annotation Format) estándar.

## ¿Qué es to_maf?

Es un método que convierte los datos internos de pyMut (formato VCF-like) de vuelta al formato MAF estándar, permitiendo la interoperabilidad con otras herramientas de análisis de mutaciones.

## Características Principales

- **Conversión automática de formato**: Convierte de formato interno pyMut a MAF estándar
- **Preservación de columnas**: Mantiene todas las columnas MAF originales cuando están disponibles
- **Orden de columnas estándar**: Utiliza el orden de columnas MAF oficial
- **Manejo de cromosomas**: Convierte automáticamente el formato de cromosomas
- **Validación de datos**: Verifica que los datos sean compatibles con el formato MAF
- **Soporte para compresión**: Puede generar archivos `.maf` o `.maf.gz`

## Uso Básico

```python
from pyMut.input import read_maf

# Cargar datos MAF
py_mut = read_maf("original.maf")

# Aplicar algunos filtros o análisis
filtered = py_mut.filter_by_chrom_sample(chrom="chr17")

# Exportar de vuelta a MAF
filtered.to_maf("filtered_chr17.maf")
```

## Parámetros

### output_path (str | Path) [requerido]
- **Descripción**: Ruta donde guardar el archivo MAF de salida
- **Formatos soportados**: `.maf`, `.maf.gz`
- **Ejemplo**: `"results/filtered_mutations.maf"`

## Conversión de Formato

### Formato pyMut (entrada)
```
CHROM | POS | REF | ALT | SAMPLE_001 | SAMPLE_002 | Hugo_Symbol | Variant_Classification
chr17 | 100 | A   | G   | A|G        | A|A        | TP53        | Missense_Mutation
chr17 | 200 | C   | T   | C|C        | C|T        | BRCA1       | Nonsense_Mutation
```

### Formato MAF (salida)
```
Hugo_Symbol | Variant_Classification | Tumor_Sample_Barcode | Chromosome | Start_Position | Reference_Allele | Tumor_Seq_Allele2
TP53        | Missense_Mutation      | SAMPLE_001          | 17         | 100           | A               | G
BRCA1       | Nonsense_Mutation      | SAMPLE_002          | 200        | 200           | C               | T
```

## Ejemplo Completo

```python
from pyMut.input import read_maf
import os

# Cargar datos originales
py_mut = read_maf("src/pyMut/data/examples/tcga_laml.maf.gz")
print(f"Mutaciones originales: {len(py_mut.data)}")

# Aplicar filtros
# Filtrar por cromosoma específico
chr17_mutations = py_mut.filter_by_chrom_sample(chrom="chr17")
print(f"Mutaciones en chr17: {len(chr17_mutations.data)}")

# Filtrar por muestra específica
sample_mutations = py_mut.filter_by_chrom_sample(sample="TCGA-AB-2802")
print(f"Mutaciones en muestra específica: {len(sample_mutations.data)}")

# Crear directorio de resultados
os.makedirs("results", exist_ok=True)

# Exportar datos filtrados
chr17_mutations.to_maf("results/chr17_mutations.maf")
sample_mutations.to_maf("results/sample_mutations.maf.gz")  # Comprimido

print("✅ Archivos MAF exportados exitosamente!")

# Verificar archivos generados
print(f"Archivo chr17: {os.path.exists('results/chr17_mutations.maf')}")
print(f"Archivo muestra: {os.path.exists('results/sample_mutations.maf.gz')}")
```

## Orden de Columnas MAF

El método utiliza el orden estándar de columnas MAF:

```python
# Columnas principales (siempre presentes)
required_columns = [
    'Hugo_Symbol',
    'Entrez_Gene_Id', 
    'Center',
    'NCBI_Build',
    'Chromosome',
    'Start_Position',
    'End_Position',
    'Strand',
    'Variant_Classification',
    'Variant_Type',
    'Reference_Allele',
    'Tumor_Seq_Allele1',
    'Tumor_Seq_Allele2',
    'dbSNP_RS',
    'dbSNP_Val_Status',
    'Tumor_Sample_Barcode',
    'Matched_Norm_Sample_Barcode',
    'Match_Norm_Seq_Allele1',
    'Match_Norm_Seq_Allele2',
    'Tumor_Validation_Allele1',
    'Tumor_Validation_Allele2',
    'Match_Norm_Validation_Allele1',
    'Match_Norm_Validation_Allele2',
    'Verification_Status',
    'Validation_Status',
    'Mutation_Status',
    'Sequencing_Phase',
    'Sequence_Source',
    'Validation_Method',
    'Score',
    'BAM_File',
    'Sequencer',
    'Tumor_Sample_UUID',
    'Matched_Norm_Sample_UUID'
]
```

## Manejo de Cromosomas

La función convierte automáticamente el formato de cromosomas:

```python
# Formato pyMut (entrada) -> Formato MAF (salida)
"chr1"  -> "1"
"chr2"  -> "2"
"chrX"  -> "X"
"chrY"  -> "Y"
"chrM"  -> "MT"
"chr23" -> "X"
"chr24" -> "Y"
```

## Expansión de Muestras

Convierte el formato ancho de pyMut al formato largo de MAF:

### pyMut (formato ancho)
```
CHROM | POS | REF | ALT | SAMPLE_001 | SAMPLE_002 | Hugo_Symbol
chr1  | 100 | A   | G   | A|G        | A|A        | GENE1
```

### MAF (formato largo)
```
Hugo_Symbol | Tumor_Sample_Barcode | Chromosome | Start_Position | Reference_Allele | Tumor_Seq_Allele2
GENE1       | SAMPLE_001          | 1          | 100           | A               | G
```

## Validación de Datos

El método realiza validaciones automáticas:

```python
# Verificaciones realizadas:
# 1. Presencia de columnas requeridas
# 2. Formato válido de cromosomas
# 3. Posiciones numéricas válidas
# 4. Alelos válidos (A, T, G, C, -, N)
# 5. Muestras con datos válidos
```

## Manejo de Errores

### Datos insuficientes
```python
try:
    py_mut.to_maf("output.maf")
except ValueError as e:
    print(f"❌ Error en datos: {e}")
    # Posibles causas:
    # - Faltan columnas requeridas
    # - Formato de cromosoma inválido
    # - Posiciones no numéricas
```

### Problemas de escritura
```python
try:
    py_mut.to_maf("/ruta/protegida/output.maf")
except PermissionError:
    print("❌ Sin permisos para escribir en el directorio")
except FileNotFoundError:
    print("❌ Directorio de destino no existe")
```

## Preservación de Metadatos

```python
# Los metadatos se preservan como comentarios en el archivo MAF
py_mut = read_maf("original.maf")
py_mut.to_maf("copy.maf")

# El archivo copy.maf incluirá:
# #version 2.4
# #source_format: MAF
# #original_file: original.maf
# #exported_at: 2024-01-15 10:30:00
# #filters_applied: chromosome:chr17
```

## Compatibilidad

### Versiones MAF soportadas
- **MAF 2.4**: Formato completo (recomendado)
- **MAF 2.3**: Compatible con limitaciones menores
- **Formatos personalizados**: Soporte básico

### Herramientas compatibles
- **maftools**: Totalmente compatible
- **TCGA**: Formato estándar TCGA
- **cBioPortal**: Compatible para importación
- **MutSigCV**: Compatible para análisis

## Optimización

### Para archivos grandes
```python
# Usar compresión para archivos grandes
py_mut.to_maf("large_dataset.maf.gz")

# Filtrar antes de exportar para reducir tamaño
filtered = py_mut.filter_by_chrom_sample(chrom=["chr1", "chr2"])
filtered.to_maf("subset.maf")
```

### Verificación de integridad
```python
# Exportar y verificar
py_mut.to_maf("test_export.maf")

# Reimportar para verificar
reimported = read_maf("test_export.maf")
print(f"Original: {len(py_mut.data)} mutaciones")
print(f"Reimportado: {len(reimported.data)} mutaciones")
```

## Casos de Uso Comunes

### Workflow de análisis
```python
# 1. Cargar datos
py_mut = read_maf("raw_data.maf")

# 2. Aplicar filtros de calidad
filtered = py_mut.filter_by_chrom_sample(chrom=["chr1", "chr17", "chrX"])

# 3. Realizar análisis TMB
tmb_results = filtered.calculate_tmb_analysis()

# 4. Exportar datos procesados
filtered.to_maf("processed_data.maf")

# 5. Exportar para herramientas externas
filtered.to_maf("for_maftools.maf.gz")
```

### Intercambio de datos
```python
# Preparar datos para colaboradores
py_mut = read_maf("internal_format.maf")

# Filtrar datos sensibles si es necesario
public_data = py_mut.filter_by_chrom_sample(
    sample=approved_samples
)

# Exportar en formato estándar
public_data.to_maf("public_dataset.maf.gz")
```