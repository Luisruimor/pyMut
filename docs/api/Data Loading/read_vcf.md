# read_vcf - Lectura de Archivos VCF

La función **read_vcf** permite cargar archivos VCF (Variant Call Format) en pyMut y convertirlos en objetos PyMutation para análisis de mutaciones.

## ¿Qué es read_vcf?

Es una función que lee archivos VCF (formato estándar para datos de variantes genómicas) y los convierte en objetos PyMutation que pueden ser utilizados para análisis y visualización de mutaciones.

## Características Principales

- **Soporte para archivos comprimidos**: Lee archivos `.vcf` y `.vcf.gz`
- **Indexación automática**: Puede crear índices tabix automáticamente
- **Optimización con pyarrow**: Utiliza pyarrow para mejor rendimiento cuando está disponible
- **Manejo de genotipos**: Convierte genotipos VCF a formato alélico
- **Parsing de INFO**: Extrae información del campo INFO de manera vectorizada
- **Cache inteligente**: Sistema de cache para acelerar cargas repetidas
- **Monitoreo de recursos**: Logging de información del sistema durante la carga

## Uso Básico

```python
from pyMut.input import read_vcf

# Cargar archivo VCF simple
py_mut = read_vcf("variants.vcf")

# Cargar archivo VCF comprimido
py_mut = read_vcf("variants.vcf.gz")

# Cargar con archivo FASTA de referencia
py_mut = read_vcf("variants.vcf", fasta="reference.fasta")

# Crear índice tabix automáticamente
py_mut = read_vcf("variants.vcf.gz", create_index=True)
```

## Parámetros

### path (str | Path) [requerido]
- **Descripción**: Ruta al archivo VCF
- **Formatos soportados**: `.vcf`, `.vcf.gz`
- **Ejemplo**: `"data/variants.vcf.gz"`

### fasta (str | Path, opcional)
- **Descripción**: Ruta al archivo FASTA de referencia
- **Uso**: Se incluye en los metadatos del objeto PyMutation resultante
- **Ejemplo**: `"reference/hg38.fasta"`

### create_index (bool, default=False)
- **Descripción**: Si crear un índice tabix para el archivo VCF
- **Uso**: Mejora el rendimiento para consultas posteriores
- **Nota**: Solo funciona con archivos `.vcf.gz`

### cache_dir (str | Path, opcional)
- **Descripción**: Directorio para almacenar archivos de cache
- **Default**: Directorio temporal del sistema
- **Uso**: Acelera cargas repetidas del mismo archivo

## Valor de Retorno

Retorna un objeto **PyMutation** que contiene:
- **data**: DataFrame con las variantes en formato VCF-like expandido
- **samples**: Lista de muestras detectadas en el archivo VCF
- **metadata**: Información sobre el archivo fuente y configuración

## Formato de Datos VCF

### Columnas Estándar VCF
```
CHROM | POS | ID | REF | ALT | QUAL | FILTER | INFO | FORMAT | SAMPLE_001 | SAMPLE_002
chr1  | 100 | .  | A   | G   | 60   | PASS   | ...  | GT:DP  | 0/1:30    | 1/1:25
```

### Conversión a Formato pyMut
```
CHROM | POS | ID | REF | ALT | QUAL | FILTER | SAMPLE_001 | SAMPLE_002 | INFO_parsed
chr1  | 100 | .  | A   | G   | 60   | PASS   | A|G        | G|G        | {...}
```

## Ejemplo Completo

```python
from pyMut.input import read_vcf
import logging

# Configurar logging para ver el progreso
logging.basicConfig(level=logging.INFO)

# Cargar datos VCF con todas las opciones
py_mut = read_vcf(
    path="src/pyMut/data/examples/ALL.chr10.vcf.gz",
    fasta="reference/hg38.fasta",
    create_index=True,
    cache_dir="cache/"
)

# Verificar carga exitosa
print(f"Muestras cargadas: {len(py_mut.samples)}")
print(f"Variantes totales: {len(py_mut.data)}")
print(f"Cromosomas únicos: {py_mut.data['CHROM'].unique()}")

# Información sobre genotipos
print(f"Columnas de muestras: {py_mut.samples[:5]}...")  # Primeras 5 muestras

# Verificar metadatos
print(f"Formato fuente: {py_mut.metadata.source_format}")
print(f"Archivo FASTA: {py_mut.metadata.fasta}")
```

## Manejo de Genotipos

La función convierte automáticamente los genotipos VCF al formato alélico de pyMut:

### Genotipos VCF (entrada)
```
FORMAT: GT:DP:GQ
SAMPLE_001: 0/1:30:99    # Heterocigoto
SAMPLE_002: 1/1:25:99    # Homocigoto alternativo
SAMPLE_003: 0/0:35:99    # Homocigoto referencia
```

### Formato pyMut (salida)
```
SAMPLE_001: A|G    # REF|ALT
SAMPLE_002: G|G    # ALT|ALT  
SAMPLE_003: A|A    # REF|REF
```

## Parsing del Campo INFO

El campo INFO se procesa automáticamente y se expande en columnas separadas:

```python
# INFO original: "AC=2;AF=0.5;AN=4;DP=100"
# Se convierte en columnas:
py_mut.data['INFO_AC']  # [2]
py_mut.data['INFO_AF']  # [0.5]
py_mut.data['INFO_AN']  # [4]
py_mut.data['INFO_DP']  # [100]
```

## Sistema de Cache

```python
# Primera carga - crea cache
py_mut1 = read_vcf("large_file.vcf.gz", cache_dir="cache/")

# Segunda carga - usa cache (mucho más rápido)
py_mut2 = read_vcf("large_file.vcf.gz", cache_dir="cache/")
```

## Optimización y Rendimiento

### Con pyarrow (recomendado)
```python
# Instalar pyarrow para mejor rendimiento
# pip install pyarrow

py_mut = read_vcf("variants.vcf.gz")  # Usa pyarrow automáticamente
```

### Sin pyarrow (fallback)
```python
# Si pyarrow no está disponible, usa cyvcf2 o pandas
# El rendimiento será menor pero funcional
```

## Manejo de Errores

### FileNotFoundError
```python
try:
    py_mut = read_vcf("archivo_inexistente.vcf")
except FileNotFoundError:
    print("❌ Archivo VCF no encontrado")
```

### Archivo VCF malformado
```python
try:
    py_mut = read_vcf("archivo_corrupto.vcf")
except ValueError as e:
    print(f"❌ Error en formato VCF: {e}")
```

### Problemas de memoria
```python
# Para archivos muy grandes
import psutil
print(f"Memoria disponible: {psutil.virtual_memory().available / 1e9:.1f} GB")

# Considerar filtrar por región primero
py_mut = read_vcf("huge_file.vcf.gz")
filtered = py_mut.region("chr1", 1000000, 2000000)
```

## Indexación Tabix

```python
# Crear índice para consultas rápidas
py_mut = read_vcf("variants.vcf.gz", create_index=True)

# El índice se guarda como variants.vcf.gz.tbi
# Permite consultas rápidas por región posteriormente
```

## Solución de Problemas

### Archivo muy grande
```python
# Monitorear progreso con logging
import logging
logging.basicConfig(level=logging.INFO)

py_mut = read_vcf("archivo_grande.vcf.gz")
```

### Muestras con nombres especiales
```python
# Los nombres de muestras se normalizan automáticamente
# Caracteres especiales se reemplazan por guiones bajos
```

### Campos INFO complejos
```python
# Los campos INFO anidados se procesan automáticamente
# Arrays se convierten en listas de Python
```

## Compatibilidad

- **Versiones VCF**: Compatible con VCF 4.0, 4.1, 4.2, 4.3
- **Compresión**: bgzip, gzip estándar
- **Índices**: tabix (.tbi), CSI (.csi)
- **Genotipos**: Diploides y haploides
- **Campos especiales**: Soporte completo para INFO, FORMAT, FILTER