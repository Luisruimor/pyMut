# create_subset - Extracción de Subconjuntos VCF

La función **extract_vcf_subset** permite crear archivos VCF más pequeños a partir de archivos VCF grandes, facilitando el trabajo con subconjuntos de datos para pruebas y análisis específicos.

## ¿Qué es extract_vcf_subset?

Es una función utilitaria que extrae un subconjunto de variantes de un archivo VCF sin requerir herramientas externas como tabix. Funciona tanto con archivos comprimidos (.gz) como sin comprimir.

## Características Principales

- **Sin dependencias externas**: No requiere tabix u otras herramientas de indexación
- **Soporte para compresión**: Maneja archivos .vcf y .vcf.gz automáticamente
- **Filtrado flexible**: Por cromosoma, posición y número máximo de variantes
- **Preservación de headers**: Mantiene todos los headers VCF originales
- **Logging detallado**: Proporciona información sobre el progreso y estadísticas
- **Manejo de memoria eficiente**: Procesa archivos línea por línea

## Uso Básico

```python
from pyMut.utils.create_subset import extract_vcf_subset

# Extraer primeras 10,000 variantes del cromosoma 10
success = extract_vcf_subset(
    input_vcf_path="large_file.vcf.gz",
    output_vcf_path="subset_chr10.vcf",
    chromosome="10",
    max_variants=10000
)

if success:
    print("✅ Subconjunto creado exitosamente")
else:
    print("❌ Error al crear subconjunto")
```

## Parámetros

### input_vcf_path (str | Path) [requerido]
- **Descripción**: Ruta al archivo VCF de entrada
- **Formatos soportados**: `.vcf`, `.vcf.gz`
- **Ejemplo**: `"data/large_variants.vcf.gz"`

### output_vcf_path (str | Path) [requerido]
- **Descripción**: Ruta donde guardar el archivo VCF de subconjunto
- **Formato**: Generalmente `.vcf` (sin comprimir para facilitar inspección)
- **Ejemplo**: `"subset/chr10_10k_variants.vcf"`

### max_variants (int, opcional)
- **Descripción**: Número máximo de variantes a extraer
- **Comportamiento**: Extrae las primeras N variantes que cumplan los criterios
- **Ejemplo**: `50000` para extraer 50,000 variantes

### chromosome (str, opcional)
- **Descripción**: Cromosoma específico a extraer
- **Formatos aceptados**: `"10"`, `"chr10"`, `"X"`, `"chrX"`
- **Normalización**: Maneja automáticamente formatos con y sin prefijo "chr"

### start_pos (int, opcional)
- **Descripción**: Posición inicial para la extracción
- **Uso**: Debe usarse junto con `chromosome`
- **Ejemplo**: `1000000` para comenzar desde la posición 1Mb

### end_pos (int, opcional)
- **Descripción**: Posición final para la extracción
- **Uso**: Debe usarse junto con `chromosome` y `start_pos`
- **Ejemplo**: `2000000` para terminar en la posición 2Mb

## Valor de Retorno

Retorna un **boolean**:
- **True**: Extracción exitosa
- **False**: Error durante la extracción

## Ejemplos Detallados

### Extracción por Cromosoma

```python
from pyMut.utils.create_subset import extract_vcf_subset
import logging

# Configurar logging para ver progreso
logging.basicConfig(level=logging.INFO)

# Extraer todas las variantes del cromosoma X
success = extract_vcf_subset(
    input_vcf_path="src/pyMut/data/examples/ALL.chr10.vcf.gz",
    output_vcf_path="chrX_variants.vcf",
    chromosome="X"
)

print(f"Extracción cromosoma X: {'✅ Exitosa' if success else '❌ Falló'}")
```

### Extracción por Número de Variantes

```python
# Crear archivo de prueba con primeras 1000 variantes
test_subset = extract_vcf_subset(
    input_vcf_path="large_dataset.vcf.gz",
    output_vcf_path="test_1k_variants.vcf",
    max_variants=1000
)

# Crear archivo de desarrollo con 50k variantes
dev_subset = extract_vcf_subset(
    input_vcf_path="large_dataset.vcf.gz",
    output_vcf_path="dev_50k_variants.vcf",
    max_variants=50000
)

print(f"Archivo de prueba: {'✅' if test_subset else '❌'}")
print(f"Archivo de desarrollo: {'✅' if dev_subset else '❌'}")
```

### Extracción por Región Específica

```python
# Extraer región específica del cromosoma 17 (región TP53)
tp53_region = extract_vcf_subset(
    input_vcf_path="whole_genome.vcf.gz",
    output_vcf_path="tp53_region.vcf",
    chromosome="17",
    start_pos=7571720,
    end_pos=7590868,
    max_variants=10000  # Limitar a 10k variantes máximo
)

print(f"Región TP53 extraída: {'✅' if tp53_region else '❌'}")
```

### Pipeline de Creación de Subconjuntos

```python
def crear_subconjuntos_desarrollo(input_file, output_dir):
    """
    Crea múltiples subconjuntos para diferentes propósitos de desarrollo
    """
    from pathlib import Path
    
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    subconjuntos = {
        "test_small": {"max_variants": 1000, "desc": "Pruebas rápidas"},
        "test_medium": {"max_variants": 10000, "desc": "Pruebas intermedias"},
        "dev_large": {"max_variants": 100000, "desc": "Desarrollo completo"},
        "chr17_only": {"chromosome": "17", "desc": "Solo cromosoma 17"},
        "chrX_only": {"chromosome": "X", "desc": "Solo cromosoma X"}
    }
    
    resultados = {}
    
    for nombre, config in subconjuntos.items():
        output_file = output_dir / f"{nombre}.vcf"
        
        print(f"\nCreando {nombre} ({config['desc']})...")
        
        success = extract_vcf_subset(
            input_vcf_path=input_file,
            output_vcf_path=output_file,
            **{k: v for k, v in config.items() if k != 'desc'}
        )
        
        resultados[nombre] = success
        status = "✅ Exitoso" if success else "❌ Falló"
        print(f"{nombre}: {status}")
    
    return resultados

# Ejemplo de uso
resultados = crear_subconjuntos_desarrollo(
    "large_dataset.vcf.gz",
    "subsets/"
)

print(f"\nResumen: {sum(resultados.values())}/{len(resultados)} subconjuntos creados")
```

## Casos de Uso Comunes

### Desarrollo y Pruebas

```python
# Crear archivo pequeño para pruebas unitarias
extract_vcf_subset(
    "production_data.vcf.gz",
    "test_data.vcf",
    max_variants=500
)

# Crear archivo mediano para pruebas de integración
extract_vcf_subset(
    "production_data.vcf.gz",
    "integration_test.vcf",
    max_variants=5000
)
```

### Análisis de Cromosomas Específicos

```python
# Cromosomas de interés en cáncer
cancer_chromosomes = ["17", "13", "3", "7", "12"]

for chrom in cancer_chromosomes:
    extract_vcf_subset(
        "whole_genome.vcf.gz",
        f"chr{chrom}_variants.vcf",
        chromosome=chrom,
        max_variants=20000
    )
```

### Creación de Datasets de Entrenamiento

```python
# Dataset pequeño para pruebas de algoritmos
extract_vcf_subset(
    "training_data.vcf.gz",
    "algorithm_test.vcf",
    max_variants=1000
)

# Dataset mediano para validación
extract_vcf_subset(
    "training_data.vcf.gz",
    "validation_set.vcf",
    max_variants=10000
)
```

## Información de Logging

La función proporciona logging detallado durante la ejecución:

```
INFO - Extracting VCF subset
INFO - Input: large_file.vcf.gz
INFO - Output: subset.vcf
INFO - Max variants: 10000
INFO - Chromosome filter: 10
INFO - Position range: None-None
INFO - Processed 10000 variants...
INFO - VCF subset extraction completed successfully
INFO - Variants extracted: 10000
INFO - Input file size: 1250.45 MB
INFO - Output file size: 12.34 MB
INFO - Size reduction: 101.3x
```

## Manejo de Errores

### Archivo de entrada no existe

```python
try:
    success = extract_vcf_subset("nonexistent.vcf", "output.vcf")
    if not success:
        print("❌ Error: Archivo de entrada no encontrado")
except Exception as e:
    print(f"❌ Error inesperado: {e}")
```

### Archivo VCF malformado

```python
# La función maneja automáticamente archivos malformados
# Registra warnings y continúa con las líneas válidas
success = extract_vcf_subset(
    "potentially_corrupted.vcf",
    "cleaned_output.vcf",
    max_variants=1000
)
```

### Problemas de permisos

```python
import os
from pathlib import Path

output_path = Path("protected_dir/output.vcf")

# Verificar permisos antes de la extracción
if not os.access(output_path.parent, os.W_OK):
    print("❌ Sin permisos de escritura en directorio de destino")
else:
    success = extract_vcf_subset("input.vcf", output_path)
```

## Optimización y Rendimiento

### Para archivos muy grandes

```python
# Procesar por lotes para archivos extremadamente grandes
def extraer_por_lotes(input_file, batch_size=50000):
    """
    Extrae múltiples lotes de un archivo muy grande
    """
    for i in range(5):  # 5 lotes
        start_variant = i * batch_size
        output_file = f"batch_{i+1}.vcf"
        
        # Nota: Esta función no soporta skip de variantes
        # Se necesitaría implementación adicional para esto
        extract_vcf_subset(
            input_file,
            output_file,
            max_variants=batch_size
        )
```

### Monitoreo de recursos

```python
import psutil
import time

def extraer_con_monitoreo(input_file, output_file, **kwargs):
    """
    Extrae subconjunto monitoreando recursos del sistema
    """
    start_time = time.time()
    start_memory = psutil.virtual_memory().used / 1024 / 1024  # MB
    
    print(f"Memoria inicial: {start_memory:.1f} MB")
    
    success = extract_vcf_subset(input_file, output_file, **kwargs)
    
    end_time = time.time()
    end_memory = psutil.virtual_memory().used / 1024 / 1024  # MB
    
    print(f"Tiempo transcurrido: {end_time - start_time:.2f} segundos")
    print(f"Memoria final: {end_memory:.1f} MB")
    print(f"Incremento de memoria: {end_memory - start_memory:.1f} MB")
    
    return success

# Ejemplo de uso
extraer_con_monitoreo(
    "large_file.vcf.gz",
    "subset.vcf",
    chromosome="17",
    max_variants=10000
)
```

## Integración con pyMut

```python
from pyMut.utils.create_subset import extract_vcf_subset
from pyMut.input import read_vcf

# 1. Crear subconjunto para análisis
extract_vcf_subset(
    "large_dataset.vcf.gz",
    "analysis_subset.vcf",
    chromosome="17",
    max_variants=5000
)

# 2. Cargar subconjunto en pyMut
py_mut = read_vcf("analysis_subset.vcf")

# 3. Realizar análisis
summary_fig = py_mut.summary_plot(title="Análisis Cromosoma 17")
summary_fig.savefig("chr17_analysis.png")

print(f"Análisis completado con {len(py_mut.data)} variantes")
```

## Solución de Problemas

### Archivos de salida vacíos

```python
# Verificar que los criterios de filtrado no sean demasiado restrictivos
success = extract_vcf_subset(
    "input.vcf.gz",
    "output.vcf",
    chromosome="99",  # Cromosoma inexistente
    max_variants=1000
)

# Si success es True pero el archivo está vacío,
# verificar los criterios de filtrado
```

### Rendimiento lento

```python
# Para mejorar rendimiento, usar filtros más específicos
# Menos eficiente: extraer todo y luego filtrar
extract_vcf_subset("huge_file.vcf.gz", "all.vcf", max_variants=100000)

# Más eficiente: filtrar por cromosoma primero
extract_vcf_subset(
    "huge_file.vcf.gz", 
    "chr17.vcf", 
    chromosome="17",
    max_variants=10000
)
```