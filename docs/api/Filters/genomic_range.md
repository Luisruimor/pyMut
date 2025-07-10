# Filtros de Rango Genómico - region y gen_region

Los métodos **region** y **gen_region** permiten filtrar datos de PyMutation por ubicación genómica específica, ya sea por coordenadas cromosómicas o por nombre de gen.

## ¿Qué son los Filtros de Rango Genómico?

Son métodos que permiten extraer subconjuntos de mutaciones basados en su ubicación en el genoma, facilitando el análisis de regiones específicas de interés.

## Características Principales

- **Filtrado por coordenadas**: Especifica cromosoma, posición inicial y final
- **Filtrado por gen**: Busca automáticamente las coordenadas de un gen específico
- **Optimización con pyarrow**: Utiliza pyarrow para consultas rápidas cuando está disponible
- **Preservación de metadatos**: Mantiene la información original y registra los filtros aplicados
- **Validación automática**: Verifica que las coordenadas sean válidas
- **Logging detallado**: Proporciona información sobre el proceso de filtrado

## Método region - Filtrado por Coordenadas

### Uso Básico

```python
from pyMut.input import read_maf

# Cargar datos
py_mut = read_maf("mutations.maf")

# Filtrar por región específica
# Cromosoma 17, posiciones 7,500,000 a 7,600,000
tp53_region = py_mut.region("chr17", 7500000, 7600000)

print(f"Mutaciones en región: {len(tp53_region.data)}")
```

### Parámetros de region

#### chrom (str) [requerido]
- **Descripción**: Cromosoma a filtrar
- **Formatos aceptados**: `"chr17"`, `"17"`, `"X"`, `"Y"`, `"chrX"`, `"chrY"`
- **Ejemplo**: `"chr17"`

#### start (int) [requerido]
- **Descripción**: Posición inicial de la región (inclusiva)
- **Coordenadas**: Basadas en 1 (estándar genómico)
- **Ejemplo**: `7500000`

#### end (int) [requerido]
- **Descripción**: Posición final de la región (inclusiva)
- **Coordenadas**: Basadas en 1 (estándar genómico)
- **Ejemplo**: `7600000`

### Ejemplos de region

```python
# Región del gen TP53 (cromosoma 17)
tp53_mutations = py_mut.region("chr17", 7571720, 7590868)

# Región del gen BRCA1 (cromosoma 17)
brca1_mutations = py_mut.region("17", 43044295, 43125483)

# Cromosoma X completo (ejemplo de región grande)
chrx_mutations = py_mut.region("X", 1, 156040895)

# Región específica en cromosoma Y
chry_region = py_mut.region("chrY", 2781479, 2781479)  # Posición específica
```

## Método gen_region - Filtrado por Nombre de Gen

### Uso Básico

```python
# Filtrar por nombre de gen (automático)
tp53_mutations = py_mut.gen_region("TP53")
brca1_mutations = py_mut.gen_region("BRCA1")

print(f"Mutaciones en TP53: {len(tp53_mutations.data)}")
print(f"Mutaciones en BRCA1: {len(brca1_mutations.data)}")
```

### Parámetros de gen_region

#### gen_name (str) [requerido]
- **Descripción**: Nombre del gen a buscar
- **Formato**: Símbolo oficial del gen (HUGO)
- **Ejemplos**: `"TP53"`, `"BRCA1"`, `"EGFR"`, `"KRAS"`

### Base de Datos de Genes

El método utiliza una base de datos interna con información genómica:

```python
# Genes soportados incluyen:
genes_disponibles = [
    "TP53",      # Cromosoma 17: 7,571,720-7,590,868
    "BRCA1",     # Cromosoma 17: 43,044,295-43,125,483
    "BRCA2",     # Cromosoma 13: 32,315,086-32,400,268
    "EGFR",      # Cromosoma 7: 55,019,017-55,211,628
    "KRAS",      # Cromosoma 12: 25,205,246-25,250,929
    "PIK3CA",    # Cromosoma 3: 179,148,114-179,240,093
    # ... y muchos más
]
```

## Ejemplo Completo

```python
from pyMut.input import read_maf
import logging

# Configurar logging para ver detalles
logging.basicConfig(level=logging.INFO)

# Cargar datos TCGA
py_mut = read_maf("src/pyMut/data/examples/tcga_laml.maf.gz")
print(f"Mutaciones totales: {len(py_mut.data)}")

# Análisis de genes específicos
genes_interes = ["TP53", "KRAS", "PIK3CA", "EGFR"]

for gen in genes_interes:
    try:
        # Filtrar por gen
        gen_mutations = py_mut.gen_region(gen)
        print(f"\n=== Análisis de {gen} ===")
        print(f"Mutaciones encontradas: {len(gen_mutations.data)}")
        
        if len(gen_mutations.data) > 0:
            # Tipos de mutaciones más comunes
            tipos = gen_mutations.data['Variant_Classification'].value_counts()
            print(f"Tipos de mutación:")
            for tipo, count in tipos.head(3).items():
                print(f"  - {tipo}: {count}")
            
            # Muestras afectadas
            muestras = gen_mutations.data['Tumor_Sample_Barcode'].nunique()
            print(f"Muestras afectadas: {muestras}")
            
    except Exception as e:
        print(f"❌ Error procesando {gen}: {e}")

# Análisis de región específica
print(f"\n=== Análisis de Región Cromosómica ===")
# Región rica en genes de cáncer en cromosoma 17
chr17_region = py_mut.region("chr17", 7000000, 8000000)
print(f"Mutaciones en chr17:7M-8M: {len(chr17_region.data)}")

# Genes en esta región
if len(chr17_region.data) > 0:
    genes_region = chr17_region.data['Hugo_Symbol'].value_counts()
    print("Genes más mutados en la región:")
    for gen, count in genes_region.head(5).items():
        print(f"  - {gen}: {count} mutaciones")
```

## Optimización con PyArrow

```python
# El filtrado utiliza pyarrow automáticamente para mejor rendimiento
try:
    # Conversión automática a tipos pyarrow
    filtered = py_mut.region("chr17", 7500000, 7600000)
    print("✅ Filtrado optimizado con pyarrow")
except ImportError:
    print("⚠️ pyarrow no disponible, usando pandas estándar")
```

## Validación de Coordenadas

```python
# Validaciones automáticas realizadas:
try:
    # Coordenadas válidas
    valid_region = py_mut.region("chr1", 1000000, 2000000)
    
    # Error: start > end
    invalid_region = py_mut.region("chr1", 2000000, 1000000)
    
except ValueError as e:
    print(f"❌ Coordenadas inválidas: {e}")

# Error: cromosoma inexistente
try:
    invalid_chr = py_mut.region("chr99", 1000000, 2000000)
except KeyError as e:
    print(f"❌ Cromosoma no encontrado: {e}")
```

## Manejo de Metadatos

```python
# Los filtros se registran automáticamente
original = py_mut
filtered = py_mut.region("chr17", 7500000, 7600000)

print("Filtros aplicados:")
for filtro in filtered.metadata.filters:
    print(f"  - {filtro}")

# Ejemplo de salida:
# - region:chr17:7500000-7600000
```

## Casos de Uso Comunes

### Análisis de Genes Candidatos

```python
# Lista de genes de interés oncológico
cancer_genes = ["TP53", "KRAS", "PIK3CA", "EGFR", "BRAF", "APC"]

# Analizar cada gen individualmente
gene_analysis = {}
for gene in cancer_genes:
    mutations = py_mut.gen_region(gene)
    gene_analysis[gene] = {
        'mutations': len(mutations.data),
        'samples': mutations.data['Tumor_Sample_Barcode'].nunique() if len(mutations.data) > 0 else 0
    }

# Resumen
for gene, stats in gene_analysis.items():
    print(f"{gene}: {stats['mutations']} mutaciones en {stats['samples']} muestras")
```

### Análisis de Regiones Cromosómicas

```python
# Análisis de brazo cromosómico
# Cromosoma 17p (brazo corto)
chr17p = py_mut.region("chr17", 1, 22300000)

# Cromosoma 17q (brazo largo)  
chr17q = py_mut.region("chr17", 22300001, 83257441)

print(f"Mutaciones en 17p: {len(chr17p.data)}")
print(f"Mutaciones en 17q: {len(chr17q.data)}")
```

### Análisis de Hotspots

```python
# Regiones hotspot conocidas
hotspots = {
    "TP53_DBD": ("chr17", 7571720, 7590868),      # Dominio de unión a DNA
    "KRAS_G12": ("chr12", 25245274, 25245276),    # Codón 12
    "PIK3CA_E545": ("chr3", 179218303, 179218305) # Codón 545
}

for name, (chrom, start, end) in hotspots.items():
    hotspot_muts = py_mut.region(chrom, start, end)
    print(f"{name}: {len(hotspot_muts.data)} mutaciones")
```

## Combinación con Otros Filtros

```python
# Combinar filtros genómicos con filtros de muestra
# 1. Filtrar por región
tp53_region = py_mut.region("chr17", 7571720, 7590868)

# 2. Filtrar por muestras específicas
specific_samples = ["TCGA-AB-2802", "TCGA-AB-2803"]
tp53_samples = tp53_region.filter_by_chrom_sample(sample=specific_samples)

print(f"Mutaciones TP53 en muestras específicas: {len(tp53_samples.data)}")
```

## Rendimiento y Optimización

### Para datasets grandes
```python
# Filtrar primero por cromosoma para reducir datos
chr17_only = py_mut.filter_by_chrom_sample(chrom="chr17")

# Luego aplicar filtro de región más específico
tp53_mutations = chr17_only.region("chr17", 7571720, 7590868)
```

### Monitoreo de rendimiento
```python
import time

start_time = time.time()
filtered = py_mut.region("chr17", 7500000, 7600000)
end_time = time.time()

print(f"Filtrado completado en {end_time - start_time:.2f} segundos")
print(f"Mutaciones filtradas: {len(filtered.data)}")
```