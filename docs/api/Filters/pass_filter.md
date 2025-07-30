# pass_filter - Filtro de Calidad PASS

El método **pass_filter** permite verificar si una variante específica tiene el estado de calidad "PASS" en el campo FILTER de datos VCF.

## ¿Qué es pass_filter?

Es un método que verifica si una variante específica (definida por cromosoma, posición, alelo de referencia y alelo alternativo) tiene el valor "PASS" en el campo FILTER, indicando que pasó todos los filtros de calidad aplicados durante el proceso de llamada de variantes.

## Características Principales

- **Verificación específica**: Busca una variante exacta por coordenadas y alelos
- **Optimización con pyarrow**: Utiliza pyarrow para consultas rápidas cuando está disponible
- **Fallback automático**: Si pyarrow falla, usa operaciones pandas estándar
- **Normalización de cromosomas**: Maneja diferentes formatos de nomenclatura cromosómica
- **Logging detallado**: Proporciona información sobre el proceso de búsqueda
- **Manejo de duplicados**: Detecta y maneja registros duplicados

## Uso Básico

```python
from pyMut.input import read_vcf

# Cargar datos VCF
py_mut = read_vcf("variants.vcf.gz")

# Verificar si una variante específica tiene FILTER=PASS
is_pass = py_mut.pass_filter(
    chrom="chr17",
    pos=7577121,
    ref="C",
    alt="T"
)

print(f"¿La variante pasó los filtros de calidad? {is_pass}")
```

## Parámetros

### chrom (str) [requerido]
- **Descripción**: Cromosoma de la variante a verificar
- **Formatos aceptados**: `"chr17"`, `"17"`, `"X"`, `"Y"`, `"chrX"`, `"chrY"`
- **Normalización**: Se normaliza automáticamente al formato estándar
- **Ejemplo**: `"chr17"`

### pos (int) [requerido]
- **Descripción**: Posición genómica de la variante
- **Coordenadas**: Basadas en 1 (estándar VCF)
- **Ejemplo**: `7577121`

### ref (str) [requerido]
- **Descripción**: Alelo de referencia
- **Formato**: Secuencia de nucleótidos (A, T, G, C, N, -)
- **Ejemplo**: `"C"`

### alt (str) [requerido]
- **Descripción**: Alelo alternativo
- **Formato**: Secuencia de nucleótidos (A, T, G, C, N, -)
- **Ejemplo**: `"T"`

## Valor de Retorno

Retorna un **boolean**:
- **True**: La variante existe y tiene FILTER="PASS"
- **False**: La variante no existe o no tiene FILTER="PASS"

## Ejemplos Detallados

### Verificación de Variantes Específicas

```python
from pyMut.input import read_vcf
import logging

# Configurar logging para ver detalles
logging.basicConfig(level=logging.INFO)

# Cargar datos VCF
py_mut = read_vcf("src/pyMut/data/examples/ALL.chr10.vcf.gz")

# Lista de variantes a verificar
variantes_interes = [
    {"chrom": "chr10", "pos": 60523, "ref": "T", "alt": "G"},
    {"chrom": "10", "pos": 60803, "ref": "T", "alt": "C"},
    {"chrom": "chr10", "pos": 61023, "ref": "C", "alt": "T"},
]

print("=== Verificación de Calidad de Variantes ===")
for i, var in enumerate(variantes_interes, 1):
    print(f"\nVariante {i}: {var['chrom']}:{var['pos']} {var['ref']}>{var['alt']}")
    
    try:
        is_pass = py_mut.pass_filter(
            chrom=var['chrom'],
            pos=var['pos'],
            ref=var['ref'],
            alt=var['alt']
        )
        
        status = "✅ PASS" if is_pass else "❌ NO PASS"
        print(f"Estado de calidad: {status}")
        
    except Exception as e:
        print(f"❌ Error: {e}")
```

### Verificación Masiva de Variantes

```python
# Verificar múltiples variantes de una lista
def verificar_calidad_variantes(py_mut, lista_variantes):
    """
    Verifica la calidad de una lista de variantes
    """
    resultados = []
    
    for var in lista_variantes:
        try:
            is_pass = py_mut.pass_filter(
                chrom=var['chrom'],
                pos=var['pos'],
                ref=var['ref'],
                alt=var['alt']
            )
            
            resultados.append({
                'variante': f"{var['chrom']}:{var['pos']} {var['ref']}>{var['alt']}",
                'pass_filter': is_pass,
                'status': 'success'
            })
            
        except Exception as e:
            resultados.append({
                'variante': f"{var['chrom']}:{var['pos']} {var['ref']}>{var['alt']}",
                'pass_filter': False,
                'status': f'error: {e}'
            })
    
    return resultados

# Ejemplo de uso
variantes_candidatas = [
    {"chrom": "chr10", "pos": 60523, "ref": "T", "alt": "G"},
    {"chrom": "chr10", "pos": 60803, "ref": "T", "alt": "C"},
    {"chrom": "chr10", "pos": 61023, "ref": "C", "alt": "T"},
]

resultados = verificar_calidad_variantes(py_mut, variantes_candidatas)

# Mostrar resultados
print("\n=== Resumen de Verificación ===")
for resultado in resultados:
    print(f"{resultado['variante']}: {resultado['pass_filter']} ({resultado['status']})")
```

## Optimización con PyArrow

```python
# El método utiliza pyarrow automáticamente para mejor rendimiento
try:
    # Conversión automática a tipos pyarrow para optimización
    is_pass = py_mut.pass_filter("chr10", 60523, "T", "G")
    print("✅ Verificación optimizada con pyarrow")
except ImportError:
    print("⚠️ pyarrow no disponible, usando pandas estándar")
except Exception as e:
    print(f"⚠️ pyarrow falló ({e}), usando pandas estándar")
```

## Manejo de Diferentes Formatos de Cromosoma

```python
# El método normaliza automáticamente los formatos de cromosoma
variante_chr = py_mut.pass_filter("chr10", 60523, "T", "G")
variante_num = py_mut.pass_filter("10", 60523, "T", "G")

# Ambas consultas son equivalentes
print(f"Con 'chr10': {variante_chr}")
print(f"Con '10': {variante_num}")
print(f"¿Son iguales? {variante_chr == variante_num}")
```

## Casos de Uso Comunes

### Validación de Variantes de Interés

```python
# Verificar variantes reportadas en literatura
variantes_literatura = [
    # Variante en TP53 reportada como patogénica
    {"chrom": "chr17", "pos": 7577121, "ref": "C", "alt": "T"},
    # Variante en BRCA1 de significado clínico
    {"chrom": "chr17", "pos": 43094077, "ref": "A", "alt": "C"},
]

print("=== Validación de Variantes de Literatura ===")
for var in variantes_literatura:
    is_pass = py_mut.pass_filter(
        chrom=var['chrom'],
        pos=var['pos'],
        ref=var['ref'],
        alt=var['alt']
    )
    
    if is_pass:
        print(f"✅ {var['chrom']}:{var['pos']} {var['ref']}>{var['alt']} - Calidad PASS")
    else:
        print(f"⚠️ {var['chrom']}:{var['pos']} {var['ref']}>{var['alt']} - No PASS o no encontrada")
```

### Control de Calidad en Pipeline

```python
def filtrar_variantes_alta_calidad(py_mut, lista_variantes):
    """
    Filtra variantes que pasan los controles de calidad
    """
    variantes_pass = []
    variantes_fail = []
    
    for var in lista_variantes:
        is_pass = py_mut.pass_filter(
            chrom=var['chrom'],
            pos=var['pos'],
            ref=var['ref'],
            alt=var['alt']
        )
        
        if is_pass:
            variantes_pass.append(var)
        else:
            variantes_fail.append(var)
    
    return variantes_pass, variantes_fail

# Ejemplo de uso en pipeline
todas_variantes = [
    {"chrom": "chr10", "pos": 60523, "ref": "T", "alt": "G"},
    {"chrom": "chr10", "pos": 60803, "ref": "T", "alt": "C"},
    {"chrom": "chr10", "pos": 61023, "ref": "C", "alt": "T"},
]

pass_vars, fail_vars = filtrar_variantes_alta_calidad(py_mut, todas_variantes)

print(f"Variantes que pasan filtros: {len(pass_vars)}")
print(f"Variantes que fallan filtros: {len(fail_vars)}")
```

### Análisis de Calidad por Región

```python
# Verificar calidad de variantes en una región específica
def analizar_calidad_region(py_mut, chrom, start, end):
    """
    Analiza la calidad de variantes en una región genómica
    """
    # Primero filtrar por región
    region_data = py_mut.region(chrom, start, end)
    
    if len(region_data.data) == 0:
        print(f"No hay variantes en {chrom}:{start}-{end}")
        return
    
    # Verificar calidad de cada variante en la región
    total_variantes = 0
    pass_variantes = 0
    
    for _, row in region_data.data.iterrows():
        total_variantes += 1
        
        is_pass = py_mut.pass_filter(
            chrom=row['CHROM'],
            pos=row['POS'],
            ref=row['REF'],
            alt=row['ALT']
        )
        
        if is_pass:
            pass_variantes += 1
    
    porcentaje_pass = (pass_variantes / total_variantes) * 100
    
    print(f"=== Análisis de Calidad: {chrom}:{start}-{end} ===")
    print(f"Total de variantes: {total_variantes}")
    print(f"Variantes PASS: {pass_variantes}")
    print(f"Porcentaje PASS: {porcentaje_pass:.1f}%")

# Ejemplo de uso
analizar_calidad_region(py_mut, "chr10", 60000, 70000)
```

## Manejo de Errores

### Variante no encontrada

```python
try:
    # Buscar variante que probablemente no existe
    is_pass = py_mut.pass_filter("chr99", 999999999, "A", "T")
    print(f"Resultado: {is_pass}")  # Será False
except KeyError as e:
    print(f"❌ Error de columnas: {e}")
```

### Datos malformados

```python
try:
    # Verificar con datos inválidos
    is_pass = py_mut.pass_filter("chr10", "invalid_pos", "A", "T")
except (ValueError, TypeError) as e:
    print(f"❌ Error en formato de datos: {e}")
```

### Columnas faltantes

```python
# El método verifica automáticamente las columnas requeridas
required_columns = ["CHROM", "POS", "REF", "ALT", "FILTER"]

# Si faltan columnas, se lanza KeyError con información detallada
try:
    is_pass = py_mut.pass_filter("chr10", 60523, "T", "G")
except KeyError as e:
    print(f"❌ Columnas faltantes: {e}")
```

## Rendimiento y Optimización

### Para consultas múltiples

```python
import time

# Medir rendimiento de consultas múltiples
variantes_test = [
    {"chrom": "chr10", "pos": 60523, "ref": "T", "alt": "G"},
    {"chrom": "chr10", "pos": 60803, "ref": "T", "alt": "C"},
    {"chrom": "chr10", "pos": 61023, "ref": "C", "alt": "T"},
] * 100  # Repetir 100 veces para test de rendimiento

start_time = time.time()

resultados = []
for var in variantes_test:
    is_pass = py_mut.pass_filter(
        chrom=var['chrom'],
        pos=var['pos'],
        ref=var['ref'],
        alt=var['alt']
    )
    resultados.append(is_pass)

end_time = time.time()

print(f"Verificadas {len(variantes_test)} variantes en {end_time - start_time:.2f} segundos")
print(f"Promedio: {(end_time - start_time) / len(variantes_test) * 1000:.2f} ms por variante")
```

### Optimización de memoria

```python
# Para datasets muy grandes, considerar filtrar primero por región
# antes de verificar variantes específicas

# Menos eficiente: buscar en todo el dataset
is_pass = py_mut.pass_filter("chr10", 60523, "T", "G")

# Más eficiente: filtrar región primero, luego buscar
region_subset = py_mut.region("chr10", 60000, 61000)
is_pass = region_subset.pass_filter("chr10", 60523, "T", "G")
```