# Tissue Expression Filter - Filtro de Expresión Tisular

El método **filter_by_tissue_expression** permite filtrar datos de mutaciones basándose en la expresión génica en tejidos específicos utilizando datos de consenso de RNA de cáncer.

## ¿Qué es el Filtro de Expresión Tisular?

Este filtro permite identificar y filtrar mutaciones en genes que están suficientemente expresados (o no expresados) en tejidos/tipos de cáncer específicos. Utiliza datos de consenso de RNA de cáncer para determinar si un gen está "activado" en un tejido particular según un umbral de expresión definido.

## Características Principales

- **Filtrado por múltiples tejidos**: Permite especificar múltiples tejidos con umbrales independientes
- **Detección automática de columnas**: Identifica automáticamente columnas de símbolos génicos
- **Caching inteligente**: Los datos de expresión se cargan una vez y se mantienen en caché
- **Filtrado bidireccional**: Puede filtrar genes expresados o no expresados
- **Códigos TCGA**: Utiliza códigos estándar de TCGA para tipos de cáncer
- **Validación robusta**: Validación completa de parámetros de entrada

## Uso Básico

```python
from pyMut.io import read_maf

# Cargar datos
py_mut = read_maf("mutations.maf")

# Filtrar genes expresados en cáncer de vejiga (umbral 5)
filtered_mut = py_mut.filter_by_tissue_expression([('BLCA', 5)])

# Filtrar genes NO expresados en adenocarcinoma de pulmón
not_expressed_mut = py_mut.filter_by_tissue_expression(
    [('LUAD', 4)], 
    keep_expressed=False
)

# Filtrar genes expresados en múltiples tejidos con diferentes umbrales
multi_tissue_mut = py_mut.filter_by_tissue_expression([
    ('BLCA', 5),    # Cáncer de vejiga, umbral 5
    ('BRCA', 3),    # Cáncer de mama, umbral 3
    ('LUAD', 4),    # Adenocarcinoma de pulmón, umbral 4
    ('COAD', 6)     # Adenocarcinoma de colon, umbral 6
])

print(f"Mutaciones originales: {len(py_mut.data)}")
print(f"Mutaciones filtradas: {len(filtered_mut.data)}")
```

## Parámetros

### tissues (List[Tuple[str, float]], requerido)
- **Descripción**: Lista de tuplas con especificaciones de tejido y umbral
- **Formato**: `[('código_tejido', umbral), ...]`
- **Códigos TCGA**: Utiliza códigos estándar como 'BLCA', 'BRCA', 'LUAD', etc.
- **Ejemplo**: `[('BLCA', 5), ('BRCA', 3)]`

### keep_expressed (bool, default=True)
- **Descripción**: Determina qué mutaciones mantener
- **True**: Mantiene genes expresados en al menos uno de los tejidos especificados
- **False**: Mantiene genes NO expresados en ninguno de los tejidos especificados

## Valor de Retorno

Retorna un nuevo objeto **PyMutation** con los datos filtrados según los criterios de expresión tisular especificados.

```python
# El objeto retornado mantiene la misma estructura
filtered_mut = py_mut.filter_by_tissue_expression([('BLCA', 5)])
print(type(filtered_mut))  # <class 'pyMut.core.PyMutation'>
```

## Códigos de Tejido TCGA Soportados

El filtro utiliza códigos estándar de TCGA para diferentes tipos de cáncer:

```python
# Ejemplos de códigos TCGA comunes
tissue_codes = {
    'BLCA': 'Bladder Urothelial Carcinoma',
    'BRCA': 'Breast Invasive Carcinoma', 
    'COAD': 'Colon Adenocarcinoma',
    'LUAD': 'Lung Adenocarcinoma',
    'LUSC': 'Lung Squamous Cell Carcinoma',
    'PRAD': 'Prostate Adenocarcinoma',
    'THCA': 'Thyroid Carcinoma',
    'LIHC': 'Liver Hepatocellular Carcinoma',
    'STAD': 'Stomach Adenocarcinoma',
    'SKCM': 'Skin Cutaneous Melanoma'
    # ... y muchos más
}
```

## Función Auxiliar: tissue_expression

Además del método de filtrado, está disponible la función auxiliar `tissue_expression` para verificaciones individuales:

```python
from pyMut.filters.tissue_expression import tissue_expression
import pandas as pd

# Verificar expresión usando símbolo génico directamente
is_expressed = tissue_expression("TSPAN6", ["BLCA", 5])
print(f"TSPAN6 expresado en BLCA (umbral 5): {is_expressed}")

# Verificar expresión usando una fila de datos
row = pd.Series({'Hugo_Symbol': 'TSPAN6', 'Chromosome': 'X'})
is_expressed = tissue_expression(row, ["BRCA", 10])
print(f"TSPAN6 expresado en BRCA (umbral 10): {is_expressed}")
```

### Parámetros de tissue_expression

- **data**: `Union[str, pd.Series]` - Símbolo génico o fila de datos
- **tissue**: `List[Union[str, float]]` - `[código_tejido, umbral]`

## Ejemplos Avanzados

### Filtrado Condicional por Tipo de Cáncer

```python
# Filtrar mutaciones relevantes para cáncer de mama
breast_cancer_mut = py_mut.filter_by_tissue_expression([
    ('BRCA', 5)  # Genes expresados en cáncer de mama
])

# Filtrar genes silenciados en cáncer de pulmón
lung_silenced_mut = py_mut.filter_by_tissue_expression([
    ('LUAD', 3),
    ('LUSC', 3)
], keep_expressed=False)
```

### Análisis Multi-Tejido

```python
# Genes expresados en al menos uno de múltiples cánceres digestivos
digestive_cancers_mut = py_mut.filter_by_tissue_expression([
    ('COAD', 4),    # Colon
    ('STAD', 4),    # Estómago  
    ('LIHC', 5),    # Hígado
    ('PAAD', 6)     # Páncreas
])

print(f"Mutaciones en genes expresados en cánceres digestivos: {len(digestive_cancers_mut.data)}")
```

### Combinación con Otros Filtros

```python
# Combinar filtro de expresión tisular con otros filtros
filtered_mut = (py_mut
    .filter_by_tissue_expression([('BRCA', 5)])
    .filter_by_pass()
    .filter_by_chromosome(['1', '2', '3']))

print(f"Mutaciones después de filtros combinados: {len(filtered_mut.data)}")
```

## Manejo de Errores

El método incluye validación robusta y manejo de errores:

```python
try:
    # Error: lista vacía
    py_mut.filter_by_tissue_expression([])
except ValueError as e:
    print(f"Error: {e}")

try:
    # Error: formato incorrecto de tupla
    py_mut.filter_by_tissue_expression([('BLCA', 5, 'extra')])
except ValueError as e:
    print(f"Error: {e}")

try:
    # Error: código de tejido no es string
    py_mut.filter_by_tissue_expression([(123, 5)])
except ValueError as e:
    print(f"Error: {e}")
```

## Consideraciones de Rendimiento

- **Caching**: Los datos de expresión se cargan una vez y se mantienen en caché
- **Memoria**: El filtrado crea una copia de los datos, no modifica el objeto original
- **Escalabilidad**: Eficiente para datasets grandes gracias al caching inteligente

## Casos de Uso Comunes

1. **Análisis específico de cáncer**: Filtrar mutaciones relevantes para un tipo de cáncer específico
2. **Genes silenciados**: Identificar mutaciones en genes no expresados (posibles supresores tumorales)
3. **Análisis comparativo**: Comparar patrones de mutación entre diferentes tipos de cáncer
4. **Priorización de variantes**: Enfocar análisis en genes activos en tejidos relevantes

## Notas Técnicas

- Los datos de expresión provienen de `rna_cancer_consensus.json`
- Se detectan automáticamente columnas como 'Hugo_Symbol', 'Gene_Symbol', etc.
- El filtrado preserva todos los metadatos del objeto PyMutation original
- Compatible con datos en formato MAF y VCF convertidos