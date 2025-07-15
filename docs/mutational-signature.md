# Análisis de Firmas Mutacionales

## Descripción General

El análisis de firmas mutacionales es una técnica poderosa para identificar patrones de mutación recurrentes en datos genómicos del cáncer. pyMut implementa un análisis completo que extrae firmas mutacionales usando Factorización de Matrices No Negativas (NMF) y genera múltiples visualizaciones para su interpretación.

## Componentes de la Visualización

La visualización de firmas mutacionales consta de 5 componentes principales:

### A. Perfiles de Firmas (Signature Profiles)
Gráficos de barras que muestran la distribución de las 96 categorías de sustitución trinucleotídica para cada firma identificada. Las barras están coloreadas según el tipo de sustitución (C>A, C>G, C>T, T>A, T>C, T>G).

### B. Similitud del Coseno (Cosine Similarity)
Mapa de calor que muestra la similitud entre las firmas identificadas y las firmas de referencia COSMIC, permitiendo la identificación de procesos mutacionales conocidos.

### C. Contribución por Muestra (Sample Contributions)
Mapa de calor que visualiza la contribución de cada firma a cada muestra individual.

### D. Distribución Relativa (Relative Distribution)
Gráfico de barras apiladas que muestra la proporción relativa de cada firma en cada muestra.

### E. Proporción Global (Overall Proportions)
Gráfico de donut que muestra la contribución total de cada firma en todo el conjunto de datos.

## Uso Básico

```python
from pyMut import PyMutation
import pandas as pd

# Configurar alta calidad
PyMutation.configure_high_quality_plots()

# Cargar datos
data = pd.read_csv("mutations.maf", sep="\t")
py_mut = PyMutation(data)

# Análisis básico con 3 firmas
fig = py_mut.mutational_signature_plot(n_signatures=3)
fig.savefig("mutational_signatures.png")
```

## Parámetros

### Parámetros principales:

- **n_signatures** (int, default=3): Número de firmas a extraer
- **sample_column** (str): Columna con identificadores de muestra
- **ref_column** (str, default="REF"): Columna con alelos de referencia
- **alt_column** (str, default="ALT"): Columna con alelos alternativos
- **context_column** (str, opcional): Columna con contexto trinucleotídico
- **cosmic_signatures** (DataFrame, opcional): Firmas COSMIC de referencia
- **figsize** (tuple, default=(20, 24)): Tamaño de la figura
- **title** (str): Título principal del análisis
- **show_interactive** (bool, default=False): Mostrar interactivamente

## Ejemplos Avanzados

### Con firmas COSMIC

```python
# Cargar firmas COSMIC
cosmic_df = pd.read_csv("COSMIC_signatures.tsv", sep='\t', index_col=0)

# Análisis con comparación COSMIC
fig = py_mut.mutational_signature_plot(
    n_signatures=4,
    cosmic_signatures=cosmic_df,
    figsize=(24, 28),
    title="Análisis con Comparación COSMIC"
)
```

### Con contexto específico

```python
# Especificar columna de contexto
fig = py_mut.mutational_signature_plot(
    n_signatures=3,
    context_column='Reference_Context',
    title="Análisis con Contexto Trinucleotídico"
)
```

### Análisis interactivo

```python
# Visualización interactiva
fig = py_mut.mutational_signature_plot(
    n_signatures=5,
    show_interactive=True
)
```

## Requisitos de Datos

### Columnas esenciales:
- **REF**: Alelo de referencia (A, C, G, T)
- **ALT**: Alelo alternativo (A, C, G, T)
- **Tumor_Sample_Barcode**: Identificador de muestra

### Columnas recomendadas:
- **Reference_Context**: Contexto trinucleotídico (ej: "ACG", "TTG")
- **Hugo_Symbol**: Símbolo del gen
- **Variant_Classification**: Clasificación de la variante

## Interpretación

### Perfiles de Firma
- Cada firma representa un proceso mutacional distinto
- Los picos altos indican contextos preferenciales
- Los colores agrupan tipos de sustitución similares

### Similitud COSMIC
- Valores cercanos a 1 indican alta similitud
- Permite identificar procesos conocidos (ej: envejecimiento, tabaco, UV)

### Contribución por Muestra
- Identifica qué firmas dominan en cada muestra
- Útil para clasificación de subtipos tumorales

## Notas Técnicas

### Extracción de Firmas
- Utiliza sklearn.decomposition.NMF
- Las matrices se normalizan automáticamente
- Se añade un valor pequeño (1e-6) para evitar divisiones por cero

### Contexto Trinucleotídico
- Si no se especifica, pyMut intenta detectar la columna automáticamente
- Busca columnas con nombres como 'context', 'trinuc', 'reference_context'
- El contexto debe incluir la base mutada en el centro

### Rendimiento
- Para datasets grandes (>10,000 mutaciones), considere:
  - Reducir el número de muestras visualizadas
  - Usar n_signatures menor
  - Procesar por lotes

## Solución de Problemas

### "No mutations found in data"
- Verificar que REF y ALT contengan nucleótidos válidos
- Asegurar que las mutaciones sean SNVs (no indels)
- Verificar el formato del contexto trinucleotídico

### Firmas muy similares
- Puede indicar sobreajuste
- Pruebe con menos firmas (n_signatures menor)
- Verifique que hay suficientes mutaciones por muestra

### Memoria insuficiente
- Reduzca el número de muestras
- Procese en lotes
- Use una máquina con más RAM

## Referencias

- Alexandrov et al. (2013) "Signatures of mutational processes in human cancer"
- COSMIC Mutational Signatures: https://cancer.sanger.ac.uk/signatures/
- Documentación de sklearn NMF: https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.NMF.html 