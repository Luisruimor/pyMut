# pyMut

Una librería Python para visualizar mutaciones genéticas a partir de archivos .TSV.

## Descripción

pyMut es una herramienta de visualización para datos de mutaciones genéticas. Permite generar visualizaciones de resumen para entender mejor los patrones y distribuciones de las mutaciones en los datos.

## Características

- Visualizaciones de resumen estadístico:
  - Variant Classification: Distribución de clasificaciones de variantes
  - Variant Type: Distribución de tipos de variantes (SNP, INS, DEL, etc.)
  - SNV Class: Distribución de clases de SNV (cambios nucleotídicos como A>G, C>T, etc.)

## Instalación

```bash
pip install pyMut
```

## Uso básico

```python
from pyMut import PyMutation
import pandas as pd

# Cargar datos de mutaciones desde un archivo TSV
data = pd.read_csv("mutations.tsv", sep="\t")

# Crear un objeto PyMutation
pyMut = PyMutation(data)

# Generar un gráfico de resumen
fig = pyMut.summary_plot()
fig.savefig("summary.png")
```

## Proyección futura

En el futuro, pyMut se ampliará con más visualizaciones y funcionalidades para un análisis más completo de mutaciones genéticas, siempre manteniendo la simplicidad y facilidad de uso.
