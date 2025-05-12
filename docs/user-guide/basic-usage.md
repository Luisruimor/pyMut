# Uso Básico

## Importar la biblioteca

Para empezar a utilizar pyMut, primero debes importar la clase principal, `PyMutation`:

```python
from pyMut import PyMutation
import pandas as pd
```

## Cargar datos

pyMut trabaja con datos de mutaciones en formato DataFrame de pandas. Puedes cargar tus datos desde un archivo TSV:

```python
# Cargar datos desde un archivo TSV
data = pd.read_csv("path/to/mutations.tsv", sep="\t")

# Crear un objeto PyMutation
pyMut = PyMutation(data)
```

## Generar visualizaciones

### Gráfico de Resumen

Un gráfico de resumen muestra estadísticas generales de las mutaciones. Esta visualización incluye:

- Variant Classification: Distribución de clasificaciones de variantes
- Variant Type: Distribución de tipos de variantes (SNP, INS, DEL, etc.)
- SNV Class: Distribución de clases de SNV (cambios nucleotídicos como A>G, C>T, etc.)

Para generar un gráfico de resumen:

```python
# Generar un gráfico de resumen
fig = pyMut.summary_plot(title="Resumen de Mutaciones")

# Guardar el gráfico como un archivo PNG
fig.savefig("summary.png")
```

### Personalización

Puedes personalizar el gráfico de resumen cambiando los parámetros:

```python
# Generar un gráfico de resumen personalizado
fig = pyMut.summary_plot(
    figsize=(15, 12),     # Tamaño de la figura
    title="Análisis de Mutaciones en Pacientes con Cáncer"  # Título personalizado
)

# Guardar el gráfico como un archivo PNG de alta resolución
fig.savefig("summary_personalizado.png", dpi=300)
```

### Visualización interactiva

También puedes mostrar el gráfico de forma interactiva:

```python
# Generar y mostrar el gráfico interactivamente
fig = pyMut.summary_plot()
import matplotlib.pyplot as plt
plt.show()
``` 