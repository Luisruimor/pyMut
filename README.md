# pyMut

Una librería Python para visualizar mutaciones genéticas a partir de archivos .TSV.

## Descripción

pyMut es una herramienta de visualización para datos de mutaciones genéticas, inspirada en herramientas como Maftools y Mutscape. Permite generar visualizaciones de resumen para entender las mutaciones genéticas en un dataset de formato TSV.

## Características

- Visualizaciones de resumen de mutaciones:
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
fig = pyMut.summary_plot(title="Resumen de mutaciones")
fig.savefig("summary.png")
```

## Requisitos

- Python 3.7+
- pandas
- matplotlib
- numpy
- seaborn

## Documentación

Para obtener más información, consulte la [documentación completa](https://pymut.readthedocs.io/).

## Contribuir

Las contribuciones son bienvenidas. Por favor, abra un issue primero para discutir lo que le gustaría cambiar.

## Licencia

[MIT](LICENSE)