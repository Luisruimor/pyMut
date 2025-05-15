# pyMut

Una librería Python para visualizar mutaciones genéticas a partir de archivos .TSV.

## Descripción

pyMut es una herramienta de visualización para datos de mutaciones genéticas, inspirada en herramientas como Maftools y Mutscape. Permite generar visualizaciones de resumen para entender las mutaciones genéticas en un dataset de formato TSV.

## Características

- Visualizaciones de resumen de mutaciones:
  - **Variant Classification**: Distribución de clasificaciones de variantes
  - **Variant Type**: Distribución de tipos de variantes (SNP, INS, DEL, etc.)
  - **SNV Class**: Distribución de clases de SNV (cambios nucleotídicos como A>G, C>T, etc.)
  - **Variants per Sample**: Distribución de variantes por muestra y mediana (TMB)
  - **Variant Classification Summary**: Diagrama de cajas y bigotes que muestra la distribución de cada tipo de variante entre las muestras, permitiendo identificar qué clasificaciones presentan mayor variabilidad entre pacientes

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

# Generar un gráfico de resumen completo
fig = pyMut.summary_plot(title="Resumen de mutaciones")
fig.savefig("summary.png")

# Generar visualizaciones individuales
tmb_fig = pyMut.variants_per_sample_plot(title="Carga Mutacional por Muestra")
tmb_fig.savefig("variants_per_sample.png")

boxplot_fig = pyMut.variant_classification_summary_plot(title="Distribución de Variantes por Muestra")
boxplot_fig.savefig("variant_classification_summary.png")
```

## Formatos de Datos Soportados

pyMut puede trabajar con datos en dos formatos principales:

- **Formato Largo**: Cada fila representa una mutación, con columnas como `Variant_Classification` y `Tumor_Sample_Barcode`.
- **Formato Ancho**: Las muestras se representan como columnas (por ejemplo, columnas con nombres como `TCGA-XX-YYYY`).

La biblioteca detecta automáticamente el formato y adapta las visualizaciones según corresponda. Para el formato ancho, pyMut puede interpretar datos en diferentes notaciones:
- Genotipos separados por pipe (`|`): como "A|B" donde A y B son diferentes alelos
- Genotipos separados por barra (`/`): como "A/B"
- Otros formatos: valores numéricos o textuales que indican presencia de variante

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