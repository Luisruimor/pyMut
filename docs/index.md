# pyMut Documentation

pyMut es una librería Python para visualizar mutaciones genéticas, inspirada en herramientas como Maftools y Mutscape.

## Instalación

```bash
pip install pyMut
```

## Uso Rápido

```python
from pyMut import PyMutation
import pandas as pd

# Cargar datos
data = pd.read_csv("mutations.tsv", sep='\t')

# Crear visualización
py_mut = PyMutation(data)
PyMutation.configure_high_quality_plots()

# Generar summary plot
fig = py_mut.summary_plot(title="Mi Análisis")
fig.savefig("summary.png")
```

## Documentación Disponible

- **[Summary Plot](summary-plot.md)** - Visualización principal que combina múltiples análisis

## Formato de Datos

pyMut requiere datos con estas columnas mínimas:

- `Hugo_Symbol` - Símbolo del gen
- `Variant_Classification` - Tipo de mutación (Missense_Mutation, Silent, etc.)

Opcionalmente:
- `Tumor_Sample_Barcode` - Identificador de muestra
- `REF` / `ALT` - Alelos de referencia y alternativo

## Ejemplo con Datos de Prueba

```python
# Usar datos de ejemplo incluidos
data = pd.read_csv("src/pyMut/data/examples/tcga_laml_converted.tsv", sep='\t')
py_mut = PyMutation(data)

# Generar todas las visualizaciones disponibles
summary_fig = py_mut.summary_plot()
vc_fig = py_mut.variant_classification_plot()
tmb_fig = py_mut.variants_per_sample_plot()
genes_fig = py_mut.top_mutated_genes_plot(mode="variants")
```

## Soporte

- **GitHub**: [https://github.com/your-username/pyMut](https://github.com/your-username/pyMut)
- **Issues**: Reporta problemas o solicita features
- **Tests**: Ejecuta `./run_clean_tests.sh` para verificar instalación
