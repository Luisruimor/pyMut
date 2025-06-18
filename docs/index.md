# pyMut - Librería de Visualización de Mutaciones Genómicas

**pyMut** es una librería Python especializada en la visualización de datos de mutaciones genómicas, diseñada siguiendo principios de programación orientada a objetos (OOP) y mejores prácticas de desarrollo.

## Características Principales

- 🧬 **Visualizaciones especializadas** para datos genómicos
- 📊 **Plot de resumen** con múltiples subgráficos informativos
- 🎨 **Oncoplot** para visualizar paisajes mutacionales  
- 🏗️ **Arquitectura OOP** limpia y extensible
- 📈 **Alta calidad** automática (DPI 300, formatos vectoriales)
- 🔧 **API intuitiva** y bien documentada
- ✅ **Test suite completo** con >95% cobertura
- 📚 **Documentación exhaustiva** con ejemplos

## Instalación

```bash
# Clonar el repositorio
git clone https://github.com/usuario/pyMut.git
cd pyMut

# Instalar dependencias
pip install -r requirements.txt

# Instalación en modo desarrollo
pip install -e .
```

## Uso Rápido

```python
from pyMut import PyMutation
import pandas as pd

# Configurar alta calidad (recomendado)
PyMutation.configure_high_quality_plots()

# Cargar datos
data = pd.read_csv("mutations.tsv", sep='\t')

# Crear objeto PyMutation  
py_mut = PyMutation(data)

# Generar visualizaciones
summary_fig = py_mut.summary_plot()
oncoplot_fig = py_mut.oncoplot()

# Guardar (automáticamente alta calidad)
summary_fig.savefig("summary.png")
oncoplot_fig.savefig("oncoplot.png")
```

## Visualizaciones Disponibles

### 1. Plot de Resumen (`summary_plot`)

Un panel completo con múltiples visualizaciones:
- Clasificación de variantes
- Tipos de variantes  
- Clases de SNV
- Variantes por muestra (TMB)
- Resumen de clasificaciones
- Genes más mutados

### 2. Oncoplot (`oncoplot`)

Visualización heatmap de paisajes mutacionales:
- Genes ordenados por frecuencia de mutación
- Muestras ordenadas por carga mutacional
- Colores por tipo de mutación
- Detección automática de Multi_Hit
- Soporte para formatos TCGA y .GT

### 3. Visualizaciones Individuales

Cada componente del summary plot está disponible individualmente:
- `variant_classification_plot()`
- `variant_type_plot()`
- `snv_class_plot()`
- `variants_per_sample_plot()`
- `variant_classification_summary_plot()`
- `top_mutated_genes_plot()`

## Estructura del Proyecto

```
pyMut/
├── src/pyMut/
│   ├── core.py                    # API principal (PyMutation class)
│   │   ├── summary.py             # Visualizaciones de resumen
│   │   └── oncoplot.py            # Oncoplot especializado
│   ├── utils/
│   │   ├── data_processing.py     # Procesamiento de datos
│   │   └── constants.py           # Constantes del proyecto
│   └── data/examples/             # Datos de ejemplo
├── examples/
│   ├── example_summary.py         # Ejemplo completo de summary plots
│   └── example_oncoplot.py        # Ejemplo completo de oncoplot
├── tests/
│   ├── test_core.py              # Tests de funcionalidad principal
│   └── test_oncoplot.py          # Tests específicos de oncoplot
├── docs/
│   ├── index.md                  # Esta documentación
│   ├── summary-plot.md           # Guía detallada de summary plots
│   └── oncoplot.md               # Guía detallada de oncoplot
└── .cursor/rules/                # Reglas para asistentes AI
```

## Casos de Uso

### Análisis de Datasets TCGA

```python
# Cargar datos TCGA
tcga_data = pd.read_csv("TCGA_mutations.maf", sep='\t')

# Análisis completo
py_mut = PyMutation(tcga_data)

# Panel de resumen
summary = py_mut.summary_plot(
    title="TCGA-LAML Mutation Analysis",
    max_samples=100,
    top_genes_count=20
)

# Oncoplot enfocado  
oncoplot = py_mut.oncoplot(
    title="Top Mutated Genes",
    top_genes_count=15,
    max_samples=200
)
```

### Análisis de Genes Específicos

```python
# Filtrar por genes de interés
oncogenes = ['TP53', 'KRAS', 'PIK3CA', 'EGFR', 'BRAF']
filtered_data = data[data['Hugo_Symbol'].isin(oncogenes)]

py_mut = PyMutation(filtered_data)
focused_plot = py_mut.oncoplot(
    title="Oncogenes Principales",
    figsize=(12, 6)
)
```

### Exportación Multi-formato

```python
fig = py_mut.summary_plot()

# Múltiples formatos de alta calidad
for fmt in ['png', 'pdf', 'svg']:
    fig.savefig(f'summary.{fmt}', dpi=300 if fmt != 'svg' else None)
```

## Formatos de Datos Soportados

### Estructura Requerida

```python
# Columnas mínimas requeridas
columns = [
    'Hugo_Symbol',           # Símbolo del gen
    'Variant_Classification', # Tipo de mutación
    'REF',                   # Alelo de referencia
    'ALT',                   # Alelo alternativo
    'TCGA-XX-XXXX'          # Columnas de muestras (auto-detectadas)
]
```

### Formatos de Genotipo

```python
# Formatos soportados:
"A|G"    # Pipe-separated (estándar)
"A/G"    # Slash-separated
"A:G"    # Colon-separated
"A;G"    # Semicolon-separated

# Casos especiales:
"A|A"    # Sin mutación (homocigoto referencia)
"G|G"    # Mutación homocigota
"."      # Datos faltantes
""       # Datos vacíos
```

## API Reference

### Clase Principal: `PyMutation`

```python
class PyMutation:
    def __init__(self, data: pd.DataFrame)
    
    @staticmethod
    def configure_high_quality_plots()
    
    def summary_plot(self, figsize=(16, 12), title="Mutation Summary", 
                     max_samples=None, top_genes_count=10, 
                     show_interactive=False) -> plt.Figure
    
    def oncoplot(self, figsize=(16, 10), title="Oncoplot",
                 top_genes_count=10, max_samples=180,
                 show_interactive=False) -> plt.Figure
    
    def variant_classification_plot(self, **kwargs) -> plt.Figure
    def variant_type_plot(self, **kwargs) -> plt.Figure  
    def snv_class_plot(self, **kwargs) -> plt.Figure
    def variants_per_sample_plot(self, **kwargs) -> plt.Figure
    def variant_classification_summary_plot(self, **kwargs) -> plt.Figure
    def top_mutated_genes_plot(self, **kwargs) -> plt.Figure
```

### Parámetros Comunes

| Parámetro | Descripción | Valor por Defecto |
|-----------|-------------|-------------------|
| `figsize` | Tamaño de figura (ancho, alto) | Varía por plot |
| `title` | Título del gráfico | Varía por plot |
| `show_interactive` | Mostrar ventana interactiva | `False` |
| `max_samples` | Límite de muestras a mostrar | `None` (todas) |
| `top_genes_count` | Número de genes principales | `10` |

## Guías Detalladas

- **[Guía de Summary Plot](summary-plot.md)**: Documentación completa de las visualizaciones de resumen
- **[Guía de Oncoplot](oncoplot.md)**: Documentación detallada del oncoplot, interpretación y casos de uso

## Ejemplos Completos

- **[Ejemplo de Summary Plots](../examples/example_summary.py)**: Script completo con todos los tipos de summary plots
- **[Ejemplo de Oncoplot](../examples/example_oncoplot.py)**: Script completo con múltiples configuraciones de oncoplot

## Desarrollo y Testing

### Ejecutar Tests

```bash
# Tests completos
python -m pytest tests/ -v

# Tests específicos de oncoplot
python -m pytest tests/test_oncoplot.py -v

# Tests con cobertura
python -m pytest tests/ --cov=src/pyMut --cov-report=html
```

### Ejecutar Ejemplos

```bash
# Ejemplo de summary plots
cd examples && python example_summary.py

# Ejemplo de oncoplot  
cd examples && python example_oncoplot.py
```

## Contribuir

1. Fork el proyecto
2. Crear rama de feature (`git checkout -b feature/nueva-caracteristica`)
3. Hacer commit de cambios (`git commit -am 'Agregar nueva característica'`)
4. Push a la rama (`git push origin feature/nueva-caracteristica`)
5. Crear Pull Request

### Estándares de Código

- **PEP 8**: Estilo de código Python
- **Type hints**: Obligatorios para todas las funciones públicas
- **Docstrings**: Documentación completa en estilo Google
- **Tests**: Cobertura mínima del 90%

## Licencia

Este proyecto está bajo la licencia MIT. Ver `LICENSE` para más detalles.

## Contacto

- **Repositorio**: https://github.com/usuario/pyMut
- **Documentación**: https://pymut.readthedocs.io
- **Issues**: https://github.com/usuario/pyMut/issues

---

Para comenzar, revisa los [ejemplos completos](../examples/) o las [guías detalladas](summary-plot.md) de cada tipo de visualización.
