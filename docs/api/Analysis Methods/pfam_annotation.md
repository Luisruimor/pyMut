# Anotación de Dominios Pfam - annotate_pfam y pfam_domains

Los métodos **annotate_pfam** y **pfam_domains** permiten anotar mutaciones con información de dominios proteicos Pfam y analizar la distribución de mutaciones en estos dominios.

## ¿Qué es la Anotación Pfam?

Pfam es una base de datos de familias de proteínas que contiene alineamientos múltiples y modelos ocultos de Markov (HMMs) para dominios proteicos. La anotación Pfam permite identificar qué dominios funcionales están afectados por las mutaciones.

## Características Principales

- **Anotación automática**: Mapea mutaciones a dominios Pfam usando bases de datos internas
- **Múltiples fuentes**: Utiliza tanto bases de datos locales como anotaciones VEP
- **Extracción automática**: Identifica automáticamente IDs de UniProt y posiciones de aminoácidos
- **Análisis estadístico**: Genera resúmenes y estadísticas de dominios afectados
- **Visualización opcional**: Puede generar gráficos de dominios más afectados
- **Filtrado flexible**: Permite incluir o excluir mutaciones sinónimas

## Método annotate_pfam - Anotación de Dominios

### Uso Básico

```python
from pyMut.input import read_maf

# Cargar datos con anotaciones VEP
py_mut = read_maf("mutations_vep_annotated.maf")

# Anotar con dominios Pfam
annotated = py_mut.annotate_pfam()

print(f"Mutaciones anotadas: {len(annotated.data)}")
print(f"Columnas Pfam añadidas: {[col for col in annotated.data.columns if 'Pfam' in col]}")
```

### Parámetros de annotate_pfam

#### db_conn (duckdb.DuckDBPyConnection, opcional)
- **Descripción**: Conexión a base de datos DuckDB con datos Pfam
- **Default**: Se crea automáticamente una conexión interna
- **Uso**: Para usar bases de datos personalizadas

#### aa_column (str, default='aa_pos')
- **Descripción**: Nombre de la columna que contiene información de posición de aminoácido
- **Formato esperado**: Debe contener UniProt ID y posición (ej: "P53_R273H")
- **Ejemplo**: `"HGVSp"`, `"aa_change"`, `"protein_change"`

#### auto_extract (bool, default=True)
- **Descripción**: Si extraer automáticamente UniProt ID y posición de otras columnas
- **True**: Busca en columnas VEP como SWISSPROT, HGVSp
- **False**: Usa solo la columna especificada en aa_column

#### prefer_database (bool, default=True)
- **Descripción**: Preferencia entre base de datos local y anotaciones VEP
- **True**: Prioriza base de datos DuckDB interna
- **False**: Prioriza anotaciones VEP existentes

### Ejemplos de annotate_pfam

```python
from pyMut.input import read_maf
import logging

# Configurar logging para ver detalles
logging.basicConfig(level=logging.INFO)

# Cargar datos con anotaciones VEP
py_mut = read_maf("src/pyMut/data/examples/tcga_laml_VEP_annotated.maf.gz")
print(f"Mutaciones originales: {len(py_mut.data)}")

# Verificar columnas disponibles para anotación
print("Columnas disponibles:")
for col in py_mut.data.columns:
    if any(term in col.lower() for term in ['uniprot', 'swissprot', 'hgvs', 'protein']):
        print(f"  - {col}")

# Anotación básica con configuración automática
annotated = py_mut.annotate_pfam()

# Verificar resultados
pfam_columns = [col for col in annotated.data.columns if 'Pfam' in col]
print(f"\nColumnas Pfam añadidas: {pfam_columns}")

# Contar mutaciones anotadas
annotated_count = annotated.data['PfamDomain'].notna().sum()
print(f"Mutaciones con anotación Pfam: {annotated_count}")

# Mostrar ejemplos de anotaciones
if annotated_count > 0:
    examples = annotated.data[annotated.data['PfamDomain'].notna()].head(3)
    print("\nEjemplos de anotaciones:")
    for _, row in examples.iterrows():
        print(f"  Gen: {row.get('Hugo_Symbol', 'N/A')}")
        print(f"  Dominio: {row.get('PfamDomain', 'N/A')}")
        print(f"  Descripción: {row.get('PfamDescription', 'N/A')}")
        print("  ---")
```

## Método pfam_domains - Análisis de Dominios

### Uso Básico

```python
# Analizar dominios Pfam más afectados
domain_analysis = py_mut.pfam_domains()

print("Dominios más afectados:")
print(domain_analysis)
```

### Parámetros de pfam_domains

#### aa_column (str, default='aa_pos')
- **Descripción**: Columna con información de posición de aminoácido
- **Mismo uso que en annotate_pfam**

#### summarize_by (str, default='PfamDomain')
- **Descripción**: Campo por el cual agrupar el análisis
- **Opciones**: `'PfamDomain'`, `'PfamFamily'`, `'PfamClan'`
- **Ejemplo**: `'PfamFamily'` para análisis por familia de proteínas

#### top_n (int, default=10)
- **Descripción**: Número de dominios principales a mostrar
- **Rango**: 1-100 (recomendado: 5-20)

#### include_synonymous (bool, default=False)
- **Descripción**: Si incluir mutaciones sinónimas en el análisis
- **False**: Solo mutaciones con impacto funcional
- **True**: Todas las mutaciones

#### plot (bool, default=False)
- **Descripción**: Si generar gráfico de barras de los resultados
- **Requiere**: matplotlib instalado

### Ejemplos Completos

```python
from pyMut.input import read_maf
import matplotlib.pyplot as plt

# Cargar y anotar datos
py_mut = read_maf("mutations_vep.maf")
annotated = py_mut.annotate_pfam()

print("=== Análisis de Dominios Pfam ===")

# Análisis básico de dominios
print("\n1. Dominios más afectados:")
domain_summary = annotated.pfam_domains(
    summarize_by='PfamDomain',
    top_n=10,
    include_synonymous=False
)
print(domain_summary)

# Análisis por familia de proteínas
print("\n2. Familias de proteínas más afectadas:")
family_summary = annotated.pfam_domains(
    summarize_by='PfamFamily',
    top_n=8,
    include_synonymous=False
)
print(family_summary)

# Análisis incluyendo mutaciones sinónimas
print("\n3. Análisis incluyendo sinónimas:")
all_mutations = annotated.pfam_domains(
    summarize_by='PfamDomain',
    top_n=10,
    include_synonymous=True
)
print(all_mutations)

# Generar gráfico
print("\n4. Generando visualización...")
domain_plot = annotated.pfam_domains(
    summarize_by='PfamDomain',
    top_n=15,
    include_synonymous=False,
    plot=True
)

# Guardar gráfico
plt.savefig("pfam_domains_analysis.png", dpi=300, bbox_inches='tight')
print("Gráfico guardado como 'pfam_domains_analysis.png'")
```

## Análisis Avanzado de Dominios

### Análisis por Gen Específico

```python
def analizar_dominios_gen(py_mut, gene_name):
    """
    Analiza dominios Pfam afectados en un gen específico
    """
    # Filtrar por gen
    gene_data = py_mut.data[py_mut.data['Hugo_Symbol'] == gene_name].copy()
    
    if len(gene_data) == 0:
        print(f"No se encontraron mutaciones en {gene_name}")
        return
    
    print(f"=== Análisis de Dominios: {gene_name} ===")
    print(f"Mutaciones totales: {len(gene_data)}")
    
    # Crear objeto temporal para análisis
    from pyMut.core import PyMutation
    gene_pymut = PyMutation(gene_data, py_mut.metadata, py_mut.samples)
    
    # Anotar dominios
    annotated_gene = gene_pymut.annotate_pfam()
    
    # Analizar dominios
    if len(annotated_gene.data) > 0:
        domain_analysis = annotated_gene.pfam_domains(
            top_n=20,
            include_synonymous=False
        )
        
        print("\nDominios afectados:")
        print(domain_analysis)
        
        # Estadísticas adicionales
        total_with_domains = annotated_gene.data['PfamDomain'].notna().sum()
        print(f"\nMutaciones con dominios anotados: {total_with_domains}")
        
        if total_with_domains > 0:
            unique_domains = annotated_gene.data['PfamDomain'].nunique()
            print(f"Dominios únicos afectados: {unique_domains}")

# Ejemplo de uso
analizar_dominios_gen(py_mut, "TP53")
analizar_dominios_gen(py_mut, "KRAS")
```

### Comparación de Dominios entre Grupos

```python
def comparar_dominios_grupos(py_mut, grupo1_samples, grupo2_samples):
    """
    Compara dominios afectados entre dos grupos de muestras
    """
    # Filtrar grupos
    grupo1 = py_mut.filter_by_chrom_sample(sample=grupo1_samples)
    grupo2 = py_mut.filter_by_chrom_sample(sample=grupo2_samples)
    
    # Anotar dominios
    g1_annotated = grupo1.annotate_pfam()
    g2_annotated = grupo2.annotate_pfam()
    
    print("=== Comparación de Dominios entre Grupos ===")
    
    # Análisis grupo 1
    print(f"\nGrupo 1 ({len(grupo1_samples)} muestras):")
    g1_domains = g1_annotated.pfam_domains(top_n=10)
    print(g1_domains)
    
    # Análisis grupo 2
    print(f"\nGrupo 2 ({len(grupo2_samples)} muestras):")
    g2_domains = g2_annotated.pfam_domains(top_n=10)
    print(g2_domains)
    
    # Dominios únicos por grupo
    g1_unique = set(g1_annotated.data['PfamDomain'].dropna().unique())
    g2_unique = set(g2_annotated.data['PfamDomain'].dropna().unique())
    
    print(f"\nDominios únicos Grupo 1: {len(g1_unique - g2_unique)}")
    print(f"Dominios únicos Grupo 2: {len(g2_unique - g1_unique)}")
    print(f"Dominios compartidos: {len(g1_unique & g2_unique)}")

# Ejemplo de uso
grupo_a = ["TCGA-AB-2802", "TCGA-AB-2803", "TCGA-AB-2869"]
grupo_b = ["TCGA-AB-2988", "TCGA-AB-2989", "TCGA-AB-2990"]

comparar_dominios_grupos(py_mut, grupo_a, grupo_b)
```

## Tipos de Anotaciones Pfam

### Información Disponible

```python
# Después de la anotación, las siguientes columnas están disponibles:
annotated = py_mut.annotate_pfam()

pfam_info = {
    'PfamDomain': 'Identificador del dominio Pfam (ej: PF00046)',
    'PfamDescription': 'Descripción del dominio (ej: Homeobox domain)',
    'PfamFamily': 'Familia de proteínas',
    'PfamClan': 'Clan de dominios relacionados',
    'PfamStart': 'Posición de inicio del dominio',
    'PfamEnd': 'Posición final del dominio',
    'UniProtID': 'Identificador UniProt de la proteína',
    'AAPosition': 'Posición del aminoácido mutado'
}

print("Información Pfam disponible:")
for column, description in pfam_info.items():
    if column in annotated.data.columns:
        print(f"✅ {column}: {description}")
    else:
        print(f"❌ {column}: {description}")
```

### Ejemplos de Dominios Comunes

```python
# Dominios frecuentemente mutados en cáncer
cancer_domains = {
    'PF00046': 'Homeobox domain',
    'PF00096': 'Zinc finger, C2H2 type',
    'PF00169': 'PH domain',
    'PF00018': 'SH3 domain',
    'PF00017': 'SH2 domain',
    'PF00069': 'Protein kinase domain',
    'PF00481': 'Protein phosphatase 2C',
    'PF00125': 'Core histone H2A/H2B/H3/H4'
}

# Verificar presencia en datos
annotated = py_mut.annotate_pfam()
found_domains = annotated.data['PfamDomain'].value_counts()

print("Dominios de cáncer encontrados:")
for domain_id, description in cancer_domains.items():
    if domain_id in found_domains.index:
        count = found_domains[domain_id]
        print(f"✅ {domain_id} ({description}): {count} mutaciones")
    else:
        print(f"❌ {domain_id} ({description}): No encontrado")
```

## Casos de Uso Específicos

### Análisis de Genes Supresores de Tumores

```python
# Genes supresores de tumores conocidos
tumor_suppressors = ["TP53", "RB1", "APC", "BRCA1", "BRCA2", "VHL"]

print("=== Análisis de Dominios en Genes Supresores ===")

for gene in tumor_suppressors:
    gene_data = py_mut.data[py_mut.data['Hugo_Symbol'] == gene]
    
    if len(gene_data) > 0:
        # Crear objeto temporal
        gene_pymut = PyMutation(gene_data, py_mut.metadata, py_mut.samples)
        annotated_gene = gene_pymut.annotate_pfam()
        
        print(f"\n{gene}:")
        print(f"  Mutaciones: {len(gene_data)}")
        
        # Dominios afectados
        domains = annotated_gene.data['PfamDomain'].value_counts()
        if len(domains) > 0:
            print(f"  Dominios únicos: {len(domains)}")
            print(f"  Dominio más afectado: {domains.index[0]} ({domains.iloc[0]} mut)")
        else:
            print("  Sin anotaciones de dominios")
```

### Análisis de Hotspots en Dominios

```python
def analizar_hotspots_dominios(py_mut, min_mutations=5):
    """
    Identifica hotspots de mutaciones en dominios específicos
    """
    # Anotar dominios
    annotated = py_mut.annotate_pfam()
    
    # Filtrar mutaciones con dominios anotados
    with_domains = annotated.data[annotated.data['PfamDomain'].notna()]
    
    if len(with_domains) == 0:
        print("No hay mutaciones con dominios anotados")
        return
    
    print("=== Análisis de Hotspots en Dominios ===")
    
    # Agrupar por dominio y posición
    hotspots = with_domains.groupby(['PfamDomain', 'AAPosition']).size()
    hotspots = hotspots[hotspots >= min_mutations].sort_values(ascending=False)
    
    print(f"Hotspots encontrados (≥{min_mutations} mutaciones):")
    
    for (domain, position), count in hotspots.head(10).items():
        # Obtener información adicional
        domain_info = with_domains[
            (with_domains['PfamDomain'] == domain) & 
            (with_domains['AAPosition'] == position)
        ].iloc[0]
        
        gene = domain_info.get('Hugo_Symbol', 'Unknown')
        description = domain_info.get('PfamDescription', 'Unknown')
        
        print(f"  {domain} (pos {position}): {count} mutaciones")
        print(f"    Gen: {gene}")
        print(f"    Descripción: {description}")
        print("    ---")

# Ejemplo de uso
analizar_hotspots_dominios(py_mut, min_mutations=3)
```

## Integración con Otros Análisis

### Combinación con Análisis TMB

```python
# Analizar dominios en muestras con alta carga mutacional
tmb_results = py_mut.calculate_tmb_analysis()
high_tmb_samples = tmb_results['analysis'][
    tmb_results['analysis']['TMB_Total_Normalized'] > 2.0
]['Sample'].tolist()

# Filtrar muestras de alta TMB
high_tmb_data = py_mut.filter_by_chrom_sample(sample=high_tmb_samples)

# Analizar dominios en estas muestras
high_tmb_annotated = high_tmb_data.annotate_pfam()
domain_analysis = high_tmb_annotated.pfam_domains(
    top_n=15,
    include_synonymous=False
)

print("Dominios más afectados en muestras de alta TMB:")
print(domain_analysis)
```

### Preparación para Análisis Funcional

```python
# Exportar datos anotados para análisis downstream
annotated = py_mut.annotate_pfam()

# Crear tabla de resumen para análisis funcional
functional_summary = annotated.data[
    annotated.data['PfamDomain'].notna()
].groupby(['Hugo_Symbol', 'PfamDomain', 'PfamDescription']).agg({
    'Tumor_Sample_Barcode': 'count',
    'Variant_Classification': lambda x: x.value_counts().index[0]
}).rename(columns={
    'Tumor_Sample_Barcode': 'MutationCount',
    'Variant_Classification': 'MostCommonType'
}).reset_index()

# Guardar para análisis externo
functional_summary.to_csv("pfam_functional_summary.csv", index=False)
print("Resumen funcional guardado en 'pfam_functional_summary.csv'")
```

## Solución de Problemas

### Datos sin Anotaciones VEP

```python
# Si los datos no tienen anotaciones VEP previas
try:
    annotated = py_mut.annotate_pfam()
    print("✅ Anotación exitosa")
except Exception as e:
    print(f"❌ Error en anotación: {e}")
    print("Posibles soluciones:")
    print("1. Verificar que los datos tengan columnas UniProt/SWISSPROT")
    print("2. Usar datos pre-anotados con VEP")
    print("3. Especificar columna correcta con aa_column")
```

### Optimización para Datasets Grandes

```python
# Para datasets muy grandes, procesar por lotes
def anotar_por_lotes(py_mut, batch_size=10000):
    """
    Anota dominios Pfam procesando datos por lotes
    """
    total_rows = len(py_mut.data)
    annotated_batches = []
    
    for i in range(0, total_rows, batch_size):
        print(f"Procesando lote {i//batch_size + 1}...")
        
        # Crear lote
        batch_data = py_mut.data.iloc[i:i+batch_size].copy()
        batch_pymut = PyMutation(batch_data, py_mut.metadata, py_mut.samples)
        
        # Anotar lote
        batch_annotated = batch_pymut.annotate_pfam()
        annotated_batches.append(batch_annotated.data)
    
    # Combinar resultados
    final_data = pd.concat(annotated_batches, ignore_index=True)
    return PyMutation(final_data, py_mut.metadata, py_mut.samples)

# Uso para datasets grandes
if len(py_mut.data) > 50000:
    annotated = anotar_por_lotes(py_mut, batch_size=5000)
else:
    annotated = py_mut.annotate_pfam()
```