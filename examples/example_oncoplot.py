"""
Ejemplo de uso del oncoplot en pyMut.

Este script demuestra cómo generar oncoplots usando la librería pyMut,
mostrando diferentes parámetros de configuración y casos de uso.

Un oncoplot (también conocido como waterfall plot en algunos contextos) es una
visualización fundamental en genómica del cáncer que muestra patrones de
mutación a través de muestras y genes en formato heatmap.

Autor: pyMut Development Team
Fecha: 2024
"""

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

# Añadir src al path para poder importar pyMut
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

try:
    from pyMut import PyMutation
    print("✓ PyMutation importado correctamente")
except ImportError as e:
    print(f"✗ Error al importar PyMutation: {e}")
    sys.exit(1)


def main():
    """Función principal que ejecuta todos los ejemplos de oncoplot."""
    print("=" * 60)
    print("EJEMPLO DE ONCOPLOT - pyMut")
    print("=" * 60)
    print()
    
    try:
        # Cargar datos desde el archivo de ejemplo
        data_path = os.path.join('src', 'pyMut', 'data', 'examples', 'tcga_laml_converted.tsv')
        if not os.path.exists(data_path):
            # Ruta alternativa si se ejecuta desde otro directorio
            data_path = os.path.join('..', 'src', 'pyMut', 'data', 'examples', 'tcga_laml_converted.tsv')
        
        data = pd.read_csv(data_path, sep='\t')
        print(f"✓ Datos cargados desde: {data_path}")
        print(f"  Dimensiones: {data.shape}")
        
        # Analizar estructura de datos
        print("\n" + "=" * 40)
        print("ANÁLISIS DE ESTRUCTURA DE DATOS")
        print("=" * 40)
        analyze_data_structure(data)
        
        # Configurar salida de alta calidad para todos los plots
        PyMutation.configure_high_quality_plots()
        print("\n✓ Configuración de alta calidad activada")
        
        # Crear objeto PyMutation
        py_mut = PyMutation(data)
        print("✓ Objeto PyMutation creado correctamente")
        
        # Crear directorio de salidas dentro de examples/
        script_dir = os.path.dirname(os.path.abspath(__file__))
        output_dir = os.path.join(script_dir, 'oncoplot_outputs')
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            print(f"✓ Directorio de salida creado: {output_dir}")
        
        print("\n" + "=" * 40)
        print("GENERANDO ONCOPLOTS")
        print("=" * 40)
        
        # 1. Oncoplot básico
        print("\n1. Generando oncoplot básico...")
        basic_fig = py_mut.oncoplot(
            title="Oncoplot Básico - TCGA LAML",
            show_interactive=False
        )
        basic_output = os.path.join(output_dir, 'oncoplot_basico.png')
        basic_fig.savefig(basic_output)
        print(f"   ✓ Guardado: {basic_output}")
        plt.close(basic_fig)
        
        # 2. Oncoplot personalizado con más genes
        print("\n2. Generando oncoplot personalizado...")
        custom_fig = py_mut.oncoplot(
            title="Oncoplot Personalizado - Top 15 Genes",
            top_genes_count=15,
            max_samples=50,
            figsize=(18, 10),
            show_interactive=False
        )
        custom_output = os.path.join(output_dir, 'oncoplot_personalizado.png')
        custom_fig.savefig(custom_output)
        print(f"   ✓ Guardado: {custom_output}")
        plt.close(custom_fig)
        
        # 3. Oncoplot enfocado (pocos genes, muchas muestras)
        print("\n3. Generando oncoplot enfocado...")
        focused_fig = py_mut.oncoplot(
            title="Oncoplot Enfocado - Top 5 Genes",
            top_genes_count=5,
            max_samples=100,
            figsize=(20, 6),
            show_interactive=False
        )
        focused_output = os.path.join(output_dir, 'oncoplot_enfocado.png')
        focused_fig.savefig(focused_output)
        print(f"   ✓ Guardado: {focused_output}")
        plt.close(focused_fig)
        
        # 4. Demostrar modo interactivo (solo si el script se ejecuta directamente)
        if __name__ == "__main__":
            print("\n4. Ejemplo de modo interactivo...")
            print("   (Se abrirá una ventana interactiva - ciérrala para continuar)")
            interactive_fig = py_mut.oncoplot(
                title="Oncoplot Interactivo",
                top_genes_count=8,
                show_interactive=True
            )
            plt.close(interactive_fig)
            print("   ✓ Modo interactivo completado")
        
        # 5. Múltiples formatos de salida
        print("\n5. Generando múltiples formatos...")
        multi_format_fig = py_mut.oncoplot(
            title="Oncoplot - Múltiples Formatos",
            top_genes_count=10
        )
        
        formats = ['png', 'pdf', 'svg']
        for fmt in formats:
            format_output = os.path.join(output_dir, f'oncoplot_multi_formato.{fmt}')
            if fmt == 'svg':
                multi_format_fig.savefig(format_output, format=fmt)
            else:
                multi_format_fig.savefig(format_output, format=fmt, dpi=300)
            print(f"   ✓ Guardado: {format_output}")
        
        plt.close(multi_format_fig)
        
        # 6. Demostrar validación de parámetros
        print("\n6. Demostrando validación de parámetros...")
        demonstrate_parameter_validation(py_mut)
        
        print("\n" + "=" * 60)
        print("RESUMEN DE ARCHIVOS GENERADOS")
        print("=" * 60)
        
        generated_files = []
        for file in os.listdir(output_dir):
            file_path = os.path.join(output_dir, file)
            size = os.path.getsize(file_path)
            generated_files.append((file, size))
        
        generated_files.sort()
        for filename, size in generated_files:
            size_kb = size / 1024
            print(f"  {filename:<30} ({size_kb:.1f} KB)")
        
        print(f"\n✓ Total de archivos generados: {len(generated_files)}")
        print(f"✓ Ubicación: {os.path.abspath(output_dir)}")
        
        print("\n" + "=" * 60)
        print("EJEMPLO COMPLETADO EXITOSAMENTE")
        print("=" * 60)
        print("\nTodos los oncoplots han sido generados correctamente.")
        print("Los archivos están disponibles en la carpeta examples/oncoplot_outputs/")
        print("\nPuedes revisar los diferentes estilos y configuraciones")
        print("comparando los archivos generados.")
        
    except FileNotFoundError as e:
        print(f"✗ Error: No se pudo encontrar el archivo de datos: {e}")
        print("  Asegúrate de que el archivo tcga_laml_converted.tsv existe en la ubicación esperada.")
        sys.exit(1)
    except Exception as e:
        print(f"✗ Error inesperado: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def analyze_data_structure(data: pd.DataFrame) -> None:
    """
    Analiza y muestra información sobre la estructura de los datos.
    
    Args:
        data: DataFrame con los datos de mutaciones
    """
    print(f"Columnas disponibles ({len(data.columns)}):")
    
    # Categorizar columnas
    sample_columns = [col for col in data.columns if col.startswith('TCGA-')]
    standard_columns = ['Hugo_Symbol', 'Variant_Classification', 'REF', 'ALT']
    other_columns = [col for col in data.columns 
                    if col not in sample_columns and col not in standard_columns]
    
    print(f"  Columnas estándar ({len(standard_columns)}): {', '.join(standard_columns)}")
    print(f"  Columnas de muestras ({len(sample_columns)}): {sample_columns[:5]}{'...' if len(sample_columns) > 5 else ''}")
    if other_columns:
        print(f"  Otras columnas ({len(other_columns)}): {', '.join(other_columns[:3])}{'...' if len(other_columns) > 3 else ''}")
    
    print(f"\nInformación de genes:")
    unique_genes = data['Hugo_Symbol'].nunique()
    print(f"  Genes únicos: {unique_genes}")
    if unique_genes > 0:
        top_genes = data['Hugo_Symbol'].value_counts().head(3)
        print(f"  Top 3 genes: {', '.join(top_genes.index.tolist())}")
    
    print(f"\nInformación de variantes:")
    if 'Variant_Classification' in data.columns:
        unique_variants = data['Variant_Classification'].nunique()
        variant_counts = data['Variant_Classification'].value_counts()
        print(f"  Tipos únicos: {unique_variants}")
        print(f"  Top 3 tipos: {', '.join(variant_counts.head(3).index.tolist())}")
    
    print(f"\nInformación de muestras:")
    print(f"  Total de muestras detectadas: {len(sample_columns)}")
    
    # Estimar número de mutaciones
    mutation_count = 0
    for col in sample_columns[:10]:  # Solo revisar las primeras 10 para eficiencia
        if col in data.columns:
            mutations_in_col = data[col].apply(
                lambda x: '|' in str(x) and str(x) not in ['A|A', 'T|T', 'C|C', 'G|G', '.', '']
            ).sum()
            mutation_count += mutations_in_col
    
    print(f"  Estimación de mutaciones (primeras 10 muestras): {mutation_count}")


def demonstrate_parameter_validation(py_mut: PyMutation) -> None:
    """
    Demuestra la validación de parámetros del oncoplot.
    
    Args:
        py_mut: Objeto PyMutation para testing
    """
    print("   Probando parámetros válidos...")
    
    # Parámetros válidos
    try:
        fig = py_mut.oncoplot(top_genes_count=5, max_samples=10)
        plt.close(fig)
        print("   ✓ Parámetros válidos aceptados")
    except Exception as e:
        print(f"   ✗ Error inesperado con parámetros válidos: {e}")
    
    # Parámetros inválidos
    print("   Probando parámetros inválidos...")
    
    invalid_params = [
        ({"top_genes_count": 0}, "top_genes_count = 0"),
        ({"max_samples": -1}, "max_samples = -1"),
        ({"figsize": (0, 5)}, "figsize inválido"),
    ]
    
    for params, description in invalid_params:
        try:
            py_mut.oncoplot(**params)
            print(f"   ✗ ERROR: {description} debería haber fallado")
        except ValueError:
            print(f"   ✓ {description} correctamente rechazado")
        except Exception as e:
            print(f"   ? {description} falló con error inesperado: {e}")


if __name__ == "__main__":
    main() 