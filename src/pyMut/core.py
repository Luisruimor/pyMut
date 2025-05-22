import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Dict, Union, Optional, Tuple
from .utils.constants import (
    VARIANT_CLASSIFICATION_COLUMN, VARIANT_TYPE_COLUMN, SAMPLE_COLUMN, 
    GENE_COLUMN, REF_COLUMN, ALT_COLUMN, FUNCOTATION_COLUMN, 
    DEFAULT_SUMMARY_FIGSIZE, DEFAULT_PLOT_FIGSIZE, DEFAULT_PLOT_TITLE,
    DEFAULT_TOP_GENES_COUNT, MODE_VARIANTS, VALID_PLOT_MODES
)

class PyMutation:
    """
    Clase principal para visualizar mutaciones genéticas a partir de datos en formato TSV.
    
    Esta clase proporciona métodos para generar visualizaciones de resumen de datos
    de mutaciones genéticas, mostrando estadísticas generales como distribuciones
    de tipos de variantes, clasificaciones y cambios de nucleótidos.
    
    Attributes:
        data (pd.DataFrame): DataFrame que contiene los datos de mutaciones.
    """
    
    def __init__(self, data: pd.DataFrame):
        """
        Inicializa un objeto PyMutation con un DataFrame de pandas.
        
        Args:
            data (pd.DataFrame): DataFrame conteniendo datos de mutaciones.
        
        Raises:
            ValueError: Si el DataFrame está vacío o no es un DataFrame válido.
        """
        if not isinstance(data, pd.DataFrame):
            raise ValueError("The 'data' parameter must be a pandas DataFrame.")
        
        if data.empty:
            raise ValueError("The provided DataFrame is empty. No data to analyze.")
        
        self.data = data
    
    def save_figure(self, figure: plt.Figure, filename: str, 
                   dpi: int = 300, bbox_inches: str = 'tight', **kwargs) -> None:
        """
        Guarda una figura con configuración de alta calidad por defecto.
        
        Este método centraliza el guardado de figuras para asegurar que todas
        las visualizaciones se guarden con la mejor calidad posible.
        
        Args:
            figure: La figura de matplotlib a guardar.
            filename: Nombre del archivo donde guardar la figura.
            dpi: Resolución en puntos por pulgada (300 = alta calidad).
            bbox_inches: Ajuste de márgenes ('tight' = sin espacios innecesarios).
            **kwargs: Parámetros adicionales para matplotlib.savefig().
        
        Examples:
            >>> py_mut = PyMutation(data)
            >>> fig = py_mut.summary_plot()
            >>> py_mut.save_figure(fig, 'mi_resumen.png')  # Alta calidad automática
            >>> py_mut.save_figure(fig, 'mi_resumen.pdf', dpi=600)  # Muy alta calidad
        """
        figure.savefig(filename, dpi=dpi, bbox_inches=bbox_inches, **kwargs)
        print(f"📁 Figura guardada: {filename} (DPI: {dpi}, márgenes: {bbox_inches})")
    
    @staticmethod
    def configure_high_quality_plots():
        """
        Configura matplotlib para generar gráficos de alta calidad por defecto.
        
        Esta función modifica la configuración global de matplotlib para que
        TODAS las figuras se guarden automáticamente con alta calidad, sin
        necesidad de especificar parámetros cada vez.
        
        Configuraciones aplicadas:
        - DPI: 300 (alta resolución)
        - bbox_inches: 'tight' (márgenes optimizados)
        - Formato: PNG con compresión optimizada
        
        Examples:
            >>> PyMutation.configure_high_quality_plots()  # Configurar una vez
            >>> py_mut = PyMutation(data)
            >>> fig = py_mut.summary_plot()
            >>> fig.savefig('plot.png')  # ¡Automáticamente alta calidad!
        
        Note:
            Esta configuración afecta a TODAS las figuras de matplotlib en la sesión.
            Se recomienda llamar esta función al inicio del script.
        """
        import matplotlib as mpl
        
        # Configurar DPI por defecto para alta resolución
        mpl.rcParams['figure.dpi'] = 300
        mpl.rcParams['savefig.dpi'] = 300
        
        # Configurar márgenes automáticos
        mpl.rcParams['savefig.bbox'] = 'tight'
        
        # Configurar formato y compresión
        mpl.rcParams['savefig.format'] = 'png'
        mpl.rcParams['savefig.transparent'] = False
        
        # Mejorar calidad de texto
        mpl.rcParams['savefig.facecolor'] = 'white'
        mpl.rcParams['savefig.edgecolor'] = 'none'
        
        print("✅ Configuración de alta calidad activada para matplotlib")
        print("   • DPI: 300 (alta resolución)")
        print("   • Márgenes: automáticos (tight)")
        print("   • Formato: PNG optimizado")
        print("   ℹ️  Ahora todas las figuras se guardarán automáticamente en alta calidad")
    
    def summary_plot(self, 
                   figsize: Tuple[int, int] = DEFAULT_SUMMARY_FIGSIZE,
                   title: str = DEFAULT_PLOT_TITLE,
                   max_samples: Optional[int] = 200,
                   top_genes_count: int = DEFAULT_TOP_GENES_COUNT,
                   show_interactive: bool = False) -> plt.Figure:
        """
        Genera un gráfico de resumen con estadísticas generales de las mutaciones.
        
        Esta visualización incluye múltiples gráficos:
        - Variant Classification: Distribución de clasificaciones de variantes
        - Variant Type: Distribución de tipos de variantes (SNP, INS, DEL, etc.)
        - SNV Class: Distribución de clases de SNV (cambios nucleotídicos como A>G, C>T, etc.)
        - Variants per Sample: Distribución de variantes por muestra y mediana (TMB)
        - Top Mutated Genes: Genes más frecuentemente mutados
        
        Args:
            figsize: Tamaño de la figura.
            title: Título del gráfico.
            max_samples: Número máximo de muestras a mostrar en el gráfico de variantes por muestra.
                        Si es None, se muestran todas las muestras.
            top_genes_count: Número de genes a mostrar en el gráfico de genes más mutados.
                        Si hay menos genes que este número, se mostrarán todos.
            show_interactive: Si es True, muestra la visualización en modo interactivo.
            
        Returns:
            Figura de matplotlib con el gráfico de resumen.
        """
        from .visualizations.summary import create_summary_plot
        from .utils.data_processing import extract_variant_classifications, extract_variant_types
        
        # Preprocesar los datos para asegurar que tenemos las columnas necesarias
        processed_data = extract_variant_classifications(
            self.data, 
            variant_column=VARIANT_CLASSIFICATION_COLUMN,
            funcotation_column=FUNCOTATION_COLUMN
        )
        
        processed_data = extract_variant_types(
            processed_data,
            variant_column=VARIANT_TYPE_COLUMN,
            funcotation_column=FUNCOTATION_COLUMN
        )
        
        # Generar el gráfico de resumen
        fig = create_summary_plot(processed_data, figsize, title, max_samples, top_genes_count)
        
        # Si se solicita mostrar interactivamente
        if show_interactive:
            self._show_figure_interactive(fig)
        
        return fig
    
    def variant_classification_plot(self,
                                    figsize: Tuple[int, int] = DEFAULT_PLOT_FIGSIZE,
                                    title: str = "Variant Classification",
                                    show_interactive: bool = False) -> plt.Figure:
        """
        Genera un gráfico de barras horizontal mostrando la distribución de clasificaciones de variantes.
        
        Args:
            figsize: Tamaño de la figura.
            title: Título del gráfico.
            show_interactive: Si es True, muestra la visualización en modo interactivo.
            
        Returns:
            Figura de matplotlib con el gráfico de clasificación de variantes.
        """
        from .visualizations.summary import create_variant_classification_plot
        from .utils.data_processing import extract_variant_classifications
        
        # Preprocesar los datos para asegurar que tenemos la columna necesaria
        processed_data = extract_variant_classifications(
            self.data, 
            variant_column="Variant_Classification",
            funcotation_column="FUNCOTATION"
        )
        
        # Crear figura y ejes
        fig, ax = plt.subplots(figsize=figsize)
        
        # Generar el gráfico, pasando set_title=False para evitar título duplicado
        create_variant_classification_plot(processed_data, ax=ax, set_title=False)
        
        # Configurar título
        if title:
            fig.suptitle(title, fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        
        # Si se solicita mostrar interactivamente
        if show_interactive:
            self._show_figure_interactive(fig)
        
        return fig
    
    def variant_type_plot(self,
                          figsize: Tuple[int, int] = DEFAULT_PLOT_FIGSIZE,
                          title: str = "Variant Type",
                          show_interactive: bool = False) -> plt.Figure:
        """
        Genera un gráfico de barras horizontal mostrando la distribución de tipos de variantes.
        
        Args:
            figsize: Tamaño de la figura.
            title: Título del gráfico.
            show_interactive: Si es True, muestra la visualización en modo interactivo.
            
        Returns:
            Figura de matplotlib con el gráfico de tipos de variantes.
        """
        from .visualizations.summary import create_variant_type_plot
        from .utils.data_processing import extract_variant_types
        
        # Preprocesar los datos para asegurar que tenemos la columna necesaria
        processed_data = extract_variant_types(
            self.data,
            variant_column="Variant_Type",
            funcotation_column="FUNCOTATION"
        )
        
        # Crear figura y ejes
        fig, ax = plt.subplots(figsize=figsize)
        
        # Generar el gráfico, pasando set_title=False para evitar título duplicado
        create_variant_type_plot(processed_data, ax=ax, set_title=False)
        
        # Configurar título
        if title:
            fig.suptitle(title, fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        
        # Si se solicita mostrar interactivamente
        if show_interactive:
            self._show_figure_interactive(fig)
        
        return fig
    
    def snv_class_plot(self,
                        figsize: Tuple[int, int] = DEFAULT_PLOT_FIGSIZE,
                        title: str = "SNV Class",
                        ref_column: str = "REF",
                        alt_column: str = "ALT",
                        show_interactive: bool = False) -> plt.Figure:
        """
        Genera un gráfico de barras horizontal mostrando la distribución de clases de SNV.
        
        Args:
            figsize: Tamaño de la figura.
            title: Título del gráfico.
            ref_column: Nombre de la columna que contiene el alelo de referencia.
            alt_column: Nombre de la columna que contiene el alelo alternativo.
            show_interactive: Si es True, muestra la visualización en modo interactivo.
            
        Returns:
            Figura de matplotlib con el gráfico de clases de SNV.
        """
        from .visualizations.summary import create_snv_class_plot
        
        fig, ax = plt.subplots(figsize=figsize)
        create_snv_class_plot(
            self.data, 
            ref_column=ref_column,
            alt_column=alt_column,
            ax=ax,
            set_title=False  # Evitar título duplicado
        )
        
        if title:
            fig.suptitle(title, fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        # Si se solicita mostrar interactivamente
        if show_interactive:
            self._show_figure_interactive(fig)
        
        return fig
    
    def variants_per_sample_plot(self,
                                 figsize: Tuple[int, int] = DEFAULT_PLOT_FIGSIZE,
                                 title: str = "Variants per Sample",
                                 variant_column: str = "Variant_Classification",
                                 sample_column: str = "Tumor_Sample_Barcode",
                                 max_samples: Optional[int] = 200,
                                 show_interactive: bool = False) -> plt.Figure:
        """
        Genera un gráfico de barras apiladas mostrando el número de variantes por muestra (TMB)
        y su composición por tipo de variante.

        Args:
            figsize: Tamaño de la figura.
            title: Título del gráfico.
            variant_column: Nombre de la columna que contiene la clasificación de variante.
            sample_column: Nombre de la columna que contiene el identificador de la muestra,
                          o string que se usará para identificar columnas de muestra si las
                          muestras están como columnas.
            max_samples: Número máximo de muestras a mostrar. Si es None, se muestran todas.
                        Si hay más muestras que este número, solo se muestran las primeras max_samples.
            show_interactive: Si es True, muestra la visualización en modo interactivo.
            
        Returns:
            Figura de matplotlib con el gráfico de variantes por muestra.
        """
        from .visualizations.summary import create_variants_per_sample_plot
        from .utils.data_processing import extract_variant_classifications

        # Si variant_column no está en las columnas, intentar normalizarlo
        if variant_column not in self.data.columns:
            # Comprobar si hay una versión con diferente capitalización
            column_lower = variant_column.lower()
            for col in self.data.columns:
                if col.lower() == column_lower:
                    variant_column = col
                    break
        
        # Asegurar que la columna de clasificación de variantes existe o se extrae
        processed_data = extract_variant_classifications(
            self.data, 
            variant_column=variant_column,
            funcotation_column="FUNCOTATION"
        )
        
        fig, ax = plt.subplots(figsize=figsize)
        create_variants_per_sample_plot(
            processed_data, 
            variant_column=variant_column,
            sample_column=sample_column,
            ax=ax,
            set_title=False,  # Evitar título duplicado
            max_samples=max_samples  # Pasar el límite de muestras configurado
        )
        
        # No modificar el título si contiene la mediana
        if title and not title.startswith("Variants per Sample"):
            fig.suptitle(title, fontsize=16, fontweight='bold', y=1.02)
        elif title:
            fig.suptitle(title, fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        
        # Si se solicita mostrar interactivamente
        if show_interactive:
            self._show_figure_interactive(fig)
        
        return fig
    
    
    def variant_classification_summary_plot(self,
                                           figsize: Tuple[int, int] = DEFAULT_PLOT_FIGSIZE,
                                           title: str = "Variant Classification Summary",
                                           variant_column: str = "Variant_Classification",
                                           sample_column: str = "Tumor_Sample_Barcode",
                                           show_interactive: bool = False) -> plt.Figure:
        """
        Genera un diagrama de cajas y bigotes (boxplot) que resume, para cada clasificación de variantes,
        la distribución (entre las muestras) del número de alelos alternativos detectados.

        Este gráfico muestra la variabilidad entre muestras para cada tipo de clasificación de variante,
        permitiendo identificar cuáles presentan más diferencias entre pacientes.

        Args:
            figsize: Tamaño de la figura.
            title: Título del gráfico.
            variant_column: Nombre de la columna que contiene la clasificación de variante.
            sample_column: Nombre de la columna que contiene el identificador de la muestra.
                          Si no existe, se asume que las muestras son columnas (formato ancho).
            show_interactive: Si es True, muestra la visualización en modo interactivo.
            
        Returns:
            Figura de matplotlib con el gráfico de cajas y bigotes.
        """
        from .visualizations.summary import create_variant_classification_summary_plot
        from .utils.data_processing import extract_variant_classifications

        # Asegurar que la columna de clasificación de variantes existe o se extrae
        processed_data = extract_variant_classifications(
            self.data, 
            variant_column=variant_column,
            funcotation_column="FUNCOTATION"
        )
        
        # Verificar si estamos en formato ancho (muestras como columnas)
        is_wide_format = sample_column not in processed_data.columns
        if is_wide_format:
            # Detectar y mostrar información sobre el formato
            sample_cols = [col for col in processed_data.columns if col.startswith('TCGA-') or 
                           (isinstance(col, str) and col.count('-') >= 2)]
            if sample_cols:
                print(f"Detected wide format with {len(sample_cols)} possible sample columns.")
        
        fig, ax = plt.subplots(figsize=figsize)
        create_variant_classification_summary_plot(
            processed_data, 
            variant_column=variant_column,
            sample_column=sample_column,
            ax=ax,
            show_labels=True,  # Asegurarnos de que siempre muestre las etiquetas cuando se genera individualmente
            set_title=False  # Evitar título duplicado
        )
        
        # Configurar título
        if title:
            fig.suptitle(title, fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        
        # Si se solicita mostrar interactivamente
        if show_interactive:
            self._show_figure_interactive(fig)
        
        return fig


    def top_mutated_genes_plot(self,
                              figsize: Tuple[int, int] = DEFAULT_PLOT_FIGSIZE,
                              title: str = "Top Mutated Genes",
                              mode: str = MODE_VARIANTS, 
                              variant_column: str = VARIANT_CLASSIFICATION_COLUMN,
                              gene_column: str = GENE_COLUMN,
                              sample_column: str = SAMPLE_COLUMN,
                              count: int = DEFAULT_TOP_GENES_COUNT,
                              show_interactive: bool = False) -> plt.Figure:
        """
        Genera un diagrama de barras horizontal mostrando los genes más mutados y la distribución
        de variantes según su clasificación.

        Args:
            figsize: Tamaño de la figura.
            title: Título del gráfico.
            mode: Modo de conteo de mutaciones: "variants" (cuenta número total de variantes)
                  o "samples" (cuenta número de muestras afectadas).
            variant_column: Nombre de la columna que contiene la clasificación de variante.
            gene_column: Nombre de la columna que contiene el símbolo del gen.
            sample_column: Nombre de la columna que contiene el identificador de la muestra,
                          o prefijo para identificar columnas de muestra si están como columnas.
            count: Número de genes principales a mostrar.
            show_interactive: Si es True, muestra la visualización en modo interactivo.
            
        Returns:
            Figura de matplotlib con el gráfico de genes más mutados.
        
        Raises:
            ValueError: Si 'count' no es un número positivo o 'mode' no es un valor válido.
        """
        from .visualizations.summary import create_top_mutated_genes_plot
        from .utils.data_processing import extract_variant_classifications

        # Validar parámetros
        if not isinstance(count, int):
            raise ValueError(f"El parámetro 'count' debe ser un número entero, se recibió: {count}")
        if count <= 0:
            raise ValueError(f"El parámetro 'count' debe ser un número entero positivo, se recibió: {count}")
        
        # Verificar que el modo sea válido
        if mode not in VALID_PLOT_MODES:
            raise ValueError(f"El modo '{mode}' no es válido. Los valores permitidos son: {', '.join(VALID_PLOT_MODES)}")

        # Si variant_column no está en las columnas, intentar normalizarlo
        if variant_column not in self.data.columns:
            # Comprobar si hay una versión con diferente capitalización
            column_lower = variant_column.lower()
            for col in self.data.columns:
                if col.lower() == column_lower:
                    variant_column = col
                    break
                    
        # Si gene_column no está en las columnas, intentar normalizarlo
        if gene_column not in self.data.columns:
            # Comprobar si hay una versión con diferente capitalización
            column_lower = gene_column.lower()
            for col in self.data.columns:
                if col.lower() == column_lower:
                    gene_column = col
                    break
        
        # Asegurar que la columna de clasificación de variantes existe o se extrae
        processed_data = extract_variant_classifications(
            self.data, 
            variant_column=variant_column,
            funcotation_column=FUNCOTATION_COLUMN
        )
        
        fig, ax = plt.subplots(figsize=figsize)
        create_top_mutated_genes_plot(
            processed_data, 
            mode=mode,
            variant_column=variant_column,
            gene_column=gene_column,
            sample_column=sample_column,
            count=count,
            ax=ax,
            set_title=False  # Evitar título duplicado
        )
        
        # Ajustar título personalizado basado en el modo
        if title:
            if mode == "variants" and title == "Top Mutated Genes":
                fig.suptitle("Top mutated genes (variants)", fontsize=16, fontweight='bold', y=0.98)
            elif mode == "samples" and title == "Top Mutated Genes":
                fig.suptitle("Top mutated genes (samples)", fontsize=16, fontweight='bold', y=0.98)
            else:
                fig.suptitle(title, fontsize=16, fontweight='bold', y=0.98)
        
        # Usar tight_layout con un padding adicional para mejorar los márgenes
        plt.tight_layout(pad=1.2)
        
        # Ajustar los márgenes para tener un aspecto más consistente entre modes
        # Aumentar el margen izquierdo para evitar que el texto se salga
        plt.subplots_adjust(left=0.15, right=0.9)
        
        # Si se solicita mostrar interactivamente
        if show_interactive:
            self._show_figure_interactive(fig)
        
        return fig
    
    def _show_figure_interactive(self, figure: plt.Figure) -> None:
        """
        Muestra una figura específica en modo interactivo sin afectar otras figuras.
        
        Este método privado se utiliza internamente para mostrar solo la figura específica
        en modo interactivo cuando se usa show_interactive=True en los métodos de visualización.
        
        Args:
            figure: La figura específica a mostrar en modo interactivo.
        """
        # Guardar el estado actual del modo interactivo
        was_interactive = plt.isinteractive()
        
        # Variable para controlar cuándo se cierra la ventana
        window_closed = [False]  # Usar lista para que sea mutable en la función interna
        
        def on_close(event):
            """Callback que se ejecuta cuando se cierra la ventana."""
            window_closed[0] = True
        
        try:
            # Habilitar modo interactivo temporalmente solo si no estaba activo
            if not was_interactive:
                plt.ion()
            
            # Conectar el evento de cierre de la ventana
            figure.canvas.mpl_connect('close_event', on_close)
            
            # Mostrar solo esta figura específica
            figure.show()
            
            # Forzar el render inmediato de la figura
            figure.canvas.draw()
            figure.canvas.flush_events()
            
            # Informar al usuario
            title_text = "Sin título"
            if figure._suptitle and figure._suptitle.get_text():
                title_text = figure._suptitle.get_text()
            
            print(f"Figura '{title_text}' mostrada en modo interactivo.")
            print("Cierra la ventana para continuar con la ejecución del script.")
            
            # Esperar hasta que el usuario cierre la ventana
            # Usar un loop que verifique la variable window_closed
            import time
            while not window_closed[0] and plt.fignum_exists(figure.number):
                # Usar solo time.sleep() sin plt.pause() para evitar efectos secundarios
                time.sleep(0.1)
                # Permitir que matplotlib procese eventos de forma mínima
                figure.canvas.flush_events()
                    
        finally:
            # Restaurar el estado original del modo interactivo solo si lo cambiamos
            if not was_interactive:
                plt.ioff()  # Desactivar modo interactivo si no estaba activo antes

