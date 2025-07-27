import logging
from pathlib import Path
from typing import List, Optional, Union
import pandas as pd
import numpy as np

from .input import required_columns_MAF
from .utils.format import reverse_format_chr

# Importar PyArrow si está disponible
try:
    import pyarrow as pa
    import pyarrow.compute as pc
    HAS_PYARROW = True
except ImportError:
    HAS_PYARROW = False

logger = logging.getLogger(__name__)

# Variable global para cachear el orden de columnas
_MAF_COLUMN_ORDER_CACHE: Optional[List[str]] = None

def _load_maf_column_order() -> List[str]:
    """
    Load the MAF column order from MAF_COL_ORDER.csv file.
    Uses a cached version if available to avoid reading the file multiple times.

    Returns
    -------
    List[str]
        List of column names in the order specified in MAF_COL_ORDER.csv
    """
    global _MAF_COLUMN_ORDER_CACHE
    
    # Usar caché si ya existe
    if _MAF_COLUMN_ORDER_CACHE is not None:
        return _MAF_COLUMN_ORDER_CACHE
        
    maf_columns_file = Path(__file__).parent / "data" / "MAF_COL_ORDER.csv"

    try:
        if maf_columns_file.exists():
            columns_df = pd.read_csv(maf_columns_file)
            # Extract column names from the 'nombre' column, ordered by 'id'
            ordered_columns = columns_df.sort_values('id')['nombre'].tolist()
            logger.debug(f"Loaded {len(ordered_columns)} column names from MAF_COL_ORDER.csv")
            
            # Guardar en caché
            _MAF_COLUMN_ORDER_CACHE = ordered_columns
            return ordered_columns
        else:
            logger.warning(f"MAF_COL_ORDER.csv not found at {maf_columns_file}, using default order")
            return []
    except Exception as e:
        logger.warning(f"Error loading MAF_COL_ORDER.csv: {e}, using default order")
        return []


def to_maf(self, output_path: str | Path) -> None:
    """
    Export a PyMutation object back to MAF format.

    This function reverses the transformations done by read_maf() to recreate
    a MAF file from a PyMutation object. The output follows the column order
    specified in MAF_COL_ORDER.csv when available.

    Parameters
    ----------
    output_path : str | Path
        Path where the MAF file will be written.

    Raises
    ------
    ValueError
        If the PyMutation object doesn't contain the necessary data for MAF export.
    """
    output_path = Path(output_path)
    logger.info("Starting MAF export to: %s", output_path)

    # Get the data and samples from PyMutation object
    data = self.data.copy()
    samples = self.samples
    metadata = self.metadata
    
    # Log total number of variants to process
    total_variants = len(data)
    logger.info(f"Starting to process {total_variants} variants from {len(samples)} samples")

    # ─── 1) VALIDATE REQUIRED COLUMNS FOR EXPORT ──────────────────────
    vcf_like_cols = ["CHROM", "POS", "REF", "ALT","ID"]
    missing_vcf_cols = [col for col in vcf_like_cols if col not in data.columns]
    if missing_vcf_cols:
        raise ValueError(f"Missing required VCF-style columns for MAF export: {missing_vcf_cols}")

    missing_samples = [sample for sample in samples if sample not in data.columns]
    if missing_samples:
        raise ValueError(f"Missing sample columns for MAF export: {missing_samples}")

    # ─── 2) CONVERT TO MAF FORMAT ────────────────────────────────
    # Determinar si usar PyArrow para conjuntos de datos grandes
    use_pyarrow = HAS_PYARROW and len(data) > 10000
    
    if use_pyarrow:
        logger.info(f"Using PyArrow for processing large dataset ({len(data)} rows)")
        
        # Convertir a tabla PyArrow para procesamiento más rápido
        table = pa.Table.from_pandas(data)
        
        # Crear columnas base usando PyArrow
        # Convertir CHROM a Chromosome - solo extraer el número del cromosoma
        chrom_array = pc.cast(table['CHROM'], pa.string())
        # Eliminar el prefijo 'chr' si existe
        chrom_array = pc.utf8_replace_substring(chrom_array, 'chr', '')
        table = table.append_column('Chromosome', chrom_array)
        
        # Otras columnas base
        table = table.append_column('Start_Position', table['POS'])
        table = table.append_column('Reference_Allele', pc.cast(table['REF'], pa.string()))
        table = table.append_column('NCBI_Build', pa.array([metadata.assembly] * len(table)))
        table = table.append_column('dbSNP_RS', pc.cast(table['ID'], pa.string()))
        
        # Convertir de vuelta a pandas para operaciones específicas
        base_data = table.to_pandas()
    else:
        # Crear columnas base para todas las variantes usando pandas
        base_data = data.copy()
        # Eliminar el prefijo 'chr' si existe
        base_data['Chromosome'] = base_data['CHROM'].str.replace('chr', '', regex=False)
        base_data['Start_Position'] = base_data['POS']
        base_data['Reference_Allele'] = base_data['REF'].astype(str)
        base_data['NCBI_Build'] = metadata.assembly
        base_data['dbSNP_RS'] = base_data['ID'].astype(str)
    
    # Procesar cada muestra y crear un DataFrame por muestra
    sample_dfs = []
    processed_samples = 0
    total_samples = len(samples)
    log_frequency = max(1, total_samples // 50)  # Log cada ~2% de las muestras, mínimo 1
    
    for sample_idx, sample in enumerate(samples, 1):
        # Log progress for sample processing (solo cada log_frequency muestras)
        if sample_idx % log_frequency == 0 or sample_idx == 1 or sample_idx == total_samples:
            logger.info(f"Processing sample {sample_idx}/{total_samples}: {sample} ({sample_idx/total_samples*100:.1f}%)")
        
        # Crear una copia de los datos para esta muestra
        sample_data = base_data.copy()
        sample_data['Tumor_Sample_Barcode'] = sample
        
        # Filtrar solo variantes donde el genotipo no es REF|REF
        ref_pattern = sample_data['REF'] + '|' + sample_data['REF']
        sample_data = sample_data[sample_data[sample] != ref_pattern]
        
        # Incrementar contador de muestras procesadas
        processed_samples += 1
        
        # Log number of variants for this sample (solo cada log_frequency muestras)
        if sample_idx % log_frequency == 0 or sample_idx == 1 or sample_idx == total_samples:
            sample_variants = len(sample_data)
            logger.info(f"Sample {sample}: {sample_variants} variants found")
        
        if len(sample_data) == 0:
            continue
        
        # Procesar genotipos de manera vectorizada
        if '|' in sample_data[sample].iloc[0]:
            # Separar por '|'
            sample_data[['Tumor_Seq_Allele1', 'Tumor_Seq_Allele2']] = sample_data[sample].str.split('|', expand=True)
        elif '/' in sample_data[sample].iloc[0]:
            # Separar por '/'
            sample_data[['Tumor_Seq_Allele1', 'Tumor_Seq_Allele2']] = sample_data[sample].str.split('/', expand=True)
        else:
            # Caso especial: usar REF y ALT
            sample_data['Tumor_Seq_Allele1'] = sample_data['REF']
            sample_data['Tumor_Seq_Allele2'] = sample_data['ALT']
        
        # Calcular End_Position de manera vectorizada
        # Para SNPs (ref y alt de longitud 1)
        snp_mask = (sample_data['Reference_Allele'].str.len() == 1) & (sample_data['Tumor_Seq_Allele2'].str.len() == 1)
        sample_data.loc[snp_mask, 'End_Position'] = sample_data.loc[snp_mask, 'Start_Position'].astype(int)
        
        # Para deleciones (ref más largo que alt)
        del_mask = sample_data['Reference_Allele'].str.len() > sample_data['Tumor_Seq_Allele2'].str.len()
        sample_data.loc[del_mask, 'End_Position'] = (sample_data.loc[del_mask, 'Start_Position'] + sample_data.loc[del_mask, 'Reference_Allele'].str.len() - 1).astype(int)
        
        # Para inserciones y otros
        other_mask = ~(snp_mask | del_mask)
        sample_data.loc[other_mask, 'End_Position'] = sample_data.loc[other_mask, 'Start_Position'].astype(int)
        
        # Seleccionar columnas relevantes
        maf_cols = ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 
                    'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode', 
                    'NCBI_Build', 'dbSNP_RS']
        
        # Añadir otras columnas que no son VCF ni muestras
        other_cols = [col for col in data.columns if col not in vcf_like_cols + samples]
        
        # Combinar todas las columnas
        all_cols = maf_cols + other_cols
        
        # Seleccionar solo las columnas que existen
        existing_cols = [col for col in all_cols if col in sample_data.columns]
        
        sample_dfs.append(sample_data[existing_cols])
    
    # Combinar todos los DataFrames de muestras
    if not sample_dfs:
        logger.warning("No variant data found to export")
        maf_df = pd.DataFrame(columns=required_columns_MAF)
    else:
        maf_df = pd.concat(sample_dfs, ignore_index=True)
        
    # Log resumen después de procesar todas las muestras
    total_variants_found = len(maf_df)
    logger.info(f"Sample processing completed: {processed_samples}/{total_samples} samples processed")
    logger.info(f"Total variants found: {total_variants_found} variants")

    # ─── 3) ENSURE COLUMN ORDER ─────────────────────

    # Load the preferred column order from MAF_COL_ORDER.csv
    preferred_column_order = _load_maf_column_order()

    # Ensure required columns are present
    for col in required_columns_MAF:
        if col not in maf_df.columns:
            # Add missing required columns with default values
            maf_df[col] = "."

    # Optimized column ordering
    if preferred_column_order:
        # Filtrar solo las columnas que existen en el DataFrame
        existing_preferred_columns = [col for col in preferred_column_order if col in maf_df.columns]
        
        # Obtener columnas que no están en el orden preferido
        remaining_columns = list(set(maf_df.columns) - set(existing_preferred_columns))
        
        # Combinar las listas
        final_columns = existing_preferred_columns + remaining_columns
        
        logger.info(f"Using MAF_COL_ORDER.csv column order: {len(final_columns)} columns arranged")
        
        # Reordenar el DataFrame
        maf_df = maf_df[final_columns]
    
    # Ensure End_Position is an integer
    if 'End_Position' in maf_df.columns:
        maf_df['End_Position'] = maf_df['End_Position'].astype(int)

    # ─── 4) WRITE TO FILE WITH COMMENTS ──────────────────────────────
    # Definir tamaño de lote para escritura eficiente
    chunk_size = 10000  # Ajustar según necesidades y memoria disponible
    
    with open(output_path, 'w', encoding='utf-8') as f:
        # Write only INFO comments from metadata if they exist
        if metadata.notes:
            for line in metadata.notes.split('\n'):
                if line.strip():
                    # Only include lines that start with INFO or ##INFO
                    if line.startswith('##INFO') or line.startswith('INFO'):
                        if not line.startswith('#'):
                            f.write(f"#{line}\n")
                        else:
                            f.write(f"{line}\n")

        # Escribir encabezado
        f.write('\t'.join(maf_df.columns) + '\n')
        
        # Para conjuntos de datos pequeños, escribir directamente
        if len(maf_df) <= chunk_size:
            logger.info(f"Writing {len(maf_df)} variants to file")
            maf_df.to_csv(f, sep='\t', index=False, header=False, lineterminator='\n')
            logger.info(f"Progress: {len(maf_df)}/{len(maf_df)} variants written (100.0%)")
        else:
            # Para conjuntos de datos grandes, escribir en lotes
            total_rows = len(maf_df)
            logger.info(f"Writing large dataset ({total_rows} variants) in chunks of {chunk_size}")
            
            for i in range(0, total_rows, chunk_size):
                chunk = maf_df.iloc[i:i+chunk_size]
                chunk.to_csv(f, sep='\t', index=False, header=False, lineterminator='\n')
                
                # Calcular el número real de filas escritas (puede ser menor que i+chunk_size en el último lote)
                rows_written = min(i + chunk_size, total_rows)
                
                # Registrar progreso para cada lote
                logger.info(f"Progress: {rows_written}/{total_rows} variants written ({rows_written/total_rows*100:.1f}%)")
        
        # Eliminar el último salto de línea si existe
        if f.tell() > 0:  # Asegurarse de que el archivo no está vacío
            # Verificar si el último carácter es un salto de línea
            f.seek(f.tell() - 1)
            # No podemos leer en modo 'w', así que asumimos que hay un salto de línea
            # y lo eliminamos directamente
            f.truncate()

    # Log final summary with detailed statistics
    logger.info(f"MAF export completed successfully: {len(maf_df)} variants processed and written to {output_path}")
    logger.info(f"Conversion summary: {len(samples)} samples, {total_variants} input variants, {len(maf_df)} output variants")


# Add the method to PyMutation class
def add_to_maf_method_to_pymutation():
    from .core import PyMutation
    PyMutation.to_maf = to_maf

add_to_maf_method_to_pymutation()
