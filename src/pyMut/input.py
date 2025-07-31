import gzip
import hashlib
import io
import logging
import os
import subprocess
import time
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd

from .core import PyMutation, MutationMetadata
from .utils.format import format_rs, format_chr, normalize_variant_classification

# Optional imports
try:
    import pyarrow as pa
    import pyarrow.compute as pc

    HAS_PYARROW = True
except ImportError:
    HAS_PYARROW = False

try:
    import cyvcf2

    HAS_CYVCF2 = True
except ImportError:
    HAS_CYVCF2 = False

try:
    import psutil

    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False

# Logger configuration
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
if not logger.handlers:
    logging.basicConfig(
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        level=logging.INFO,
    )

# Column validation
required_columns_VCF: List[str] = ["CHROM", "POS", "ID", "REF", "ALT", "FILTER"]
_required_canonical_VCF = {c.lower(): c for c in required_columns_VCF}

RANK_SEVERE = [
    # HIGH
    "stop_gained", "stop_lost", "frameshift_variant",
    "start_lost",
    "splice_donor_variant", "splice_acceptor_variant",
    "splice_donor_region_variant", "splice_donor_5th_base_variant",
    "splice_region_variant", "splice_polypyrimidine_tract_variant",
    # MODERATE
    "inframe_deletion", "inframe_insertion",
    "missense_variant", "protein_altering_variant",
    "synonymous_variant", "stop_retained_variant", "start_retained_variant",
    "incomplete_terminal_codon_variant",
    # LOW/MODIFIER – UTRs, CDS, etc.
    "5_prime_UTR_variant", "3_prime_UTR_variant", "coding_sequence_variant",
    "intron_variant", "non_coding_transcript_exon_variant",
    "non_coding_transcript_variant", "NMD_transcript_variant",
    "mature_miRNA_variant",
    "upstream_gene_variant", "downstream_gene_variant",
    "intergenic_variant",
    # Reguladoras
    "regulatory_region_variant", "TF_binding_site_variant"
]
rank_ser = pd.Series(range(len(RANK_SEVERE)), index=RANK_SEVERE)

MAP_SNV = {
    "stop_gained": "Nonsense_Mutation",
    "stop_lost": "Nonstop_Mutation",
    "stop_retained_variant": "Silent",
    "missense_variant": "Missense_Mutation",
    "synonymous_variant": "Silent",
    "splice_donor_variant": "Splice_Site",
    "splice_acceptor_variant": "Splice_Site",
    "splice_donor_region_variant": "Splice_Site",
    "splice_donor_5th_base_variant": "Splice_Site",
    "splice_region_variant": "Splice_Region",
    "splice_polypyrimidine_tract_variant": "Splice_Region",
    "5_prime_UTR_variant": "5'UTR",
    "3_prime_UTR_variant": "3'UTR",
    "upstream_gene_variant": "5'Flank",
    "downstream_gene_variant": "3'Flank",
    "intron_variant": "Intron",
    "intergenic_variant": "IGR",
    # todo: RNA / non‑coding ya se gestionan más abajo
}

MAP_INDEL = {
    "stop_gained": "FrameShift",  # se completará con _Del o _Ins
    "stop_lost": "InFrame",
    "missense_variant": "InFrame",
    # los demás (UTR, intrón, etc.) toman la misma etiqueta que SNV
}

required_columns_MAF: List[str] = [
    "Chromosome",
    "Start_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode",
]
_required_canonical_MAF = {c.lower(): c for c in required_columns_MAF}


# Utility functions

def _get_cache_path(file_path: Path, cache_dir: Optional[Path] = None, **config_params) -> Path:
    """Generate cache file path based on file hash and configuration parameters."""
    if cache_dir is None:
        cache_dir = file_path.parent / ".pymut_cache"

    cache_dir.mkdir(parents=True, exist_ok=True)

    # Create hash from file path, modification time, size, and configuration parameters
    file_stat = file_path.stat()
    hash_input = f"{file_path}_{file_stat.st_mtime}_{file_stat.st_size}"
    
    # Add configuration parameters to hash
    if config_params:
        config_str = "_".join(f"{k}={v}" for k, v in sorted(config_params.items()))
        hash_input += f"_{config_str}"
    
    file_hash = hashlib.md5(hash_input.encode()).hexdigest()[:16]

    return cache_dir / f"{file_path.stem}_{file_hash}.parquet"


def _create_tabix_index(vcf_path: Path, create_index: bool = False) -> bool:
    """Create Tabix index if it doesn't exist and create_index is True."""
    if not create_index:
        return False

    tbi_path = Path(str(vcf_path) + ".tbi")
    if tbi_path.exists():
        return True

    try:
        logger.info("Creating Tabix index for %s", vcf_path)
        subprocess.run(["tabix", "-p", "vcf", str(vcf_path)], check=True, capture_output=True)
        logger.info("Tabix index created successfully")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.warning("Failed to create Tabix index: %s", e)
        return False


def _get_primary_consequence(consequence_series: pd.Series) -> pd.Series:
    """
    Get the primary (most severe) consequence from VEP Consequence column using vectorized operations.
    
    Parameters
    ----------
    consequence_series : pd.Series
        Series containing VEP Consequence values (may contain multiple consequences separated by '&')
    
    Returns
    -------
    pd.Series
        Series with the most severe consequence for each variant
    """
    # explode separa por '&' y duplica las filas manteniendo su índice original
    tmp = (
        consequence_series
        .str.split('&')
        .explode()
        .str.strip()  # quita espacios
        .to_frame('term')
    )
    tmp['rank'] = tmp['term'].map(rank_ser)

    # elegimos la consecuencia con rank más bajo para cada variante
    primary = (
        tmp.reset_index()  # 'index' ← índice original de df
        .sort_values(['index', 'rank'])
        .drop_duplicates('index')
        .set_index('index')['term']
    )

    return primary


def _classify_variant(row):
    """
    Classify a single variant based on primary consequence and variant class.
    
    Parameters
    ----------
    row : pd.Series
        Row containing 'prim_cons', 'VEP_VARIANT_CLASS', 'REF', 'ALT' columns
    
    Returns
    -------
    str
        Variant classification in MAF format
    """
    cons = row['prim_cons']
    vtype = row['VEP_VARIANT_CLASS']

    # --- SNVs ------------------------------------------------------------
    if vtype == 'SNV':
        return MAP_SNV.get(cons, 'RNA')

    # --- indels ----------------------------------------------------------
    elif vtype in ('deletion', 'insertion'):
        # longitud REF y ALT deben estar en tu DataFrame
        indel_len = abs(len(str(row['REF'])) - len(str(row['ALT'])))
        in_frame = (indel_len % 3 == 0)

        base = MAP_INDEL.get(cons, MAP_SNV.get(cons, 'RNA'))

        if base in ('FrameShift', 'Nonsense_Mutation'):  # stop_gained
            base = 'Frame_Shift' if not in_frame else 'In_Frame'
        elif base == 'InFrame':
            base = 'In_Frame'

        suffix = 'Del' if vtype == 'deletion' else 'Ins'
        return f'{base}_{suffix}' if not base.endswith('_') else base + suffix

    # --- otros (raro) ----------------------------------------------------
    return 'RNA'


def _map_vep_variant_class_to_maf_variant_type(row):
    """
    Map VEP variant_class to MAF Variant_Type according to translation rules.
    
    Translation rules:
    - SNV → SNP
    - substitution → ONP (with length rules for DNP/TNP)
    - insertion → INS
    - deletion → DEL
    - indel → IND (new category)
    - inversion → INV
    - copy_number_variation → CNV
    - If no rule matches → UNKNOWN
    
    Length rules for alleles (REF vs ALT):
    - If len(REF) == len(ALT) and length = 1 → SNP
    - If equal lengths ≥ 2 → substitution; subdivide as:
      • 2 bp → DNP  • 3 bp → TNP  • ≥ 4 bp → ONP
    - If len(REF) == 0 after VCF normalization → INS
    - If len(ALT) == 0 → DEL
    - If both lengths > 0 and different → true mixed indel → IND
    
    Parameters
    ----------
    row : pd.Series
        Row containing 'VEP_VARIANT_CLASS', 'REF', 'ALT' columns
    
    Returns
    -------
    str
        Variant_Type in MAF format
    """
    vep_variant_class = row.get('VEP_VARIANT_CLASS', '')
    ref = str(row.get('REF', ''))
    alt = str(row.get('ALT', ''))

    # Handle missing or empty values
    if not vep_variant_class or vep_variant_class in ['', '.', 'nan', 'None']:
        vep_variant_class = ''
    if not ref or ref in ['.', 'nan', 'None']:
        ref = ''
    if not alt or alt in ['.', 'nan', 'None']:
        alt = ''

    # Calculate lengths
    ref_len = len(ref)
    alt_len = len(alt)

    # First, handle specific VEP variant classes that should not be overridden by length rules
    vep_to_maf_map = {
        'SNV': 'SNP',
        'insertion': 'INS',
        'deletion': 'DEL',
        'indel': 'IND',
        'inversion': 'INV',
        'copy_number_variation': 'CNV'
    }

    # If we have a direct VEP mapping (except substitution), use it
    if vep_variant_class in vep_to_maf_map:
        return vep_to_maf_map[vep_variant_class]

    # Handle substitution with length-based rules
    if vep_variant_class == 'substitution':
        if ref_len == alt_len:
            if ref_len == 1:
                return 'SNP'
            elif ref_len == 2:
                return 'DNP'
            elif ref_len == 3:
                return 'TNP'
            else:  # >= 4
                return 'ONP'
        else:
            # Different lengths for substitution - treat as mixed indel
            return 'IND'

    # If no VEP class or unknown VEP class, apply length-based rules as fallback
    # This handles cases where VEP_VARIANT_CLASS is missing, empty, or not in our known mappings
    if ref_len == alt_len and ref_len == 1:
        return 'SNP'
    elif ref_len == alt_len and ref_len >= 2:
        if ref_len == 2:
            return 'DNP'
        elif ref_len == 3:
            return 'TNP'
        else:
            return 'ONP'
    elif ref_len == 0:
        return 'INS'
    elif alt_len == 0:
        return 'DEL'
    elif ref_len > 0 and alt_len > 0 and ref_len != alt_len:
        return 'IND'

    # If nothing matches, return UNKNOWN
    return 'UNKNOWN'


def _vectorized_gt_to_alleles(gt_series: pd.Series, ref_series: pd.Series, alt_series: pd.Series) -> pd.Series:
    """Convert genotypes to alleles using vectorized operations."""
    gt_array = gt_series.astype(str).values
    ref_array = ref_series.astype(str).values
    alt_array = alt_series.astype(str).values

    result = np.empty(len(gt_array), dtype=object)

    for i in range(len(gt_array)):
        gt = gt_array[i]
        ref = ref_array[i]
        alt = alt_array[i]

        gt_str = str(gt)
        if gt_str in {".", "./.", ".|.", "nan", "None", ""}:
            result[i] = gt if gt != "nan" else "."
            continue

        gt_core = gt_str.split(":", 1)[0]
        sep = "|" if "|" in gt_core else "/"
        alt_list = alt.split(",") if alt else []
        allele_indices = gt_core.replace("|", "/").split("/")
        translated = []

        for idx in allele_indices:
            if idx == ".":
                translated.append(".")
            else:
                try:
                    i_idx = int(idx)
                    if i_idx == 0:
                        translated.append(ref)
                    else:
                        translated.append(alt_list[i_idx - 1] if i_idx - 1 < len(alt_list) else ".")
                except (ValueError, IndexError):
                    translated.append(".")

        result[i] = sep.join(translated)

    return pd.Series(result, index=gt_series.index)


def _parse_info_column_vectorized(info_series: pd.Series) -> pd.DataFrame:
    """Parse INFO column using pyarrow if available."""
    if not HAS_PYARROW:
        return _parse_info_column(info_series)

    try:
        # Convert to pyarrow for faster string operations
        arrow_series = pa.array(info_series.fillna(""))
        split_info = pc.split_pattern(arrow_series, pattern=";")
        all_keys = set()
        parsed_data = []

        for i in range(len(split_info)):
            row_dict = {}
            if split_info[i].is_valid:
                for entry in split_info[i].as_py():
                    if not entry:
                        continue
                    if "=" in entry:
                        k, v = entry.split("=", 1)
                        row_dict[k] = v
                        all_keys.add(k)
                    else:
                        row_dict[entry] = True
                        all_keys.add(entry)
            parsed_data.append(row_dict)

        # Create DataFrame with all keys as columns
        result_df = pd.DataFrame(parsed_data)
        return result_df

    except Exception as e:
        logger.warning("PyArrow INFO parsing failed (%s), falling back to pandas", e)
        return _parse_info_column(info_series)


# MAF processing functions
def _standardise_vcf_columns(vcf: pd.DataFrame) -> pd.DataFrame:
    """Normaliza nombres de columnas (insensible a mayúsculas) para VCF."""
    rename_map = {
        col: _required_canonical_VCF[col.lstrip("#").lower()]
        for col in vcf.columns
        if col.lstrip("#").lower() in _required_canonical_VCF
    }
    vcf.rename(columns=rename_map, inplace=True)
    return vcf


def _standardise_maf_columns(maf: pd.DataFrame) -> pd.DataFrame:
    """
    Standardizes MAF column names in a case-insensitive manner.

    Parameters
    ----------
    maf : pd.DataFrame
        Input DataFrame containing MAF data.

    Returns
    -------
    pd.DataFrame
        The same DataFrame with columns renamed to their canonical names
        (as defined in `required_columns_MAF`).
    """
    rename_map = {
        col: _required_canonical_MAF[col.lower()]
        for col in maf.columns
        if col.lower() in _required_canonical_MAF
    }
    maf.rename(columns=rename_map, inplace=True)
    return maf


def _open_text_maybe_gzip(path: str | Path):
    """
    Return a text-mode file handle, handling both gzip-compressed and plain text files.

    Parameters
    ----------
    path : str | Path
        Path to the file to open. If the file extension is '.gz', it will be
        opened as gzip-compressed text; otherwise, it will be opened as plain text.

    Returns
    -------
    TextIO
        A file-like object opened in text mode with UTF-8 encoding, where
        encoding errors are replaced.

    Raises
    ------
    FileNotFoundError
        If the specified file does not exist.
    """
    path = Path(path)
    if path.suffix == ".gz":
        logger.debug("Opening %s as gzip-text.", path)
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    logger.debug("Opening %s as plain text.", path)
    return open(path, encoding="utf-8", errors="replace")


def _gt_to_alleles(gt: str, ref: str, alt: str) -> str:
    """
    Convert a numeric genotype string to actual allele sequences.

    Examples
    --------
    0|1  →  T|C
    1/1  →  C/C

    Parameters
    ----------
    gt : str
        Genotype string, may include phased ('|') or unphased ('/') separators,
        and may include a ':'-delimited genotype quality or other info.
    ref : str
        Reference allele sequence.
    alt : str
        Comma-separated list of alternate allele sequences.

    Returns
    -------
    str
        Genotype string with numeric indices replaced by the corresponding
        allele sequences. Returns an empty string if `gt` is empty, and
        preserves no-calls ('.', './.', '.|.') as-is.
    """
    if not gt:  # empty string
        return ""
    # Keep only the part before ':' if it exists
    gt_core = gt.split(":", 1)[0]

    # Case "." (no call)
    if gt_core in {".", "./.", ".|."}:
        return gt_core

    # Determine the separator used ('|' preferentially)
    sep = "|" if "|" in gt_core else "/"

    # List with alternative alleles (can be several)
    alt_list: list[str] = alt.split(",") if alt else []

    allele_indices = gt_core.replace("|", "/").split("/")
    translated: list[str] = []
    for idx in allele_indices:
        if idx == ".":
            translated.append(".")
            continue
        try:
            i = int(idx)
        except ValueError:
            translated.append(".")
            continue

        if i == 0:
            translated.append(ref)
        else:
            # i-1 because ALT[0] is allele «1»
            translated.append(alt_list[i - 1] if i - 1 < len(alt_list) else ".")
    return sep.join(translated)


def _parse_info_column(info_series: pd.Series) -> pd.DataFrame:
    """
    Expand the INFO column (semicolon‐delimited key–value pairs) into a DataFrame.

    Value‐less flags are marked as True. The result has the same number of rows
    as `info_series`, with each key becoming a column.

    Parameters
    ----------
    info_series : pd.Series
        Input Series containing INFO strings (e.g., "DP=100;PASS;AF=0.5") or NaN.

    Returns
    -------
    pd.DataFrame
        DataFrame where each column corresponds to a key from the INFO field.
        Values are strings for key–value pairs, and True for flags without values.
    """

    def _parser(cell: str | float) -> dict[str, str | bool]:
        if pd.isna(cell):
            return {}
        out: dict[str, str | bool] = {}
        for entry in str(cell).split(";"):
            if not entry:
                continue
            if "=" in entry:
                k, v = entry.split("=", 1)
                out[k] = v
            else:
                # Flag without value (e.g., "PASS")
                out[entry] = True
        return out

    parsed = info_series.apply(_parser)
    info_df = pd.DataFrame(parsed.tolist())  # type: ignore[arg-type]
    return info_df


# Main functions

def _get_mutated_genotype_vectorized(genotype_series, ref_allele):
    """
    Función auxiliar optimizada para encontrar genotipos mutados.
    
    Parameters
    ----------
    genotype_series : pd.Series
        Serie con genotipos para una muestra específica
    ref_allele : str
        Alelo de referencia
    
    Returns
    -------
    str
        Genotipo mutado si existe, sino genotipo de referencia
    """
    ref_genotype = f"{ref_allele}|{ref_allele}"
    
    # Filtrar genotipos no nulos
    valid_genotypes = genotype_series.dropna()
    
    # Encontrar genotipos mutados (diferentes de referencia)
    mutated = valid_genotypes[valid_genotypes != ref_genotype]
    
    # Retornar el primer genotipo mutado encontrado, o referencia si no hay
    return mutated.iloc[0] if len(mutated) > 0 else ref_genotype


def _consolidate_variants(maf_df, samples):
    """
    Versión optimizada que consolida variantes idénticas usando operaciones vectorizadas.
    
    Parameters
    ----------
    maf_df : pd.DataFrame
        DataFrame con variantes expandidas por muestra
    samples : list
        Lista de nombres de muestras
    
    Returns
    -------
    pd.DataFrame
        DataFrame consolidado con variantes únicas
    """
    import time
    start_time = time.time()
    
    # Crear una copia para no modificar el DataFrame original
    maf_work = maf_df.copy()
    
    # OPTIMIZACIÓN 1: Usar hash numérico en lugar de concatenación de strings
    # Crear clave de agrupación más eficiente
    maf_work['chrom_hash'] = pd.util.hash_pandas_object(maf_work['CHROM'], index=False)
    maf_work['variant_key'] = (
        maf_work['chrom_hash'].astype(str) + '_' + 
        maf_work['POS'].astype(str) + '_' + 
        pd.util.hash_pandas_object(maf_work['REF'], index=False).astype(str) + '_' +
        pd.util.hash_pandas_object(maf_work['ALT'], index=False).astype(str)
    )
    
    # Separar columnas de muestras de columnas de anotación
    annotation_cols = [col for col in maf_work.columns 
                      if col not in samples + ['variant_key', 'chrom_hash']]
    
    # OPTIMIZACIÓN 2: Usar operaciones vectorizadas con groupby.agg()
    logger = logging.getLogger(__name__)
    logger.info("Consolidating variants using vectorized operations...")
    
    # Consolidar columnas de anotación (tomar el primer valor de cada grupo)
    annotation_consolidated = maf_work.groupby('variant_key')[annotation_cols].first()
    
    # OPTIMIZACIÓN 3: Procesar todas las muestras de una vez usando operaciones vectorizadas
    sample_data = {}
    
    for sample in samples:
        # Para cada muestra, encontrar el genotipo no-referencia si existe
        sample_col = maf_work.groupby('variant_key')[sample].apply(
            lambda x: _get_mutated_genotype_vectorized(x, annotation_consolidated.loc[x.name, 'REF'])
        )
        sample_data[sample] = sample_col
    
    # Crear DataFrame de muestras consolidado
    samples_df = pd.DataFrame(sample_data)
    
    # OPTIMIZACIÓN 4: Combinar usando concat en lugar de iteración
    consolidated_df = pd.concat([annotation_consolidated, samples_df], axis=1)
    
    # Reset index y limpiar columnas temporales
    consolidated_df.reset_index(drop=True, inplace=True)
    
    processing_time = time.time() - start_time
    logger.info("Variant consolidation completed in %.2f seconds", processing_time)
    
    return consolidated_df


def _consolidate_variants_chunked(maf_df, samples, chunk_size=50000):
    """
    Versión que procesa en chunks para datasets muy grandes.
    """
    if len(maf_df) <= chunk_size:
        return _consolidate_variants(maf_df, samples)
    
    # Procesar en chunks
    chunks = []
    for i in range(0, len(maf_df), chunk_size):
        chunk = maf_df.iloc[i:i+chunk_size]
        consolidated_chunk = _consolidate_variants(chunk, samples)
        chunks.append(consolidated_chunk)
    
    # Combinar chunks y consolidar resultado final
    combined = pd.concat(chunks, ignore_index=True)
    return _consolidate_variants(combined, samples)


def read_maf(path: str | Path, assembly: str, cache_dir: Optional[str | Path] = None, consolidate_variants: bool = True) -> PyMutation:
    """
    Read a MAF file and return a PyMutation object with automatic caching.

    Parameters
    ----------
    path : str | Path
        Path to the MAF file (.maf or .maf.gz).
    assembly : str
        Version of the genome assembly. Must be either "37" or "38".
    cache_dir : str | Path, optional
        Directory for caching processed files. If None, uses a .pymut_cache
        directory next to the input file.
    consolidate_variants : bool, default True
        If True, consolidates identical variants that appear in multiple samples
        into a single row. If False, maintains the original behavior where each
        MAF row becomes a separate DataFrame row.

    Returns
    -------
    PyMutation
        Set of mutations read from the MAF file, converted to the same
        wide-format produced by `read_vcf`. Includes metadata with information
        about the source file, comments, and configuration.

    Raises
    ------
    FileNotFoundError
        If the specified MAF file path does not exist.
    ValueError
        If the MAF file is missing required columns or has an invalid format.
        If the assembly parameter is not "37" or "38".
    ImportError
        If the 'pyarrow' library cannot be imported and the 'c' engine
        alternative also fails.
    Exception
        For any other errors encountered while reading or processing the file.

    Notes
    -----
    Performance optimizations include:
    - Automatic caching of processed results
    - PyArrow-accelerated data processing when available
    """
    # Validate assembly parameter
    if assembly not in ["37", "38"]:
        raise ValueError("Assembly parameter must be either '37' or '38'")

    start_time = time.time()
    path = Path(path)
    cache_dir_path = Path(cache_dir) if cache_dir else None

    logger.info("Starting MAF reading: %s", path)

    # Check cache
    cache_path = _get_cache_path(path, cache_dir_path, consolidate_variants=consolidate_variants)
    if cache_path.exists():
        logger.info("Loading from cache: %s", cache_path)
        try:
            maf = pd.read_parquet(cache_path)

            # Load metadata from cache info file
            cache_info_path = cache_path.with_suffix('.json')
            if cache_info_path.exists():
                import json
                with open(cache_info_path) as f:
                    cache_info = json.load(f)
                comments = cache_info.get('comments', [])
                samples = cache_info.get('samples', [])
            else:
                comments = []
                samples = []

            metadata = MutationMetadata(
                source_format="MAF",
                file_path=str(path),
                filters=["."],
                assembly=assembly,
                notes="\n".join(comments) if comments else None,
            )

            logger.info("Cache loaded successfully in %.2f seconds", time.time() - start_time)
            return PyMutation(maf, metadata, samples)

        except Exception as e:
            logger.warning("Cache loading failed (%s), proceeding with fresh read", e)

    # Separate comments from body
    comments: list[str] = []
    buf = io.StringIO()
    try:
        with _open_text_maybe_gzip(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    comments.append(line.rstrip("\n"))
                else:
                    buf.write(line)
        buf.seek(0)
        logger.debug("Comments found: %d", len(comments))
    except FileNotFoundError:
        logger.error("File not found: %s", path)
        raise
    except Exception:
        logger.exception("Error reading header/comments from MAF.")
        raise

    # Load dataframe
    csv_kwargs = dict(sep="\t", dtype_backend="pyarrow", low_memory=False)
    try:
        logger.info("Reading MAF with 'pyarrow' engine…")
        maf = pd.read_csv(buf, **csv_kwargs)
        logger.info("Reading with 'pyarrow' completed.")
    except (ValueError, ImportError) as err:
        logger.warning("'pyarrow' failed (%s). Retrying with 'c' engine.", err)
        buf.seek(0)
        maf = pd.read_csv(buf, engine="c", low_memory=False, sep="\t")

    # Validate and standardize columns
    _standardise_maf_columns(maf)
    missing = [c for c in required_columns_MAF if c not in maf.columns]
    if missing:
        msg = f"Missing required columns in MAF: {', '.join(missing)}"
        logger.error(msg)
        raise ValueError(msg)
    logger.debug("Required columns present.")

    # Normalize Variant Classification column names
    maf = normalize_variant_classification(maf)

    # Generate VCF-style fields
    maf["CHROM"] = maf["Chromosome"].astype(str).map(format_chr)
    maf["POS"] = maf["Start_Position"].astype("int64")

    if "dbSNP_RS" in maf.columns:
        maf["ID"] = (
            maf["dbSNP_RS"]
            .fillna(".")
            .replace(".", np.nan)
            .map(lambda x: format_rs(str(x)) if pd.notna(x) else ".")
            .fillna(".")
        )
    else:
        maf["ID"] = "."
    maf["REF"] = maf["Reference_Allele"].astype(str)
    maf["ALT"] = maf["Tumor_Seq_Allele2"].fillna(maf["Tumor_Seq_Allele1"]).astype(str)
    # Check if QUAL and FILTER columns exist before defaulting to "."
    if "QUAL" not in maf.columns:
        maf["QUAL"] = "."

    if "FILTER" not in maf.columns:
        maf["FILTER"] = "."

    # Expand samples to columns
    samples = maf["Tumor_Sample_Barcode"].dropna().unique().tolist()
    logger.info("Detected %d unique samples.", len(samples))

    # Generate a DataFrame with as many columns as samples and with default value
    ref_ref = maf["REF"] + "|" + maf["REF"]
    default_cols = pd.DataFrame(
        {sample: ref_ref for sample in samples},
        index=maf.index,
        dtype="object"  # avoid unnecessary casting
    )

    # Add all columns at once
    maf = pd.concat([maf, default_cols], axis=1)

    ref_alt = maf["REF"] + "|" + maf["ALT"]
    for sample in samples:
        mask = maf["Tumor_Sample_Barcode"] == sample
        maf.loc[mask, sample] = ref_alt[mask]

    # ─── 5.5) CONSOLIDATE DUPLICATE VARIANTS (OPTIONAL) ----------------
    if consolidate_variants:
        logger.info("Consolidating duplicate variants across samples...")
        original_rows = len(maf)
        
        # Usar versión optimizada basada en el tamaño del dataset
        if len(maf) > 100000:  # Para datasets grandes
            maf = _consolidate_variants_chunked(maf, samples)
        else:  # Para datasets normales
            maf = _consolidate_variants(maf, samples)
        
        consolidated_rows = len(maf)
        logger.info("Consolidated %d rows into %d unique variants", 
                   original_rows, consolidated_rows)

    # ─── 5.6) REMOVE UNNECESSARY COLUMNS -------------------------------
    # We no longer need 'Tumor_Sample_Barcode' (mutations were
    # expanded into columns) nor 'Chromosome' (now it's in 'CHROM').
    maf.drop(columns=["Chromosome"], inplace=True)

    # ─── 6) ORDER COLUMNS -------------------------------------------------
    vcf_like = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
    final_cols = vcf_like + samples + [
        c for c in maf.columns if c not in vcf_like + samples
    ]
    maf = maf[final_cols]

    # ─── 7) SAVE TO CACHE ──────────────────────────────────────────────────
    try:
        logger.info("Saving to cache: %s", cache_path)
        maf.to_parquet(cache_path, index=False)

        # Ensure data is written to disk to avoid corruption
        with open(cache_path, 'r+b') as f:
            os.fsync(f.fileno())

        # Save metadata
        cache_info = {
            'comments': comments,
            'samples': samples,
            'processing_time': time.time() - start_time,
            'backend_used': 'pandas'
        }

        import json
        cache_info_path = cache_path.with_suffix('.json')
        with open(cache_info_path, 'w') as f:
            json.dump(cache_info, f, indent=2)
            os.fsync(f.fileno())

        logger.debug("Cache saved successfully with fsync")
    except Exception as e:
        logger.warning("Failed to save cache: %s", e)

    # ─── 8) BUILD PyMutation OBJECT -------------------------------------
    metadata = MutationMetadata(
        source_format="MAF",
        file_path=str(path),
        filters=["."],
        assembly=assembly,
        notes="\n".join(comments) if comments else None,
    )

    total_time = time.time() - start_time
    logger.info("MAF processed successfully: %d rows, %d columns in %.2f seconds",
                *maf.shape, total_time)
    return PyMutation(maf, metadata, samples)


def read_vcf(
        path: str | Path,
        assembly: str,
        create_index: bool = False,
        cache_dir: Optional[str | Path] = None
) -> PyMutation:
    """
    Read a VCF file and return a PyMutation object with optimized performance.

    This function provides a high-performance VCF reader with pandas + PyArrow
    optimization and automatic caching.

    Parameters
    ----------
    path : str | Path
        Path to the VCF file (.vcf or .vcf.gz).
    assembly : str
        Version of the genome assembly. Must be either "37" or "38".
    create_index : bool, default False
        Whether to create a Tabix index (.tbi) if it doesn't exist.
        Requires 'tabix' command to be available in PATH.
    cache_dir : str | Path, optional
        Directory for caching processed files. If None, uses a .pymut_cache
        directory next to the input file.

    Returns
    -------
    PyMutation
        Set of mutations read from the VCF file, converted to the same
        wide-format produced by `read_maf`. Includes metadata with information
        about the source file, comments, and configuration.

    Raises
    ------
    FileNotFoundError
        If the specified VCF file path does not exist.
    ValueError
        If the VCF file is missing required columns or has an invalid format.
        If the assembly parameter is not "37" or "38".
    Exception
        For any other errors encountered while reading or processing the file.

    Notes
    -----
    Performance optimizations include:
    - Vectorized genotype conversion (eliminates pd.apply bottleneck)
    - Direct file reading without StringIO intermediate step
    - PyArrow-accelerated INFO column parsing
    - Automatic caching of processed results
    - Automatic Tabix indexing for compressed files
    """
    # Validate assembly parameter
    if assembly not in ["37", "38"]:
        raise ValueError("Assembly parameter must be either '37' or '38'")

    start_time = time.time()
    path = Path(path)
    cache_dir_path = Path(cache_dir) if cache_dir else None

    logger.info("Starting optimized VCF reading: %s", path)

    # ─── 1) CHECK CACHE ─────────────────────────────────────────────────────
    cache_path = _get_cache_path(path, cache_dir_path)
    if cache_path.exists():
        logger.info("Loading from cache: %s", cache_path)
        try:
            vcf = pd.read_parquet(cache_path)

            # Load metadata from cache info file
            cache_info_path = cache_path.with_suffix('.json')
            if cache_info_path.exists():
                import json
                with open(cache_info_path) as f:
                    cache_info = json.load(f)
                meta_lines = cache_info.get('meta_lines', [])
                sample_columns = cache_info.get('sample_columns', [])
            else:
                meta_lines = []
                sample_columns = []

            metadata = MutationMetadata(
                source_format="VCF",
                file_path=str(path),
                filters=["."],
                assembly=assembly,
                notes="\n".join(meta_lines) if meta_lines else None,
            )

            logger.info("Cache loaded successfully in %.2f seconds", time.time() - start_time)
            return PyMutation(vcf, metadata, sample_columns)

        except Exception as e:
            logger.warning("Cache loading failed (%s), proceeding with fresh read", e)

    # ─── 3) CREATE TABIX INDEX IF NEEDED ───────────────────────────────────
    if path.suffix == '.gz':
        _create_tabix_index(path, create_index)

    # ─── 4) PARSE HEADER AND META-LINES WITHOUT STRINGIO ───────────────────
    meta_lines: list[str] = []
    header_line: str = ""

    try:
        with _open_text_maybe_gzip(path) as fh:
            for line in fh:
                if line.startswith("##"):
                    meta_lines.append(line.rstrip("\n"))
                elif line.startswith("#CHROM"):
                    header_line = line.rstrip("\n")
                    break

        logger.debug("Meta-lines found: %d", len(meta_lines))
        logger.debug("Header line found: %s", "Yes" if header_line else "No")

    except FileNotFoundError:
        logger.error("File not found: %s", path)
        raise
    except Exception:
        logger.exception("Error reading header/meta-lines from VCF.")
        raise

    if not header_line:
        msg = "VCF file is missing the header line starting with #CHROM"
        logger.error(msg)
        raise ValueError(msg)

    # ─── 5) LOAD DATAFRAME WITH PANDAS + PYARROW ────────────────────────────
    header_cols = header_line.lstrip("#").split("\t")

    logger.info("Reading VCF with pandas + pyarrow optimization...")

    # Read directly without StringIO
    csv_kwargs = {
        "sep": "\t",
        "names": header_cols,
        "dtype_backend": "pyarrow" if HAS_PYARROW else None,
        "low_memory": False,
        "skiprows": len(meta_lines) + 1  # Skip meta-lines and header
    }

    try:
        vcf = pd.read_csv(path, **csv_kwargs)
        logger.info("Pandas reading completed.")
    except Exception as err:
        logger.warning("Pandas with pyarrow failed (%s). Retrying with standard engine.", err)
        csv_kwargs.pop("dtype_backend", None)
        vcf = pd.read_csv(path, engine="c", **csv_kwargs)

    # ─── 6) VALIDATE AND STANDARDIZE COLUMNS ───────────────────────────────
    vcf_columns = vcf.columns.tolist()

    # Standardize column names
    rename_map = {
        col: _required_canonical_VCF[col.lstrip("#").lower()]
        for col in vcf_columns
        if col.lstrip("#").lower() in _required_canonical_VCF
    }

    if rename_map:
        vcf = vcf.rename(columns=rename_map)

    missing = [c for c in required_columns_VCF if c not in vcf.columns]
    if missing:
        msg = f"Missing required columns in VCF: {', '.join(missing)}"
        logger.error(msg)
        raise ValueError(msg)
    logger.debug("Required columns present.")

    # First, identify the original sample columns from the header before INFO expansion
    original_header_cols = header_cols
    standard_vcf_cols_temp = {"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"}
    original_sample_columns = [col for col in original_header_cols if col not in standard_vcf_cols_temp]

    # ─── 7) PROCESS CHROMOSOME COLUMN ───────────────────────────────────────
    vcf["CHROM"] = vcf["CHROM"].astype(str).map(format_chr)

    # ─── 7.5) VECTORIZED GENOTYPE CONVERSION ─────────────────────────────────
    logger.info("Starting vectorized genotype conversion before INFO expansion...")

    if original_sample_columns:
        # Vectorized conversion for all sample columns at once
        conversion_start = time.time()

        # Get reference data once
        ref_array = vcf["REF"].astype(str).values
        alt_array = vcf["ALT"].astype(str).values

        # Process all sample columns in batches to manage memory
        batch_size = min(100, len(original_sample_columns))
        for batch_start in range(0, len(original_sample_columns), batch_size):
            batch_end = min(batch_start + batch_size, len(original_sample_columns))
            batch_cols = original_sample_columns[batch_start:batch_end]

            if batch_start % 500 == 0:
                logger.debug("Processing genotype batch %d-%d of %d samples",
                             batch_start + 1, batch_end, len(original_sample_columns))

            # Convert batch of columns
            for sample_col in batch_cols:
                if sample_col in vcf.columns:  # Ensure column exists
                    vcf[sample_col] = _vectorized_gt_to_alleles(
                        vcf[sample_col],
                        vcf["REF"],
                        vcf["ALT"]
                    )

        conversion_time = time.time() - conversion_start
        logger.info("GT conversion completed: %.2f s", conversion_time)

    # ─── 8) EXPAND INFO COLUMN WITH VECTORIZED OPERATIONS ──────────────────
    if "INFO" in vcf.columns:
        logger.info("Expanding INFO column with vectorized operations...")
        expand_start = time.time()

        # Pandas - direct processing
        info_df = _parse_info_column_vectorized(vcf["INFO"])
        vcf = pd.concat([vcf, info_df], axis=1)

        # Remove the original INFO column after expansion
        vcf = vcf.drop(columns=["INFO"])

        logger.debug("INFO column expanded into %d columns in %.2f s",
                     info_df.shape[1], time.time() - expand_start)

    # ─── 9) PROCESS FUNCOTATOR ANNOTATIONS ──────────────────────────────────
    if "FUNCOTATION" in vcf.columns:
        logger.info("Processing Funcotator annotations...")
        try:
            from .utils.data_processing import extract_genome_changes

            funcotator_df = extract_genome_changes(vcf, funcotation_column="FUNCOTATION")
            vcf = pd.concat([vcf, funcotator_df], axis=1)

            logger.debug("Funcotator annotations processed.")
        except ImportError:
            logger.warning("Could not import data_processing module for Funcotator processing.")
        except Exception as e:
            logger.warning("Error processing Funcotator annotations: %s", e)

    # ─── 9.5) CSQ ───────────────────────────────────
    # Check if there's a sample column named "CSQ" in the original header
    has_csq_sample = "CSQ" in original_sample_columns

    csq_sample_new_name = None
    if has_csq_sample and "CSQ" in vcf.columns:
        logger.info("Detected CSQ naming conflict: sample named 'CSQ' conflicts with VEP annotations")

        # The CSQ column now contains INFO-derived data, but we need to separate sample vs INFO CSQ
        csq_sample_new_name = "CSQ_sample"
        counter = 1
        while csq_sample_new_name in vcf.columns:
            csq_sample_new_name = f"CSQ_sample_{counter}"
            counter += 1

        # Create a new column for the sample CSQ data by extracting it from the original VCF
        sample_csq_data = []
        csq_sample_index = original_header_cols.index("CSQ")

        with _open_text_maybe_gzip(path) as fh:
            line_count = 0
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) > csq_sample_index:
                    sample_csq_data.append(parts[csq_sample_index])
                else:
                    sample_csq_data.append(".")
                line_count += 1
                if line_count >= len(vcf):  # Don't read more than we have rows
                    break

        # Add the sample CSQ data as a new column
        vcf[csq_sample_new_name] = sample_csq_data[:len(vcf)]
        logger.info(f"Created sample CSQ column: {csq_sample_new_name}")

    # Expand VEP CSQ annotations
    if "CSQ" in vcf.columns:
        logger.info("Expanding VEP CSQ annotations into individual columns...")
        try:
            csq_expansion_start = time.time()

            # Get the CSQ format from meta-lines to understand the field structure
            csq_format = None
            for meta_line in meta_lines:
                if meta_line.startswith('##INFO=<ID=CSQ') and 'Format:' in meta_line:
                    # Extract format from description
                    format_start = meta_line.find('Format: ') + 8
                    format_end = meta_line.find('"', format_start)
                    if format_end == -1:
                        format_end = meta_line.find('>', format_start)
                    csq_format = meta_line[format_start:format_end]
                    break

            if csq_format:
                csq_fields = [field.strip() for field in csq_format.split('|')]
                logger.debug(f"Found CSQ format with {len(csq_fields)} fields: {csq_fields[:5]}...")

                # Expand CSQ annotations
                csq_expanded_data = []
                for idx, csq_value in enumerate(vcf["CSQ"]):
                    row_data = {}
                    if pd.notna(csq_value) and csq_value and csq_value != ".":
                        # Split multiple CSQ entries (separated by commas)
                        csq_entries = str(csq_value).split(',')

                        if csq_entries:
                            csq_values = csq_entries[0].split('|')
                            for i, field_name in enumerate(csq_fields):
                                if i < len(csq_values):
                                    value = csq_values[i].strip()
                                    row_data[f"VEP_{field_name}"] = value if value else ""
                                else:
                                    row_data[f"VEP_{field_name}"] = ""
                    else:
                        # Fill with empty strings for missing CSQ data
                        for field_name in csq_fields:
                            row_data[f"VEP_{field_name}"] = ""

                    csq_expanded_data.append(row_data)

                csq_expanded_df = pd.DataFrame(csq_expanded_data)

                vcf = pd.concat([vcf, csq_expanded_df], axis=1)

                vcf = vcf.drop(columns=["CSQ"])

                logger.info(f"CSQ expanded into {len(csq_fields)} VEP annotation columns in %.2f s",
                            time.time() - csq_expansion_start)
            else:
                logger.warning("Could not find CSQ format in meta-lines, keeping CSQ as single column")

        except Exception as e:
            logger.warning(f"Error expanding VEP CSQ annotations: {e}")

    # ─── 9.5.5) GENERATE HUGO_SYMBOL FROM VEP_SYMBOL AND VEP_NEAREST ────────
    if 'Hugo_Symbol' in vcf.columns:
        logger.info("Hugo_Symbol column already exists, skipping generation")
    else:
        logger.info("Generating Hugo_Symbol column from VEP_SYMBOL and VEP_NEAREST...")
        try:
            hugo_symbol_start = time.time()

            has_vep_symbol = 'VEP_SYMBOL' in vcf.columns
            has_vep_nearest = 'VEP_NEAREST' in vcf.columns

            if has_vep_symbol or has_vep_nearest:
                # Generate Hugo_Symbol based on VEP_SYMBOL and VEP_NEAREST availability and values
                vcf['Hugo_Symbol'] = ''

                # Get VEP_SYMBOL and VEP_NEAREST columns, handling missing columns
                vep_symbol_series = vcf['VEP_SYMBOL'] if has_vep_symbol else pd.Series([''] * len(vcf), index=vcf.index)
                vep_nearest_series = vcf['VEP_NEAREST'] if has_vep_nearest else pd.Series([''] * len(vcf),
                                                                                          index=vcf.index)

                # Convert None, NaN, or null-like values to empty strings for comparison
                vep_symbol_clean = vep_symbol_series.fillna('').astype(str)
                vep_symbol_clean = vep_symbol_clean.replace(['None', 'null', 'NULL'], '')

                vep_nearest_clean = vep_nearest_series.fillna('').astype(str)
                vep_nearest_clean = vep_nearest_clean.replace(['None', 'null', 'NULL'], '')

                # Apply the logic from the requirements
                if has_vep_symbol and not has_vep_nearest:
                    # Only VEP_SYMBOL exists, use it
                    vcf['Hugo_Symbol'] = vep_symbol_clean
                elif has_vep_symbol and has_vep_nearest:
                    # Both exist, use VEP_SYMBOL unless it's empty, then use VEP_NEAREST
                    vcf['Hugo_Symbol'] = vep_symbol_clean.where(vep_symbol_clean != '', vep_nearest_clean)
                elif not has_vep_symbol and has_vep_nearest:
                    # Only VEP_NEAREST exists, use it
                    vcf['Hugo_Symbol'] = vep_nearest_clean
                else:
                    # Neither exists (this case shouldn't happen due to the if condition above)
                    vcf['Hugo_Symbol'] = ''

                logger.info("Hugo_Symbol column generated in %.2f s",
                            time.time() - hugo_symbol_start)
                logger.debug(
                    f"Hugo_Symbol: VEP_SYMBOL available: {has_vep_symbol}, VEP_NEAREST available: {has_vep_nearest}")

                # Log some statistics
                non_empty_count = (vcf['Hugo_Symbol'] != '').sum()
                logger.debug(f"Hugo_Symbol: {non_empty_count} non-empty values out of {len(vcf)} total rows")
            else:
                logger.debug("Neither VEP_SYMBOL nor VEP_NEAREST columns found, skipping Hugo_Symbol generation")

        except Exception as e:
            logger.warning(f"Error generating Hugo_Symbol: {e}")

    # ─── 9.6) GENERATE VARIANT_CLASSIFICATION FROM VEP DATA ─────────────────
    if 'Variant_Classification' in vcf.columns:
        logger.info("Variant_Classification column already exists, skipping generation")
    elif "VEP_Consequence" in vcf.columns and "VEP_VARIANT_CLASS" in vcf.columns:
        logger.info("Generating Variant_Classification from VEP_Consequence and VEP_VARIANT_CLASS...")
        try:
            variant_classification_start = time.time()

            # Get primary consequence using vectorized operations
            vcf['prim_cons'] = _get_primary_consequence(vcf['VEP_Consequence'])

            # Apply classification function
            vcf['Variant_Classification'] = vcf.apply(_classify_variant, axis=1)

            # Fill any NaN values with 'RNA' (for non-coding transcripts)
            vcf['Variant_Classification'] = vcf['Variant_Classification'].fillna('RNA')

            # Convert to categorical for memory efficiency and ordering
            # First, let's see what unique values we actually have
            unique_values = vcf['Variant_Classification'].unique()
            logger.debug(f"Unique Variant_Classification values before categorical conversion: {unique_values}")

            ordered_categories = [
                'Nonsense_Mutation', 'Nonstop_Mutation',
                'Splice_Site', 'Splice_Region',
                'Frame_Shift_Del', 'Frame_Shift_Ins',
                'In_Frame_Del', 'In_Frame_Ins',
                'Missense_Mutation', 'Silent',
                "5'UTR", "3'UTR", "5'Flank", "3'Flank",
                'Intron', 'IGR', 'RNA'
            ]

            # Add any unique values that aren't in our predefined categories
            for val in unique_values:
                if val not in ordered_categories and pd.notna(val):
                    ordered_categories.append(val)

            vcf['Variant_Classification'] = pd.Categorical(
                vcf['Variant_Classification'],
                categories=ordered_categories,
                ordered=True
            )

            # Clean up temporary column
            vcf = vcf.drop(columns=['prim_cons'])

            logger.info("Variant_Classification generated in %.2f s",
                        time.time() - variant_classification_start)

        except Exception as e:
            logger.warning(f"Error generating Variant_Classification: {e}")
    else:
        logger.debug(
            "VEP_Consequence or VEP_VARIANT_CLASS columns not found, skipping Variant_Classification generation")

    # ─── 9.7) NORMALIZE VARIANT_CLASSIFICATION TO UPPERCASE ─────────────────
    vcf = normalize_variant_classification(vcf)
    logger.debug("Variant_Classification values normalized to uppercase")

    # ─── 9.8) GENERATE VARIANT_TYPE FROM VEP VARIANT_CLASS ──────────────────
    if 'Variant_Type' in vcf.columns:
        logger.info("Variant_Type column already exists, skipping generation")
    elif "VEP_VARIANT_CLASS" in vcf.columns:
        logger.info("Generating Variant_Type from VEP_VARIANT_CLASS...")
        try:
            variant_type_start = time.time()

            # Apply variant type mapping function
            vcf['Variant_Type'] = vcf.apply(_map_vep_variant_class_to_maf_variant_type, axis=1)

            # Fill any NaN values with 'UNKNOWN'
            vcf['Variant_Type'] = vcf['Variant_Type'].fillna('UNKNOWN')

            # Convert to categorical for memory efficiency and ordering
            unique_variant_types = vcf['Variant_Type'].unique()
            logger.debug(f"Unique Variant_Type values: {unique_variant_types}")

            # Define ordered categories for Variant_Type
            variant_type_categories = [
                'SNP', 'DNP', 'TNP', 'ONP',  # Substitutions
                'INS', 'DEL', 'IND',  # Indels
                'INV', 'CNV',  # Structural variants
                'UNKNOWN'  # Fallback
            ]

            # Add any unique values that aren't in our predefined categories
            for val in unique_variant_types:
                if val not in variant_type_categories and pd.notna(val):
                    variant_type_categories.append(val)

            vcf['Variant_Type'] = pd.Categorical(
                vcf['Variant_Type'],
                categories=variant_type_categories,
                ordered=True
            )

            logger.info("Variant_Type generated in %.2f s",
                        time.time() - variant_type_start)

        except Exception as e:
            logger.warning(f"Error generating Variant_Type: {e}")
    else:
        logger.debug("VEP_VARIANT_CLASS column not found, skipping Variant_Type generation")

    # ─── 10) DEFINE SAMPLE COLUMNS FOR LATER USE ─────────────────────────────
    # Use original_sample_columns but filter out any columns that might have been added during processing
    standard_vcf_cols = {"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"}
    all_columns = vcf.columns.tolist()
    sample_columns = [col for col in original_sample_columns if col in all_columns]

    logger.info("Using %d sample columns for output organization", len(sample_columns))

    # ─── 11) ENSURE REQUIRED COLUMNS EXIST ──────────────────────────────────
    if "QUAL" not in vcf.columns:
        vcf["QUAL"] = "."
    if "FILTER" not in vcf.columns:
        vcf["FILTER"] = "."

    # Remove FORMAT column as it's not needed
    if "FORMAT" in vcf.columns:
        vcf = vcf.drop(columns=["FORMAT"])
        logger.debug("FORMAT column removed from output DataFrame")

    # ─── 12) ORDER COLUMNS ──────────────────────────────────────────────────
    vcf_core = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
    final_cols = vcf_core + sample_columns + [
        c for c in vcf.columns if c not in vcf_core + sample_columns
    ]
    vcf = vcf[final_cols]

    # ─── 13) SAVE TO CACHE ──────────────────────────────────────────────────
    try:
        logger.info("Saving to cache: %s", cache_path)
        vcf.to_parquet(cache_path, index=False)

        # Ensure data is written to disk to avoid corruption
        with open(cache_path, 'r+b') as f:
            os.fsync(f.fileno())

        # Save metadata
        cache_info = {
            'meta_lines': meta_lines,
            'sample_columns': sample_columns,
            'processing_time': time.time() - start_time,
            'backend_used': 'pandas'
        }

        import json
        cache_info_path = cache_path.with_suffix('.json')
        with open(cache_info_path, 'w') as f:
            json.dump(cache_info, f, indent=2)
            os.fsync(f.fileno())

        logger.debug("Cache saved successfully with fsync")
    except Exception as e:
        logger.warning("Failed to save cache: %s", e)

    # ─── 14) BUILD PyMutation OBJECT ────────────────────────────────────────
    metadata = MutationMetadata(
        source_format="VCF",
        file_path=str(path),
        filters=["."],
        assembly=assembly,
        notes="\n".join(meta_lines) if meta_lines else None,
    )

    total_time = time.time() - start_time
    logger.info("VCF processed successfully: %d rows, %d columns in %.2f seconds",
                *vcf.shape, total_time)

    return PyMutation(vcf, metadata, sample_columns)
