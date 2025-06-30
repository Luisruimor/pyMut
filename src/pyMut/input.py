import gzip
import logging
import re
from typing import List
import io
from pathlib import Path
import numpy as np
import pandas as pd

from .core import PyMutation, MutationMetadata
from .utils.format import formatear_rs, formatear_chr, normalize_variant_classification

# ────────────────────────────────────────────────────────────────
# LOGGER CONFIGURATION
# ────────────────────────────────────────────────────────────────
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # Change to DEBUG for more verbosity
if not logger.handlers:
    logging.basicConfig(
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        level=logging.INFO,
    )

# ────────────────────────────────────────────────────────────────
# COLUMNS VALIDATION
# ────────────────────────────────────────────────────────────────
required_columns_VCF: List[str] = ["CHROM", "POS", "ID", "REF", "ALT", "FILTER"]
_required_canonical_VCF = {c.lower(): c for c in required_columns_VCF}

required_columns_MAF: List[str] = [
    "Chromosome",
    "Start_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode",
]
_required_canonical_MAF = {c.lower(): c for c in required_columns_MAF}


# ════════════════════════════════════════════════════════════════
# FUNCTIONS FOR MAF PROCESSING
# ════════════════════════════════════════════════════════════════
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
    if not gt:           # empty string
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

# ════════════════════════════════════════════════════════════════
# MAIN FUNCTIONS
# ════════════════════════════════════════════════════════════════
def read_vcf(path: str | Path) -> PyMutation:
    """
    Lee un archivo .vcf o .vcf.gz y devuelve un objeto `PyMutation`.
    """

    path = Path(path)
    logger.info("Iniciando lectura de VCF: %s", path)

    # ─── 1) CONTAR METALÍNEAS Y DETECTAR FUNCOTATION ─────────────────────────
    meta_lines = 0
    funcotator_fields: List[str] | None = None
    try:
        with _open_text_maybe_gzip(path) as fh:
            for ln in fh:
                if ln.startswith("##"):
                    meta_lines += 1
                    # ¿Contiene definiciones de Funcotator?
                    if ln.startswith("##INFO=<ID=FUNCOTATION"):
                        # Regex para extraer campos listados por Funcotator
                        match = re.search(r"Funcotation fields are: ([^>]+)", ln)
                        if match:
                            funcotator_fields = match.group(1).split("|")
                elif ln.startswith("#"):
                    break
                else:
                    # Ni es ## ni es # → fin de meta-líneas, comienzo de datos
                    break
        logger.debug(
            "Meta-líneas: %s | Campos FUNCOTATOR: %s",
            meta_lines,
            funcotator_fields,
        )
    except FileNotFoundError as exc:
        logger.error("Archivo no encontrado: %s", path)
        raise
    except Exception as exc:
        logger.exception("Error al leer cabecera/metalíneas.")
        raise

    # ─── 2) CARGAR DATAFRAME ────────────────────────────────────────────────
    csv_kwargs = dict(
        sep="\t",
        dtype=str,
        skiprows=meta_lines,
        header=0,
        compression="infer",
    )

    try:
        logger.info("Leyendo VCF con motor 'pyarrow'…")
        vcf_df = pd.read_csv(path, engine="pyarrow", **csv_kwargs)
        logger.info("Lectura con 'pyarrow' completada.")
    except (ValueError, ImportError) as err:
        logger.warning("Falló 'pyarrow' (%s).  Reintentando con motor 'c'.", err)
        vcf_df = pd.read_csv(path, engine="c", low_memory=False, **csv_kwargs)

    # Eliminamos el # de la columna "#CHROM"
    first_col = vcf_df.columns[0]
    if first_col.startswith("#"):
        vcf_df.rename(columns={first_col: first_col.lstrip("#")}, inplace=True)

    # ─── 3) VALIDAR COLUMNAS ────────────────────────────────────────────────
    vcf_df = _standardise_vcf_columns(vcf_df)
    missing = [c for c in required_columns_VCF if c not in vcf_df.columns]
    if missing:
        msg = f"Faltan columnas requeridas en el VCF: {', '.join(missing)}"
        logger.error(msg)
        raise ValueError(msg)

    if "FORMAT" in vcf_df.columns:
        vcf_df.drop(columns="FORMAT", inplace=True)

    # Convertir genotipos en columnas de muestras a alelos reales
    standard_cols = set(required_columns_VCF + ["INFO", "QUAL"])
    sample_cols = [col for col in vcf_df.columns if col not in standard_cols]

    if sample_cols:
        logger.info("Convirtiendo genotipos a alelos para %d muestras...", len(sample_cols))
        for col in sample_cols:
            vcf_df[col] = vcf_df.apply(
                lambda row: _gt_to_alleles(row[col], row["REF"], row["ALT"]),
                axis=1
            )


    # ─── 4) EXPANSIÓN DE LA COLUMNA INFO ────────────────────────────────────
    if "INFO" in vcf_df.columns:
        logger.info("Expandiendo columna INFO a campos individuales…")
        info_expanded = _parse_info_column(vcf_df["INFO"])
        # Evitar colisiones de nombres
        cols_duplicadas = set(vcf_df.columns) & set(info_expanded.columns)
        if cols_duplicadas:
            info_expanded.rename(
                columns={c: f"INFO_{c}" for c in cols_duplicadas}, inplace=True
            )
        vcf_df = pd.concat([vcf_df.drop(columns=["INFO"]), info_expanded], axis=1)
        logger.debug("INFO expandido; columnas añadidas: %s", list(info_expanded.columns))

    # ─── 5) EXPANSIÓN FUNCOTATOR ────────────────────────────────────────────
    if funcotator_fields and "FUNCOTATION" in vcf_df.columns:
        logger.info("Expandiendo anotaciones FUNCOTATOR…")
        # Separar múltiples anotaciones por coma
        raw_func = vcf_df["FUNCOTATION"].fillna("")
        exploded = raw_func.str.split(",", expand=False)

        # Repetir cada fila tantas veces como anotaciones tenga
        counts = exploded.str.len()
        vcf_rep = vcf_df.loc[vcf_df.index.repeat(counts)].reset_index(drop=True)

        func_series = exploded.explode().reset_index(drop=True).str.strip("[]")
        func_df = func_series.str.split("|", expand=True)
        func_df.columns = funcotator_fields

        vcf_df = pd.concat([vcf_rep.drop(columns=["FUNCOTATION"]), func_df], axis=1)
        logger.debug("FUNCOTATOR expandido en %d columnas.", len(funcotator_fields))

    # ─── 6) POST-PROCESADO FINAL ────────────────────────────────────────────
    vcf_df["CHROM"] = vcf_df["CHROM"].replace(
        {"23": "X", "24": "Y", "chr23": "X", "chr24": "Y"}
    )

    # ─── 7) CREAR OBJETO PyMutation ────────────────────────────────────────
    metadata = MutationMetadata(
        source_format="VCF",
        file_path=str(path),
        filters=["None"],
        fasta="None",
        notes=(
            "Archivo VCF cargado correctamente."
            + (
                " Se expandieron anotaciones Funcotator."
                if funcotator_fields
                else ""
            )
        ),
    )
    logger.info("VCF procesado con éxito: %d filas, %d columnas.", *vcf_df.shape)
    return PyMutation(data=vcf_df, metadata=metadata)

def read_maf(path: str | Path, fasta: str | Path | None = None) -> PyMutation:
    """
    Read a MAF file and return a PyMutation object.

    Parameters
    ----------
    path : str | Path
        Path to the MAF file (.maf or .maf.gz).
    fasta : str | Path, optional
        Path to the reference FASTA file. If provided, it will be included
        in the metadata of the resulting PyMutation object.

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
    ImportError
        If the 'pyarrow' library cannot be imported and the 'c' engine
        alternative also fails.
    Exception
        For any other errors encountered while reading or processing the file.
    """
    path = Path(path)
    logger.info("Starting MAF reading: %s", path)

    # ─── 1) SEPARATE COMMENTS (#) FROM BODY ──────────────────────────────
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

    # ─── 2) LOAD DATAFRAME -------------------------------------------------
    csv_kwargs = dict(sep="\t", dtype_backend="pyarrow", low_memory=False)
    try:
        logger.info("Reading MAF with 'pyarrow' engine…")
        maf = pd.read_csv(buf, **csv_kwargs)
        logger.info("Reading with 'pyarrow' completed.")
    except (ValueError, ImportError) as err:
        logger.warning("'pyarrow' failed (%s). Retrying with 'c' engine.", err)
        buf.seek(0)
        maf = pd.read_csv(buf, engine="c", low_memory=False, sep="\t")

    # ─── 3) VALIDATE AND STANDARDIZE COLUMNS ---------------------------------
    _standardise_maf_columns(maf)
    missing = [c for c in required_columns_MAF if c not in maf.columns]
    if missing:
        msg = f"Missing required columns in MAF: {', '.join(missing)}"
        logger.error(msg)
        raise ValueError(msg)
    logger.debug("Required columns present.")

    # Normalize Variant Classification column names
    maf = normalize_variant_classification(maf)

    # ─── 4) GENERATE VCF-STYLE FIELDS ---------------------------------------
    maf["CHROM"] = maf["Chromosome"].astype(str).map(formatear_chr)
    maf["POS"] = maf["Start_Position"].astype("int64")

    if "dbSNP_RS" in maf.columns:
        maf["ID"] = (
            maf["dbSNP_RS"]
            .fillna(".")
            .replace(".", np.nan)
            .map(lambda x: formatear_rs(str(x)) if pd.notna(x) else ".")
            .fillna(".")
        )
    else:
        maf["ID"] = "."
    maf["REF"] = maf["Reference_Allele"].astype(str)
    maf["ALT"] = maf["Tumor_Seq_Allele2"].fillna(maf["Tumor_Seq_Allele1"]).astype(str)
    maf["QUAL"] = "."
    maf["FILTER"] = "."

    # ─── 5) EXPAND SAMPLES TO COLUMNS ------------------------------------
    samples = maf["Tumor_Sample_Barcode"].dropna().unique().tolist()
    logger.info("Detected %d unique samples.", len(samples))

    # Generate a DataFrame with as many columns as samples and with default value
    ref_ref = maf["REF"] + "|" + maf["REF"]
    default_cols = pd.DataFrame(
        {sample: ref_ref for sample in samples},
        index=maf.index,
        dtype="object"            # avoid unnecessary casting
    )

    # Add all columns at once
    maf = pd.concat([maf, default_cols], axis=1)

    ref_alt = maf["REF"] + "|" + maf["ALT"]
    for sample in samples:
        mask = maf["Tumor_Sample_Barcode"] == sample
        maf.loc[mask, sample] = ref_alt[mask]

    # ─── 5.5) REMOVE UNNECESSARY COLUMNS -------------------------------
    # We no longer need 'Tumor_Sample_Barcode' (mutations were
    # expanded into columns) nor 'Chromosome' (now it's in 'CHROM').
    maf.drop(columns=["Chromosome"], inplace=True)

    # ─── 6) ORDER COLUMNS -------------------------------------------------
    vcf_like = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
    final_cols = vcf_like + samples + [
        c for c in maf.columns if c not in vcf_like + samples
    ]
    maf = maf[final_cols]

    # ─── 7) BUILD PyMutation OBJECT -------------------------------------
    metadata = MutationMetadata(
        source_format="MAF",
        file_path=str(path),
        filters=["."],
        fasta=str(fasta) if fasta else "",
        notes="\n".join(comments) if comments else None,
    )

    logger.info("MAF processed successfully: %d rows, %d columns.", *maf.shape)
    return PyMutation(maf, metadata,samples)