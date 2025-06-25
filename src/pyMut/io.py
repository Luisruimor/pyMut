import gzip
import logging
import re
from typing import List
import io
from pathlib import Path
import numpy as np
import pandas as pd

from .core import PyMutation, MutationMetadata
from .utils.format import formatear_rs, formatear_chr

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

def normalizar_variant_classification(df: pd.DataFrame) -> pd.DataFrame:
    """
    Converts values in any 'Variant Classification' column to uppercase,
    regardless of Gencode prefix, version, or capitalization.

    Examples of matched column names:
        - Gencode_43_variantClassification
        - gencode_34_variantclassification
        - variant_classification
        - Variant_Classification
        - gencode_99_VariantClassification

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.

    Returns
    -------
    pd.DataFrame
        The same DataFrame with matched columns normalized to uppercase.
        The same reference is returned to allow method chaining if desired.
    """
    # Regular expression:
    #  - ^                  : start of string
    #  - (gencode_\d+_)?    : optional prefix 'gencode_<num>_' (case insensitive)
    #  - variant[_]?classification : body of the name (allows 'variantclassification' or with '_')
    #  - $                  : end of string
    patron = re.compile(r'^(gencode_\d+_)?variant[_]?classification$', flags=re.IGNORECASE)

    # Find columns that match the pattern
    columnas_objetivo = [col for col in df.columns if patron.match(col)]

    # Convert values to uppercase for each found column
    for col in columnas_objetivo:
        # Only makes sense for object type columns (strings)
        if pd.api.types.is_string_dtype(df[col]):
            df[col] = df[col].str.upper()

    return df

# ════════════════════════════════════════════════════════════════
# MAIN FUNCTIONS
# ════════════════════════════════════════════════════════════════
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
    maf = normalizar_variant_classification(maf)

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
    maf.drop(columns=["Tumor_Sample_Barcode", "Chromosome"], inplace=True)

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