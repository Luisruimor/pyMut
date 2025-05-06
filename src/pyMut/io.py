import gzip
import logging
import re
from pathlib import Path
from typing import List

import pandas as pd

from .core import PyMutation, MutationMetadata

# ────────────────────────────────────────────────────────────────
# LOGGING DEL MÓDULO
# ────────────────────────────────────────────────────────────────
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # Cambia a DEBUG para más verbosidad
if not logger.handlers:
    logging.basicConfig(
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        level=logging.INFO,
    )

# ────────────────────────────────────────────────────────────────
# VALIDACIONES DE COLUMNAS
# ────────────────────────────────────────────────────────────────
# Columnas requeridas para VCF (según el estándar)
required_columns_VCF: List[str] = ["CHROM", "POS", "ID", "REF", "ALT", "FILTER"]
_required_canonical_VCF = {c.lower(): c for c in required_columns_VCF}

# Columnas requeridas para MAF
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
# UTILIDADES INTERNAS
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
    """Normaliza nombres de columnas (insensible a mayúsculas) para MAF."""
    rename_map = {
        col: _required_canonical_MAF[col.lower()]
        for col in maf.columns
        if col.lower() in _required_canonical_MAF
    }
    maf.rename(columns=rename_map, inplace=True)
    return maf


def _open_text_maybe_gzip(path: str | Path):
    """Devuelve un manejador de texto, aceptando .gz o texto plano."""
    path = Path(path)
    if path.suffix == ".gz":
        logger.debug("Abriendo %s como gzip-text.", path)
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    logger.debug("Abriendo %s como texto plano.", path)
    return open(path, encoding="utf-8", errors="replace")


def _parse_info_column(info_series: pd.Series) -> pd.DataFrame:
    """
    Expande la columna INFO (clave-valor separados por ';') a un DataFrame.

    Flags sin valor se marcan con True. El resultado tiene las mismas filas
    que `info_series` y las claves como columnas.
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
                # Flag sin valor (por ejemplo, "PASS")
                out[entry] = True
        return out

    parsed = info_series.apply(_parser)
    info_df = pd.DataFrame(parsed.tolist())  # type: ignore[arg-type]
    return info_df


# ════════════════════════════════════════════════════════════════
# FUNCIÓN PRINCIPAL: read_vcf
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