"""
Field mapping utilities for pyMut.

Maps canonical field names to their aliases in MAF/VEP/ANNOVAR/cBioPortal outputs.
Provides utilities to get canonical names from any alias without modifying the original DataFrame.

Example:
>>> from fields import col
>>> aa_series = col(df, "Protein_position")
>>> uniprot_series = col(df, "UNIPROT", required=True)
"""
from __future__ import annotations

from typing import Dict, List, Sequence

# Main dictionary: key = canonical name, value = list of aliases
FIELDS: Dict[str, List[str]] = {
    # Posición de aminoácido en la proteína
    "Protein_position": [
        "Protein_position",
        "VEP_Protein_position",
        "AAPos",
        "AA_position",
        "aa_pos",
        "AA_pos",
        "AaPos",
    ],
    # Símbolo de gen (HGNC)
    "Hugo_Symbol": [
        "Hugo_Symbol",
        "HUGO_SYMBOL",
        "VEP_SYMBOL",
        "SYMBOL",
        "VEP_NEAREST",
        "NEAREST",
        "Gene_Symbol",
        "GENE_SYMBOL",
    ],
    # Consecuencia funcional de la variante
    "Consequence": [
        "VEP_Consequence",
        "Consequence",
        "Func.refGene",
        "Effect",
    ],
    # Identificador de gen Ensembl
    "Gene_ID": [
        "VEP_Gene",
        "Gene_ID",
        "ENSG_ID",
    ],
    # Identificador de transcrito Ensembl
    "Transcript_ID": [
        "VEP_Feature",
        "Transcript_ID",
        "ENST",
        "Transcript",
    ],
    # Tipo de feature (Transcript, miRNA, etc.)
    "Feature_type": [
        "VEP_Feature_type",
        "Feature_type",
        "FeatureType",
    ],
    # Posición en cDNA
    "cDNA_position": [
        "VEP_cDNA_position",
        "cDNA_position",
        "CDNAPos",
        "HGVSc",
    ],
    # Posición en CDS
    "CDS_position": [
        "VEP_CDS_position",
        "CDS_position",
    ],
    # Impacto (HIGH / MODERATE / …)
    "IMPACT": [
        "VEP_IMPACT",
        "IMPACT",
        "Severity",
    ],
    # UniProt
    "UNIPROT": [
        "VEP_SWISSPROT",
        "VEP_TREMBL",
        "UNIPROT",
        "SWISSPROT",
        "TREMBL",
        "UNIPARC",
        "UNIPROT_ISOFORM",
        "uniprot",
        "Uniprot",
    ],
    # Dominios PFAM/proteína
    "Domains": [
        "VEP_DOMAINS",
        "DOMAINS",
        "Protein_domains",
        "Pfam_domain",
    ],
    # Hebra / orientación
    "Strand": [
        "VEP_STRAND",
        "STRAND",
        "Orientation",
    ],
    # Distancia a gen
    "Distance": [
        "VEP_DISTANCE",
        "DISTANCE",
        "Distance_to_gene",
    ],
    # Cambio de aminoácido (HGVSp)
    "Protein_Change": [
        "Protein_Change",
        "HGVSp",
        "HGVSp_Short",
        "Hgvsp",
        "AAChange",
        "ProteinChange",
    ],
    # Clasificación de variante
    "Variant_Classification": [
        "Variant_Classification",
        "VariantClass",
        "Classification",
        "SeqOntology",
        "SO_term",
    ],
    # Identificador PFAM
    "Pfam_ID": [
        "pfam_id",
        "PFAM_ID",
        "PfamID",
    ],
    # Nombre PFAM
    "Pfam_Name": [
        "pfam_name",
        "PFAM_NAME",
        "PfamName",
        "PFAM_Name",
    ],
    # Tipo de variante
    "Variant_Type": [
        "Variant_Type",
        "VariantType",
        "Mutation_Type",
        "MutationType",
    ],
    # Cromosoma (CONSOLIDADO)
    "Chromosome": [
        "Chromosome",
        "Chr",
        "CHROM",
        "CHR",
        "chromosome",
        "chr",
        "chrom",
    ],
    # Posición de inicio (CONSOLIDADO)
    "Start_Position": [
        "Start_Position",
        "Start_position",
        "StartPosition",
        "POS",
        "Position",
        "Pos",
        "pos",
        "start_pos",
        "start_position",
        "start",
        "Start",
    ],
    # Alelo de referencia (CONSOLIDADO)
    "Reference_Allele": [
        "Reference_Allele",
        "Ref_Allele",
        "REF",
        "Reference",
        "Ref",
        "reference_allele",
        "ref",
        "reference",
    ],
    # Alelo tumoral 1
    "Tumor_Seq_Allele1": [
        "Tumor_Seq_Allele1",
        "Tumor_Allele1",
        "TumorAllele1",
        "Allele1",
        "tumor_seq_allele1",
    ],
    # Alelo tumoral 2 (CONSOLIDADO)
    "Tumor_Seq_Allele2": [
        "Tumor_Seq_Allele2",
        "Tumor_Allele2",
        "TumorAllele2",
        "Allele2",
        "ALT",
        "Alt",
        "alt",
        "Alternative",
        "alternative",
        "Tumor_Allele",
        "tumor_seq_allele2",
    ],
    # Código de barras de muestra tumoral (CONSOLIDADO)
    "Tumor_Sample_Barcode": [
        "Tumor_Sample_Barcode",
        "TumorSampleBarcode",
        "Tumor_Sample_ID",
        "Sample_ID",
        "SampleID",
        "SAMPLE",
        "Sample",
        "sample",
        "sample_id",
        "Barcode",
        "barcode",
        "tumor_sample_barcode",
    ],
}

# Reverse alias mapping: alias → canonical name
ALIAS_TO_CANONICAL: Dict[str, str] = {
    alias: canonical for canonical, aliases in FIELDS.items() for alias in aliases
}


# Utility functions

def canonical_name(column: str) -> str:
    """Devuelve el nombre canónico asociado a *column* (alias)."""
    return ALIAS_TO_CANONICAL.get(column, column)


def find_alias(columns: Sequence[str], canonical: str) -> str | None:
    """Devuelve el *primer* alias de *canonical* presente en *columns*.

    Si no se encuentra ninguno devuelve ``None``.
    """
    for alias in FIELDS.get(canonical, []):
        if alias in columns:
            return alias
    return None


def col(df, canonical: str, *, required: bool = False):
    """Accede a la serie DataFrame correspondiente al campo *canonical*.

    Parameters
    ----------
    df : pandas.DataFrame
    canonical : str
        Nombre canónico definido en ``FIELDS``.
    required : bool, default ``False``
        Si ``True`` provoca ``KeyError`` cuando no se encuentra
        ningún alias. Si ``False`` devuelve ``None``.
    """
    alias = find_alias(df.columns, canonical)
    if alias:
        return df[alias]
    if required:
        raise KeyError(f"Ningún alias de '{canonical}' encontrado en DataFrame")
    return None


def canonicalize_columns(columns: List[str]) -> List[str]:
    """Conservado para compatibilidad, pero normalmente **no** se usa
    si prefieres dejar las columnas intactas y llamar a ``col``."""
    return [canonical_name(c) for c in columns]


# Para añadir un nuevo campo:
# 1. Edita este archivo
# 2. Añade tu clave canónica y su lista de alias, p. ej.:
#    FIELDS["ClinSig"] = ["ClinSig", "Clinical_Significance", "CLINSIG"]
# No olvides reconstruir ALIAS_TO_CANONICAL si lo haces de forma dinámica.

__all__ = [
    "FIELDS",
    "ALIAS_TO_CANONICAL",
    "canonical_name",
    "find_alias",
    "col",
    "canonicalize_columns",
]
