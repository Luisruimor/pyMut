from typing import List, Optional
import pandas as pd
from datetime import datetime

class MutationMetadata:
    """
    Clase para almacenar metadatos de mutaciones.

    Atributos:
        source_format (str): Formato de origen (VCF, MAF, etc.).
        file_path (str): Ruta del archivo de origen.
        loaded_at (datetime): Fecha y hora de carga.
        filters (List[str]): Filtros aplicados al archivo.
        fasta (str): Ruta del archivo FASTA.
        notes (Optional[str]): Notas adicionales.
    """
    def __init__(self, source_format: str, file_path: str,
                 filters: List[str],fasta: str, notes: Optional[str] = None):
        self.source_format = source_format
        self.file_path = file_path
        self.loaded_at = datetime.now()
        self.filters = filters
        self.notes = notes
        self.fasta = fasta


class PyMutation:
    def __init__(self, data: pd.DataFrame, metadata: MutationMetadata):
        self.data = data
        self.metadata = metadata