import logging
import pandas as pd
from pathlib import Path
from typing import Union

from ..utils.fields import col

# Configure logger
logger = logging.getLogger(__name__)

def export_mutation_to_tsv(data: pd.DataFrame, output_path: Union[str, Path], assembly: str) -> None:
    """
    Export mutation data to a TSV file with specific columns.
    
    This function extracts the following columns from the DataFrame and renames them:
    - Chromosome -> CHROM
    - Start_Position -> POS
    - End_Position -> END (calculated as POS + len(REF))
    - Reference_Allele -> REF
    - Tumor_Seq_Allele2 -> ALT
    
    Parameters
    ----------
    data : pd.DataFrame
        DataFrame containing mutation data.
    output_path : str | Path
        Path where the TSV file will be written.
    assembly : str
        The genome assembly to use (required).
        
    Raises
    ------
    ValueError
        If the DataFrame doesn't contain the necessary data for export.
    """
    output_path = Path(output_path)
    logger.info("Starting TSV export to: %s", output_path)
    
    # Extract required columns using the col function to handle aliases
    export_data = pd.DataFrame()
    
    # Extract Chromosome -> CHROM
    chromosome = col(data, "Chromosome", required=False)
    if chromosome is None:
        raise ValueError("Missing required column 'Chromosome' for TSV export")
    export_data["CHROM"] = chromosome
    
    # Extract Start_Position -> POS
    start_position = col(data, "Start_Position", required=False)
    if start_position is None:
        raise ValueError("Missing required column 'Start_Position' for TSV export")
    export_data["POS"] = start_position
    
    # Extract Reference_Allele -> REF
    reference_allele = col(data, "Reference_Allele", required=False)
    if reference_allele is None:
        raise ValueError("Missing required column 'Reference_Allele' for TSV export")
    export_data["REF"] = reference_allele
    
    # Calculate End_Position as POS + len(REF)
    export_data["END"] = export_data["POS"] + export_data["REF"].str.len()
    
    # Extract Tumor_Seq_Allele2 -> ALT
    tumor_seq_allele2 = col(data, "Tumor_Seq_Allele2", required=False)
    if tumor_seq_allele2 is None:
        # Try to use ALT column if Tumor_Seq_Allele2 is not present
        alt = col(data, "ALT", required=False)
        if alt is not None:
            logger.info("Using ALT column for Tumor_Seq_Allele2")
            export_data["ALT"] = alt
        else:
            raise ValueError("Missing required column 'Tumor_Seq_Allele2' for TSV export")
    else:
        export_data["ALT"] = tumor_seq_allele2
    
    # Generate referenceGenome variable
    referenceGenome = f"GRCh{assembly}"
    logger.info("Using reference genome: %s", referenceGenome)
    
    logger.info("Dataframe completed successfully: %d rows written to %s", len(export_data))