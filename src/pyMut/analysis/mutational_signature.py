"""
Mutational signature analysis module for PyMutation.

This module provides functionality for analyzing trinucleotide contexts
and generating mutational signature matrices.
"""

import logging
from typing import Tuple, Dict, Optional

import numpy as np
import pandas as pd

# Set up logging
logger = logging.getLogger(__name__)

# Define the 96 trinucleotide contexts in standard order
TRINUCLEOTIDE_CONTEXTS = [
    # C>A mutations
    "A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T",
    "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T",
    "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T",
    "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T",

    # C>G mutations
    "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T",
    "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T",
    "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T",
    "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T",

    # C>T mutations
    "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T",
    "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T",
    "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T",
    "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",

    # T>A mutations
    "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T",
    "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T",
    "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T",
    "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T",

    # T>C mutations
    "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T",
    "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T",
    "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T",
    "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T",

    # T>G mutations
    "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T",
    "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T",
    "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T",
    "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T"
]

CONTEXT_TO_INDEX = {context: idx for idx, context in enumerate(TRINUCLEOTIDE_CONTEXTS)}
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def _get_reverse_complement(sequence: str) -> str:
    """Get reverse complement of a DNA sequence."""
    return ''.join(COMPLEMENT[base] for base in reversed(sequence))


def _normalize_to_pyrimidine(ref: str, alt: str, trinuc: str) -> Tuple[str, str, str]:
    """
    Normalize mutation to pyrimidine context (C or T as reference).

    Parameters
    ----------
    ref : str
        Reference allele
    alt : str
        Alternative allele
    trinuc : str
        Trinucleotide context

    Returns
    -------
    Tuple[str, str, str]
        Normalized (ref, alt, trinuc)
    """
    if ref in ['C', 'T']:
        return ref, alt, trinuc
    else:
        # Convert to reverse complement for pyrimidine context
        new_ref = COMPLEMENT[ref]
        new_alt = COMPLEMENT[alt]
        new_trinuc = _get_reverse_complement(trinuc)
        return new_ref, new_alt, new_trinuc


def _get_trinucleotide_context(fasta, chrom: str, pos: int) -> Optional[str]:
    """
    Extract trinucleotide context from FASTA file.

    Parameters
    ----------
    fasta : pyfaidx.Fasta
        FASTA file object
    chrom : str
        Chromosome name
    pos : int
        Position (1-based)

    Returns
    -------
    Optional[str]
        Trinucleotide context or None if not available
    """
    try:
        # Handle chromosome name variations
        chrom_key = chrom
        if chrom_key not in fasta.keys():
            if not chrom_key.startswith('chr'):
                chrom_key = f'chr{chrom}'
            elif chrom_key.startswith('chr'):
                chrom_key = chrom_key[3:]

        if chrom_key not in fasta.keys():
            logger.warning(f"Chromosome {chrom} not found in FASTA file")
            return None

        # Extract trinucleotide (1-based coordinates)
        trinuc = fasta[chrom_key][pos - 2:pos + 1].seq.upper()

        if len(trinuc) != 3 or 'N' in trinuc:
            return None

        return trinuc

    except Exception as e:
        logger.warning(f"Error extracting trinucleotide for {chrom}:{pos}: {e}")
        return None


def _create_context_label(ref: str, alt: str, trinuc: str) -> str:
    """
    Create the 96-context label.

    Parameters
    ----------
    ref : str
        Reference allele
    alt : str
        Alternative allele
    trinuc : str
        Trinucleotide context

    Returns
    -------
    str
        Context label in format "X[REF>ALT]Z"
    """
    return f"{trinuc[0]}[{ref}>{alt}]{trinuc[2]}"


class MutationalSignatureMixin:
    """
    Mixin class providing mutational signature analysis functionality for PyMutation objects.
    
    This mixin adds trinucleotide context matrix generation and mutational signature
    analysis capabilities to PyMutation, following the same architectural pattern 
    as other mixins in the project.
    """

    def trinucleotideMatrix(self, fasta_file: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Generate trinucleotide context matrix for mutational signature analysis.

        This method calculates the 96 trinucleotide contexts for all SNVs in the dataset
        and returns both a contexts matrix (96 x samples) and the original data enriched
        with trinucleotide information.

        Parameters
        ----------
        fasta_file : str
            Path to the reference genome FASTA file

        Returns
        -------
        Tuple[pd.DataFrame, pd.DataFrame]
            - contexts_df: 96 x samples matrix with trinucleotide context counts
            - enriched_data: Original data with added trinuc, class96, and idx96 columns

        Raises
        ------
        ImportError
            If pyfaidx is not installed
        ValueError
            If required columns are missing or no valid SNVs found
        """
        try:
            import pyfaidx
        except ImportError:
            raise ImportError("pyfaidx is required for trinucleotide analysis. Install with: pip install pyfaidx")

        # Get required columns and detect data format
        from ..utils.fields import col

        chrom_col = col(self.data, "Chromosome", required=True)
        pos_col = col(self.data, "Start_Position", required=True)
        ref_col = col(self.data, "Reference_Allele", required=True)
        alt_col = col(self.data, "Tumor_Seq_Allele2", required=True)
        sample_col = col(self.data, "Tumor_Sample_Barcode", required=False)

        if any(col is None for col in [chrom_col, pos_col, ref_col, alt_col]):
            raise ValueError(
                "Required columns not found. Need: Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2")

        if sample_col is None:
            standard_cols = {'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
                             'Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Start_Position',
                             'End_position', 'Strand', 'Variant_Classification', 'Variant_Type',
                             'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2',
                             'Tumor_Sample_Barcode', 'Protein_Change', 'i_TumorVAF_WU', 'i_transcript_name'}

            sample_columns = [col for col in self.data.columns if col not in standard_cols]
            logger.info(f"Detected wide format with {len(sample_columns)} sample columns")
            data_format = "wide"
        else:
            sample_columns = None
            logger.info("Detected long format with Tumor_Sample_Barcode column")
            data_format = "long"

        # Filter for SNVs and load reference genome
        snv_mask = (
                (ref_col.isin(['A', 'C', 'G', 'T'])) &
                (alt_col.isin(['A', 'C', 'G', 'T'])) &
                (ref_col != alt_col)
        )

        if not snv_mask.any():
            raise ValueError("No valid SNVs found in the dataset")

        snv_data = self.data[snv_mask].copy()
        logger.info(f"Processing {len(snv_data)} SNVs from {len(self.data)} total mutations")

        try:
            fasta = pyfaidx.Fasta(fasta_file)
            logger.info(f"Loaded FASTA file: {fasta_file}")
        except Exception as e:
            raise ValueError(f"Error loading FASTA file {fasta_file}: {e}")

        # Process each SNV to extract trinucleotide context
        trinuc_contexts = []
        class96_labels = []
        idx96_values = []

        for idx, row in snv_data.iterrows():
            chrom = str(row[chrom_col.name])
            pos = int(row[pos_col.name])
            ref = str(row[ref_col.name]).upper()
            alt = str(row[alt_col.name]).upper()

            trinuc = _get_trinucleotide_context(fasta, chrom, pos)

            if trinuc is None:
                trinuc_contexts.append(None)
                class96_labels.append(None)
                idx96_values.append(None)
                continue

            norm_ref, norm_alt, norm_trinuc = _normalize_to_pyrimidine(ref, alt, trinuc)
            context_label = _create_context_label(norm_ref, norm_alt, norm_trinuc)
            context_idx = CONTEXT_TO_INDEX.get(context_label)

            trinuc_contexts.append(norm_trinuc)
            class96_labels.append(context_label)
            idx96_values.append(context_idx)

        # Add context information and filter valid entries
        snv_data = snv_data.copy()
        snv_data['trinuc'] = trinuc_contexts
        snv_data['class96'] = class96_labels
        snv_data['idx96'] = idx96_values

        valid_mask = snv_data['idx96'].notna()
        valid_data = snv_data[valid_mask].copy()

        logger.info(f"Successfully processed {len(valid_data)} SNVs with valid trinucleotide contexts")

        if len(valid_data) == 0:
            raise ValueError("No SNVs with valid trinucleotide contexts found")

        # Create 96 x samples matrix based on data format
        if data_format == "long":
            sample_names = valid_data[sample_col.name].unique()
            contexts_matrix = np.zeros((96, len(sample_names)), dtype=int)

            for sample_idx, sample in enumerate(sample_names):
                sample_data = valid_data[valid_data[sample_col.name] == sample]
                context_counts = sample_data['idx96'].value_counts()

                for context_idx, count in context_counts.items():
                    if pd.notna(context_idx):
                        contexts_matrix[int(context_idx), sample_idx] = count

        else:
            sample_names = sample_columns
            contexts_matrix = np.zeros((96, len(sample_names)), dtype=int)

            for idx, row in valid_data.iterrows():
                context_idx = row['idx96']
                if pd.notna(context_idx):
                    context_idx = int(context_idx)

                    for sample_idx, sample_name in enumerate(sample_names):
                        genotype = str(row[sample_name])

                        # Count mutations based on genotype
                        # Genotypes like "G|A", "A|G" indicate heterozygous mutation
                        # Genotypes like "A|A" indicate homozygous mutation
                        # Genotypes like "G|G" indicate no mutation
                        if '|' in genotype:
                            alleles = genotype.split('|')
                            ref_allele = str(row[ref_col.name]).upper()
                            alt_allele = str(row[alt_col.name]).upper()

                            alt_count = sum(1 for allele in alleles if allele == alt_allele)
                            contexts_matrix[context_idx, sample_idx] += alt_count

        # Create contexts DataFrame
        contexts_df = pd.DataFrame(
            contexts_matrix,
            index=TRINUCLEOTIDE_CONTEXTS,
            columns=sample_names
        )

        logger.info(f"Generated {contexts_df.shape[0]} x {contexts_df.shape[1]} trinucleotide context matrix")

        return contexts_df, valid_data


def estimateSignatures(contexts_df: pd.DataFrame, nMin: int = 2, nTry: int = 6,
                       nrun: int = 5, parallel: int = 4, pConstant: Optional[float] = None) -> Dict:
    """
    Estimate optimal number of mutational signatures using NMF decomposition.

    This method normalizes the input matrix to frequencies, performs NMF decomposition
    for different numbers of signatures (k), and calculates stability metrics to
    identify the optimal number of signatures.

    Parameters
    ----------
    contexts_df : pd.DataFrame
        96 x samples matrix with trinucleotide context counts (from trinucleotideMatrix)
    nMin : int, default 2
        Minimum number of signatures to test
    nTry : int, default 6
        Maximum number of signatures to test
    nrun : int, default 5
        Number of NMF runs per k value for stability assessment
    parallel : int, default 4
        Number of CPU cores to use for parallel processing
    pConstant : float, optional
        Small positive constant to add to matrix if NMF fails due to zeros

    Returns
    -------
    Dict
        Dictionary containing:
        - 'metrics': DataFrame with stability metrics for each k
        - 'models': List of NMF models for each k and run
        - 'optimal_k': Suggested optimal number of signatures

    Raises
    ------
    ImportError
        If required packages (scikit-learn, scipy) are not installed
    ValueError
        If input matrix is invalid or NMF consistently fails
    """
    import pandas as pd  # Local import to ensure availability
    try:
        from sklearn.decomposition import NMF
        from sklearn.metrics import mean_squared_error
        from scipy.cluster.hierarchy import linkage, cophenet
        from scipy.spatial.distance import pdist
        import concurrent.futures
    except ImportError as e:
        missing_pkg = str(e).split("'")[1] if "'" in str(e) else "required package"
        raise ImportError(f"{missing_pkg} is required for signature estimation. "
                          f"Install with: pip install scikit-learn scipy")

    # Validate input parameters
    if not isinstance(contexts_df, pd.DataFrame):
        raise ValueError("contexts_df must be a pandas DataFrame")

    if contexts_df.shape[0] != 96:
        raise ValueError("contexts_df must have 96 rows (trinucleotide contexts)")

    if nMin < 1 or nTry < nMin:
        raise ValueError("nMin must be >= 1 and nTry must be >= nMin")

    logger.info(f"Starting signature estimation for k={nMin} to k={nTry} with {nrun} runs each")

    # Normalize matrix to frequencies
    matrix = contexts_df.values.astype(float)
    col_sums = matrix.sum(axis=0, keepdims=True)
    col_sums[col_sums == 0] = 1
    normalized_matrix = matrix / col_sums

    logger.info(f"Normalized matrix shape: {normalized_matrix.shape}")

    # Apply pseudocount if matrix has too many zeros
    zero_fraction = (normalized_matrix == 0).sum() / normalized_matrix.size
    if zero_fraction > 0.8 and pConstant is not None:
        logger.warning(f"Matrix has {zero_fraction:.2%} zeros, adding pConstant={pConstant}")
        normalized_matrix += pConstant
        col_sums = normalized_matrix.sum(axis=0, keepdims=True)
        normalized_matrix = normalized_matrix / col_sums

    def _run_nmf_single(k, run_idx, matrix):
        """Run a single NMF decomposition."""
        try:
            nmf = NMF(n_components=k, init='random', random_state=run_idx,
                      max_iter=1000, tol=1e-4)
            W = nmf.fit_transform(matrix)
            H = nmf.components_

            reconstructed = np.dot(W, H)
            rss = mean_squared_error(matrix, reconstructed) * matrix.size

            return {
                'k': k,
                'run': run_idx,
                'model': nmf,
                'W': W,
                'H': H,
                'rss': rss,
                'success': True
            }
        except Exception as e:
            logger.warning(f"NMF failed for k={k}, run={run_idx}: {e}")
            return {
                'k': k,
                'run': run_idx,
                'model': None,
                'W': None,
                'H': None,
                'rss': np.inf,
                'success': False
            }

    # Execute NMF runs in parallel and collect results
    all_results = []
    k_values = range(nMin, nTry + 1)

    with concurrent.futures.ThreadPoolExecutor(max_workers=parallel) as executor:
        futures = []
        for k in k_values:
            for run_idx in range(nrun):
                future = executor.submit(_run_nmf_single, k, run_idx, normalized_matrix)
                futures.append(future)

        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            all_results.append(result)

    results_by_k = {}
    for result in all_results:
        k = result['k']
        if k not in results_by_k:
            results_by_k[k] = []
        results_by_k[k].append(result)

    # Calculate metrics for each k
    metrics_data = []
    models_data = []

    for k in k_values:
        k_results = results_by_k[k]
        successful_results = [r for r in k_results if r['success']]

        if len(successful_results) == 0:
            logger.warning(f"All NMF runs failed for k={k}")
            continue

        # Calculate stability metrics
        rss_values = [r['rss'] for r in successful_results]
        mean_rss = np.mean(rss_values)
        std_rss = np.std(rss_values)

        # Calculate cophenetic correlation
        if len(successful_results) >= 2:
            H_matrices = [r['H'] for r in successful_results]
            n_samples = H_matrices[0].shape[1]

            consensus_matrix = np.zeros((n_samples, n_samples))

            for i in range(len(H_matrices)):
                for j in range(i + 1, len(H_matrices)):
                    H1, H2 = H_matrices[i], H_matrices[j]
                    corr_matrix = np.corrcoef(H1.T, H2.T)[:n_samples, n_samples:]
                    consensus_matrix += np.abs(corr_matrix)

            consensus_matrix /= (len(H_matrices) * (len(H_matrices) - 1) / 2)

            try:
                distance_matrix = 1 - consensus_matrix
                condensed_dist = pdist(distance_matrix)
                linkage_matrix = linkage(condensed_dist, method='average')
                cophenetic_corr, _ = cophenet(linkage_matrix, condensed_dist)
            except:
                cophenetic_corr = np.nan
        else:
            cophenetic_corr = np.nan

        # Calculate dispersion (coefficient of variation of RSS)
        dispersion = std_rss / mean_rss if mean_rss > 0 else np.inf

        metrics_data.append({
            'k': k,
            'mean_rss': mean_rss,
            'std_rss': std_rss,
            'cophenetic_corr': cophenetic_corr,
            'dispersion': dispersion,
            'successful_runs': len(successful_results),
            'total_runs': len(k_results)
        })

        models_data.extend(successful_results)

    if not metrics_data:
        raise ValueError("All NMF decompositions failed. Try adjusting pConstant or input data.")

    metrics_df = pd.DataFrame(metrics_data)

    # Determine optimal number of signatures
    optimal_k = nMin
    if len(metrics_df) > 1:
        valid_coph = metrics_df.dropna(subset=['cophenetic_corr'])
        if len(valid_coph) > 1:
            coph_diff = valid_coph['cophenetic_corr'].diff().abs()
            if not coph_diff.isna().all():
                max_drop_idx = coph_diff.idxmax()
                optimal_k = valid_coph.loc[max_drop_idx, 'k']

    logger.info(f"Signature estimation completed. Suggested optimal k: {optimal_k}")

    return {
        'metrics': metrics_df,
        'models': models_data,
        'optimal_k': optimal_k,
        'normalized_matrix': normalized_matrix,
        'original_matrix': matrix
    }


def extract_signatures(contexts_df: pd.DataFrame, k: int, nrun: int = 30,
                       pseudocount: float = 1e-4, random_seed: Optional[int] = None) -> Dict:
    """
    Extract mutational signatures using Non-negative Matrix Factorization (NMF).

    This function decomposes a 96 × samples matrix into two interpretable components:
    - Matrix W (96 × k): Signature profiles (each column sums to 1)
    - Matrix H (k × samples): Signature contributions to each sample

    The function converts rows to relative frequencies, adds a small pseudocount to avoid
    zeros, uses NNDSVD initialization with Kullback-Leibler divergence cost function,
    runs the algorithm multiple times with different random seeds, and selects the
    solution with the lowest residual for stability.

    Parameters
    ----------
    contexts_df : pd.DataFrame
        96 × samples matrix with trinucleotide context counts (from trinucleotideMatrix)
    k : int
        Number of signatures to extract (typically obtained from estimate_signatures)
    nrun : int, default 30
        Number of NMF runs with different random seeds to ensure stability
    pseudocount : float, default 1e-4
        Small positive constant added to avoid zeros in the matrix
    random_seed : int, optional
        Base random seed for reproducibility. If None, uses random initialization

    Returns
    -------
    Dict
        Dictionary containing:
        - 'W': numpy.ndarray (96 × k) - Signature profiles (columns sum to 1)
        - 'H': numpy.ndarray (k × samples) - Signature contributions
        - 'reconstruction_error': float - Final reconstruction error
        - 'best_run': int - Index of the best run selected
        - 'all_errors': list - Reconstruction errors from all runs

    Raises
    ------
    ImportError
        If required packages (scikit-learn) are not installed
    ValueError
        If input parameters are invalid or NMF consistently fails

    Notes
    -----
    The W matrix columns are renormalized so each signature sums to 1, facilitating
    comparison with external catalogs. The H matrix indicates non-negative contributions
    of each signature to each tumor, ready for clinical analysis or visualization.
    """
    try:
        from sklearn.decomposition import NMF
        from sklearn.metrics import mean_squared_error
    except ImportError as e:
        missing_pkg = str(e).split("'")[1] if "'" in str(e) else "scikit-learn"
        raise ImportError(f"{missing_pkg} is required for signature extraction. "
                          f"Install with: pip install scikit-learn")

    # Validate input parameters
    if not isinstance(contexts_df, pd.DataFrame):
        raise ValueError("contexts_df must be a pandas DataFrame")

    if contexts_df.shape[0] != 96:
        raise ValueError("contexts_df must have 96 rows (trinucleotide contexts)")

    if k < 1:
        raise ValueError("k must be >= 1")

    if k > contexts_df.shape[1]:
        raise ValueError("k cannot be larger than the number of samples")

    if nrun < 1:
        raise ValueError("nrun must be >= 1")

    if pseudocount <= 0:
        raise ValueError("pseudocount must be > 0")

    logger.info(f"Extracting {k} signatures with {nrun} runs using NMF with KL divergence")

    # Prepare and normalize matrix
    matrix = contexts_df.values.astype(float)

    col_sums = matrix.sum(axis=0, keepdims=True)
    col_sums[col_sums == 0] = 1
    normalized_matrix = matrix / col_sums

    normalized_matrix += pseudocount

    col_sums = normalized_matrix.sum(axis=0, keepdims=True)
    normalized_matrix = normalized_matrix / col_sums

    logger.info(f"Matrix normalized to frequencies with pseudocount {pseudocount}")

    # Execute multiple NMF runs for stability
    all_results = []
    all_errors = []

    if random_seed is not None:
        np.random.seed(random_seed)
        seeds = [random_seed + i for i in range(nrun)]
    else:
        seeds = [np.random.randint(0, 10000) for _ in range(nrun)]

    for run_idx in range(nrun):
        try:
            nmf = NMF(
                n_components=k,
                init='nndsvda',
                solver='mu',
                beta_loss='kullback-leibler',
                random_state=seeds[run_idx],
                max_iter=2000,
                tol=1e-6,
                alpha_W=0.0,
                alpha_H=0.0
            )

            W = nmf.fit_transform(normalized_matrix)
            H = nmf.components_

            # Calculate KL divergence reconstruction error
            reconstructed = np.dot(W, H)

            eps = 1e-10
            original_safe = normalized_matrix + eps
            reconstructed_safe = reconstructed + eps

            kl_div = np.sum(original_safe * np.log(original_safe / reconstructed_safe)
                            - original_safe + reconstructed_safe)

            all_results.append({
                'W': W.copy(),
                'H': H.copy(),
                'error': kl_div,
                'run': run_idx,
                'nmf_model': nmf
            })
            all_errors.append(kl_div)

        except Exception as e:
            logger.warning(f"NMF run {run_idx} failed: {e}")
            all_errors.append(np.inf)
            continue

    # Select best result and normalize signatures
    successful_results = [r for r in all_results if r['error'] != np.inf]
    if not successful_results:
        raise ValueError("All NMF runs failed. Try adjusting pseudocount or input data.")

    best_result = min(successful_results, key=lambda x: x['error'])
    best_run_idx = best_result['run']

    logger.info(f"Best result from run {best_run_idx} with error {best_result['error']:.6f}")
    logger.info(f"Successful runs: {len(successful_results)}/{nrun}")

    W_best = best_result['W']
    H_best = best_result['H']

    W_normalized = W_best / W_best.sum(axis=0, keepdims=True)

    column_sums = W_best.sum(axis=0, keepdims=True)
    H_adjusted = H_best * column_sums.T

    logger.info("Signatures extracted and normalized successfully")

    return {
        'W': W_normalized,
        'H': H_adjusted,
        'reconstruction_error': best_result['error'],
        'best_run': best_run_idx,
        'all_errors': all_errors,
        'successful_runs': len(successful_results),
        'total_runs': nrun
    }


def compare_signatures(W: np.ndarray, cosmic_path: str, min_cosine: float = 0.6,
                       return_matrix: bool = False) -> Dict:
    """
    Compare extracted signatures with COSMIC catalog using cosine similarity.

    This function loads the COSMIC catalog, validates context compatibility,
    calculates cosine similarity between each signature in W and each COSMIC signature,
    and returns a summary of best matches along with optionally the full similarity matrix.

    Parameters
    ----------
    W : numpy.ndarray
        96 × k matrix from extract_signatures with normalized signature profiles
        (each column sums to 1)
    cosmic_path : str
        Path to COSMIC catalog file (e.g., "COSMIC_v3.4_SBS_GRCh38.txt")
    min_cosine : float, default 0.6
        Minimum cosine similarity threshold for considering a match
    return_matrix : bool, default False
        If True, also return the full cosine similarity matrix

    Returns
    -------
    Dict
        Dictionary containing:
        - 'summary_df': pandas.DataFrame with columns ['Signature_W', 'Best_COSMIC', 'Cosine', 'Aetiology']
        - 'cosine_matrix': numpy.ndarray (k × N) - Only if return_matrix=True

    Raises
    ------
    FileNotFoundError
        If cosmic_path file does not exist
    ValueError
        If contexts don't match between W and COSMIC catalog
    ImportError
        If required packages are not installed

    Notes
    -----
    The function expects the COSMIC catalog to have 96 rows corresponding to 
    trinucleotide contexts and columns for different signatures. The first column
    should contain context labels, and there should be an 'aetiology' column or
    similar metadata.
    """
    try:
        from sklearn.metrics.pairwise import cosine_similarity
    except ImportError as e:
        missing_pkg = str(e).split("'")[1] if "'" in str(e) else "scikit-learn"
        raise ImportError(f"{missing_pkg} is required for cosine similarity calculation. "
                          f"Install with: pip install scikit-learn")

    # Validate input and load COSMIC catalog
    if not isinstance(W, np.ndarray):
        raise ValueError("W must be a numpy array")

    if W.shape[0] != 96:
        raise ValueError("W must have 96 rows (trinucleotide contexts)")

    if len(W.shape) != 2:
        raise ValueError("W must be a 2D array")

    logger.info(f"Comparing {W.shape[1]} signatures with COSMIC catalog")

    try:
        cosmic_df = pd.read_csv(cosmic_path, sep='\t')
    except FileNotFoundError:
        raise FileNotFoundError(f"COSMIC catalog file not found: {cosmic_path}")
    except Exception as e:
        raise ValueError(f"Error reading COSMIC catalog: {e}")

    logger.info(f"Loaded COSMIC catalog with shape {cosmic_df.shape}")

    if cosmic_df.shape[0] != 96:
        raise ValueError(f"COSMIC catalog must have 96 rows (trinucleotide contexts), got {cosmic_df.shape[0]}")

    # Process COSMIC catalog structure
    context_column = cosmic_df.columns[0]  # First column should be 'Type'
    contexts = cosmic_df[context_column].values  # All rows are data (no header row to skip)
    cosmic_df = cosmic_df.set_index(context_column)

    # Normalize W matrix and align with COSMIC catalog
    W_sums = W.sum(axis=0)
    if not np.allclose(W_sums, 1.0, rtol=1e-5):
        logger.warning("W matrix columns do not sum to 1, renormalizing...")
        W = W / W_sums.reshape(1, -1)
        logger.info("W matrix renormalized")

    logger.info("Aligning COSMIC catalog with standard trinucleotide context order...")

    missing_contexts = set(TRINUCLEOTIDE_CONTEXTS) - set(contexts)
    if missing_contexts:
        raise ValueError(f"COSMIC catalog is missing required contexts: {missing_contexts}")

    try:
        cosmic_df = cosmic_df.reindex(TRINUCLEOTIDE_CONTEXTS)
        if cosmic_df.isnull().any().any():
            raise ValueError("Some standard trinucleotide contexts are missing in COSMIC catalog")
        logger.info("COSMIC catalog successfully aligned to standard context order")
    except Exception as e:
        raise ValueError(f"Failed to align COSMIC catalog contexts: {e}")

    signature_columns = list(cosmic_df.columns)

    if len(signature_columns) == 0:
        raise ValueError("No signature columns found in COSMIC catalog")

    # Filter out artifact signatures
    artifact_signatures = {
        'SBS27', 'SBS43', 'SBS45', 'SBS46', 'SBS47', 'SBS48', 'SBS49', 'SBS50',
        'SBS51', 'SBS52', 'SBS53', 'SBS54', 'SBS55', 'SBS56', 'SBS57', 'SBS58',
        'SBS59', 'SBS60', 'SBS95'
    }

    filtered_columns = []
    for col in signature_columns:
        if col.endswith('c'):
            logger.info(f"Removing artifact signature {col} (ends with 'c')")
            continue
        if 'artefact' in col.lower():
            logger.info(f"Removing artifact signature {col} (contains 'artefact')")
            continue
        if col in artifact_signatures:
            logger.info(f"Removing artifact signature {col} (specified artifact)")
            continue
        filtered_columns.append(col)

    signature_columns = filtered_columns
    cosmic_df = cosmic_df[signature_columns]

    cosmic_matrix = cosmic_df.values.astype(float)

    column_sums = cosmic_matrix.sum(axis=0)
    zero_sum_mask = column_sums == 0
    if zero_sum_mask.any():
        zero_sum_sigs = [signature_columns[i] for i in range(len(signature_columns)) if zero_sum_mask[i]]
        logger.info(f"Removing signatures with zero sum: {zero_sum_sigs}")
        non_zero_mask = ~zero_sum_mask
        cosmic_matrix = cosmic_matrix[:, non_zero_mask]
        signature_columns = [sig for i, sig in enumerate(signature_columns) if non_zero_mask[i]]

    logger.info(f"Found {len(signature_columns)} valid COSMIC signatures after filtering")

    # Validate alignment and calculate cosine similarity
    if len(TRINUCLEOTIDE_CONTEXTS) != 96:
        raise ValueError(f"Expected 96 trinucleotide contexts in standard order, got {len(TRINUCLEOTIDE_CONTEXTS)}")

    if cosmic_matrix.shape[0] != 96:
        raise ValueError(f"COSMIC matrix must have 96 rows, got {cosmic_matrix.shape[0]}")

    current_contexts = cosmic_df.index.tolist()
    if current_contexts != TRINUCLEOTIDE_CONTEXTS:
        raise ValueError("COSMIC catalog contexts are not properly aligned with standard order")

    cosmic_normalized = cosmic_matrix / cosmic_matrix.sum(axis=0, keepdims=True)

    logger.info("COSMIC signatures normalized")

    cosine_matrix = cosine_similarity(W.T, cosmic_normalized.T)

    logger.info(f"Calculated cosine similarity matrix: {cosine_matrix.shape}")

    # Create comparison summary
    summary_data = []

    for i in range(W.shape[1]):
        similarities = cosine_matrix[i, :]
        best_idx = np.argmax(similarities)
        best_cosine = similarities[best_idx]
        best_cosmic = signature_columns[best_idx]

        if best_cosine >= min_cosine:
            match_status = best_cosmic
        else:
            match_status = "No match"

        aetiology = "Unknown"

        summary_data.append({
            'Signature_W': f'Signature_{i + 1}',
            'Best_COSMIC': match_status,
            'Cosine': best_cosine,
            'Aetiology': aetiology
        })

    summary_df = pd.DataFrame(summary_data)

    logger.info(f"Created summary with {len(summary_data)} signature comparisons")

    result = {'summary_df': summary_df}

    if return_matrix:
        result['cosine_matrix'] = cosine_matrix

    return result
