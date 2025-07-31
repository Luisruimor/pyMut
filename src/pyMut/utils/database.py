import gzip
import hashlib
from datetime import datetime
from pathlib import Path

import duckdb
import pandas as pd


class PfamAnnotationError(Exception):
    """Custom exception for Pfam annotation errors."""
    pass


def get_resources_path() -> Path:
    """Get the path to the resources directory."""
    return Path(__file__).parent.parent / "data" / "resources"


def get_db_path() -> Path:
    """Get the path to the DuckDB database file."""
    return get_resources_path() / "data.duckdb"


def calculate_file_hash(filepath: str) -> str:
    """Calculate SHA-256 hash of a file."""
    hash_sha256 = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_sha256.update(chunk)
    return hash_sha256.hexdigest()


def parse_idmapping_selected(mapping_file: Path, chunk_size: int = int(1e6)) -> pd.DataFrame:
    """
    Parse UniProt ID mapping file with column-based extraction.

    This function reads the UniProt ID mapping file and extracts mappings
    for Ensembl proteins (ENSP) and RefSeq proteins (NP_) to UniProt IDs,
    as well as short names from UniProt.

    Expected format:
    - Column 1: UniProt ID
    - Column 2: Short name (e.g., 1433B_HUMAN)
    - Column 4: RefSeq protein IDs (semicolon-separated)
    - Column 21: Ensembl protein IDs (semicolon-separated)

    Args:
        mapping_file: Path to the ID mapping file
        chunk_size: Size of chunks to process at a time

    Returns:
        DataFrame with columns: prot_id, uniprot, short_name
    """
    print(f"Parsing UniProt ID mapping file: {mapping_file}")

    all_mappings = []
    total_lines = 0
    kept_lines = 0

    try:
        with gzip.open(mapping_file, 'rt') as f:
            chunk_data = []

            for line_num, line in enumerate(f, 1):
                total_lines += 1

                parts = line.strip().split('\t')
                if len(parts) >= 21:  # Need at least 21 columns
                    uniprot_id = parts[0]
                    short_name = parts[1] if len(parts) > 1 and parts[1] != "-" else ""

                    # RefSeq protein IDs from column 4
                    refseq_ids = parts[3] if len(parts) > 3 else ""
                    if refseq_ids and refseq_ids != "-":
                        for refseq_id in refseq_ids.split(';'):
                            refseq_id = refseq_id.strip()
                            if refseq_id.startswith('NP_'):
                                chunk_data.append({
                                    'prot_id': refseq_id,
                                    'uniprot': uniprot_id,
                                    'short_name': short_name
                                })
                                kept_lines += 1

                    # Ensembl protein IDs from column 21
                    ensembl_ids = parts[20] if len(parts) > 20 else ""
                    if ensembl_ids and ensembl_ids != "-":
                        for ensembl_id in ensembl_ids.split(';'):
                            ensembl_id = ensembl_id.strip()
                            if ensembl_id.startswith('ENSP'):
                                chunk_data.append({
                                    'prot_id': ensembl_id,
                                    'uniprot': uniprot_id,
                                    'short_name': short_name
                                })
                                kept_lines += 1

                    # Add row for short_name association if short_name exists
                    if short_name:
                        chunk_data.append({
                            'prot_id': short_name,
                            'uniprot': uniprot_id,
                            'short_name': short_name
                        })
                        kept_lines += 1

                if len(chunk_data) >= chunk_size:
                    all_mappings.extend(chunk_data)
                    chunk_data = []

                    if line_num % (chunk_size * 10) == 0:
                        print(f"  Processed {line_num:,} lines, kept {kept_lines:,} mappings...")

            if chunk_data:
                all_mappings.extend(chunk_data)

        print(f"Processed {total_lines:,} total lines")
        print(f"Kept {kept_lines:,} protein mappings")

        if all_mappings:
            df = pd.DataFrame(all_mappings)
            df = df.drop_duplicates()
            print(f"Final unique mappings: {len(df):,}")
            return df
        else:
            print("⚠️  No relevant mappings found")
            return pd.DataFrame(columns=['prot_id', 'uniprot', 'short_name'])

    except Exception as e:
        print(f"❌ Error parsing mapping file: {e}")
        return pd.DataFrame(columns=['prot_id', 'uniprot', 'short_name'])


def check_mapping_coverage(conn: duckdb.DuckDBPyConnection) -> bool:
    """
    Test if the xref table has reasonable coverage for common protein ID types.

    Args:
        conn: DuckDB connection

    Returns:
        True if coverage looks reasonable, False otherwise
    """
    try:
        ensp_count = conn.execute("SELECT COUNT(*) FROM xref WHERE prot_id LIKE 'ENSP%'").fetchone()[0]
        np_count = conn.execute("SELECT COUNT(*) FROM xref WHERE prot_id LIKE 'NP_%'").fetchone()[0]
        total_count = conn.execute("SELECT COUNT(*) FROM xref").fetchone()[0]

        print("Mapping coverage:")
        print(f"    Total mappings: {total_count:,}")
        print(f"    ENSP mappings: {ensp_count:,}")
        print(f"    NP_ mappings: {np_count:,}")

        min_threshold = 1000  # Minimum expected mappings
        ensp_ok = ensp_count >= min_threshold
        np_ok = np_count >= min_threshold

        if ensp_ok and np_ok:
            print("✅ Mapping coverage looks good")
            return True
        else:
            print("⚠️  Low mapping coverage detected:")
            if not ensp_ok:
                print(f"    ENSP count {ensp_count:,} < {min_threshold:,}")
            if not np_ok:
                print(f"    NP_ count {np_count:,} < {min_threshold:,}")
            return False

    except Exception as e:
        print(f"❌ Error during mapping coverage test: {e}")
        return False


def build_embedded_db(force_rebuild: bool = False) -> str:
    """
    Build embedded DuckDB database with Pfam and mapping data.

    Supported files:
    - idmapping_selected.tab.gz (global) or HUMAN_9606_idmapping_selected.tab.gz (organism-specific)
    - Detection by content prefix (ENSP/NP_), not fixed column positions
    - Pfam-A.regions.tsv.gz for domain annotations

    Args:
        force_rebuild: If True, rebuild the database even if it exists

    Returns:
        Path to the created database file
    """
    db_path = get_db_path()
    resources_path = get_resources_path()

    if db_path.exists() and not force_rebuild:
        try:
            conn = duckdb.connect(str(db_path))
            tables = conn.execute("SHOW TABLES").fetchall()
            table_names = [table[0] for table in tables]
            if 'pfam' in table_names and 'xref' in table_names and 'meta' in table_names:
                print(f"✅ Database already exists at {db_path}")
                conn.close()
                return str(db_path)
            conn.close()
        except Exception:
            pass

    print("🔨 Building embedded DuckDB database...")

    conn = duckdb.connect(str(db_path))

    # 1. Read Pfam-A.regions.tsv.gz and create pfam table
    pfam_file = resources_path / "pfam" / "Pfam-A.regions.tsv.gz"
    if not pfam_file.exists():
        raise PfamAnnotationError(f"Pfam file not found: {pfam_file}")

    print("Loading Pfam data...")

    chunk_size = int(1e6)
    total_rows = 0

    conn.execute("DROP TABLE IF EXISTS pfam")
    conn.execute("""
                 CREATE TABLE pfam
                 (
                     uniprot   VARCHAR,
                     seq_start INTEGER,
                     seq_end   INTEGER,
                     pfam_id   VARCHAR,
                     pfam_name VARCHAR
                 )
                 """)

    with gzip.open(pfam_file, 'rt') as f:
        header = f.readline().strip().split('\t')
        print(f"Pfam file columns: {header}")

        f.seek(0)

        chunk_count = 0
        for chunk in pd.read_csv(f, sep='\t', chunksize=chunk_size):
            chunk_count += 1

            # Based on file structure: ['pfamseq_acc', 'seq_version', 'crc64', 'md5', 'pfamA_acc', 'seq_start', 'seq_end', 'ali_start', 'ali_end']
            column_mapping = {
                'pfamseq_acc': 'uniprot',
                'pfamA_acc': 'pfam_id',
                'seq_start': 'seq_start',
                'seq_end': 'seq_end'
            }

            chunk = chunk.rename(columns=column_mapping)

            # Use pfam_id as name placeholder
            if 'pfam_name' not in chunk.columns and 'pfam_id' in chunk.columns:
                chunk['pfam_name'] = chunk['pfam_id']

            required_cols = ['uniprot', 'seq_start', 'seq_end', 'pfam_id', 'pfam_name']
            available_cols = [col for col in required_cols if col in chunk.columns]

            if len(available_cols) == 5:
                chunk_selected = chunk[available_cols].copy()

                conn.register(f'pfam_chunk_{chunk_count}', chunk_selected)
                conn.execute(f"INSERT INTO pfam SELECT * FROM pfam_chunk_{chunk_count}")
                conn.unregister(f'pfam_chunk_{chunk_count}')

                total_rows += len(chunk_selected)

                del chunk_selected
                del chunk
            else:
                print(f"Warning: Missing columns. Available: {available_cols}, Required: {required_cols}")

            if chunk_count % 10 == 0:
                print(f"  Processed {chunk_count} chunks...")

    print(f"Loaded {total_rows:,} Pfam domain annotations")

    # 2. Read HUMAN_9606_idmapping_selected.tab.gz and create xref table
    mapping_file = resources_path / "mappings" / "HUMAN_9606_idmapping_selected.tab.gz"
    if not mapping_file.exists():
        print(f"Mapping file not found: {mapping_file}")
        print("Creating empty xref table.")
        print("📝  Note: Without UniProt mappings, Pfam annotation will be limited to variants")
        print("    that already have UniProt protein IDs in the input data.")
        conn.execute("DROP TABLE IF EXISTS xref")
        conn.execute("""
                     CREATE TABLE xref
                     (
                         prot_id    VARCHAR,
                         uniprot    VARCHAR,
                         short_name VARCHAR
                     )
                     """)
    else:
        file_size = mapping_file.stat().st_size
        if file_size < 1000:  # Less than 1KB suggests empty or corrupted file
            print(f"⚠️  Mapping file appears to be empty or corrupted (size: {file_size} bytes)")
            print(f"⚠️  Expected file location: {mapping_file}")
            print("⚠️  Creating empty xref table.")
            print("📝  Note: To enable full Pfam annotation functionality, please download")
            print("    the UniProt ID mapping file from:")
            print(
                "    https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz")
            print("    and place it at the expected location.")
            conn.execute("DROP TABLE IF EXISTS xref")
            conn.execute("""
                         CREATE TABLE xref
                         (
                             prot_id    VARCHAR,
                             uniprot    VARCHAR,
                             short_name VARCHAR
                         )
                         """)
        else:
            try:
                xref_df = parse_idmapping_selected(mapping_file, chunk_size)

                if len(xref_df) > 0:
                    conn.execute("DROP TABLE IF EXISTS xref")
                    conn.execute("""
                                 CREATE TABLE xref
                                 (
                                     prot_id    VARCHAR,
                                     uniprot    VARCHAR,
                                     short_name VARCHAR
                                 )
                                 """)

                    conn.register('xref_temp', xref_df)
                    conn.execute("INSERT INTO xref SELECT * FROM xref_temp")
                    conn.unregister('xref_temp')

                    check_mapping_coverage(conn)
                else:
                    print("⚠️  No mapping data found after filtering")
                    print("⚠️  Creating empty xref table. Pfam annotation may be limited.")
                    conn.execute("DROP TABLE IF EXISTS xref")
                    conn.execute("""
                                 CREATE TABLE xref
                                 (
                                     prot_id    VARCHAR,
                                     uniprot    VARCHAR,
                                     short_name VARCHAR
                                 )
                                 """)
            except Exception as e:
                print(f"⚠️  Error reading mapping file: {e}")
                print("⚠️  Creating empty xref table. Pfam annotation may be limited.")
                conn.execute("DROP TABLE IF EXISTS xref")
                conn.execute("""
                             CREATE TABLE xref
                             (
                                 prot_id    VARCHAR,
                                 uniprot    VARCHAR,
                                 short_name VARCHAR
                             )
                             """)

    # 3. Create indices
    print("🔍 Creating database indices...")
    conn.execute("CREATE INDEX IF NOT EXISTS ix_pfam ON pfam(uniprot, seq_start, seq_end)")
    conn.execute("CREATE INDEX IF NOT EXISTS ix_xref_prot_id ON xref(prot_id)")
    conn.execute("CREATE INDEX IF NOT EXISTS ix_xref_uniprot ON xref(uniprot)")
    conn.execute("CREATE INDEX IF NOT EXISTS ix_xref_short_name ON xref(short_name)")

    # 4. Create metadata table
    print("📝 Creating metadata table...")
    conn.execute("DROP TABLE IF EXISTS meta")
    conn.execute("""
                 CREATE TABLE meta
                 (
                     resource     VARCHAR,
                     file_path    VARCHAR,
                     release_date VARCHAR,
                     sha256_hash  VARCHAR,
                     created_at   TIMESTAMP
                 )
                 """)

    metadata_entries = [
        {
            'resource': 'pfam',
            'file_path': str(pfam_file),
            'release_date': 'unknown',
            'sha256_hash': calculate_file_hash(str(pfam_file)),
            'created_at': datetime.now()
        },
        {
            'resource': 'uniprot_mapping',
            'file_path': str(mapping_file),
            'release_date': 'unknown',
            'sha256_hash': calculate_file_hash(str(mapping_file)),
            'created_at': datetime.now()
        }
    ]

    meta_df = pd.DataFrame(metadata_entries)
    conn.register('meta_temp', meta_df)
    conn.execute("INSERT INTO meta SELECT * FROM meta_temp")
    conn.unregister('meta_temp')

    conn.close()

    print(f"✅ Database created successfully at {db_path}")
    return str(db_path)


def connect_db() -> duckdb.DuckDBPyConnection:
    """Connect to the embedded DuckDB database."""
    db_path = get_db_path()
    if not db_path.exists():
        print("Database not found. Building it now...")
        build_embedded_db()

    return duckdb.connect(str(db_path))
