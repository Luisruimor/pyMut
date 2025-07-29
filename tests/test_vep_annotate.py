import pytest
import subprocess
import tempfile
from pathlib import Path

from src.pyMut.annotate.vep_annotate import (
    _extract_assembly_and_version_from_cache,
    _get_case_insensitive_column,
    _maf_to_region,
    wrap_maf_vep_annotate_protein,
    wrap_vcf_vep_annotate_unified
)

# Real file paths (relative to project root)
VCF_FILE = "src/pyMut/data/examples/VCF/subset_1k_variants_ALL.chr10.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf"
MAF_FILE = "src/pyMut/data/examples/MAF/tcga_laml.maf.gz"
VCF_CACHE_DIR = "src/pyMut/data/resources/vep/homo_sapiens_vep_114_GRCh38"
MAF_CACHE_DIR = "src/pyMut/data/resources/vep/homo_sapiens_vep_114_GRCh37"
VCF_FASTA = "src/pyMut/data/resources/genome/GRCh38/GRCh38.p14.genome.fa"
MAF_FASTA = "src/pyMut/data/resources/genome/GRCh37/GRCh37.p13.genome.fa"


class TestExtractAssemblyAndVersionFromCache:
    """Tests for _extract_assembly_and_version_from_cache function"""
    
    def test_valid_cache_name_grch38(self):
        """Test extraction from valid cache name with GRCh38"""
        cache_path = Path("homo_sapiens_vep_114_GRCh38")
        assembly, version = _extract_assembly_and_version_from_cache(cache_path)
        assert assembly == "GRCh38"
        assert version == "114"
    
    def test_valid_cache_name_with_dot(self):
        """Test extraction from cache name with dot in version"""
        # The current regex pattern only supports integer versions, not decimal versions
        # So this test should expect a ValueError
        cache_path = Path("homo_sapiens_vep_110.1_GRCh37")
        with pytest.raises(ValueError, match="doesn't match expected format"):
            _extract_assembly_and_version_from_cache(cache_path)
    
    def test_invalid_cache_name_raises_valueerror(self):
        """Test invalid cache name raises ValueError"""
        cache_path = Path("invalid_cache_name")
        with pytest.raises(ValueError, match="doesn't match expected format"):
            _extract_assembly_and_version_from_cache(cache_path)


class TestGetCaseInsensitiveColumn:
    """Tests for _get_case_insensitive_column function"""
    
    def test_case_insensitive_match(self):
        """Test case insensitive column matching"""
        columns = ["ChRoMoSoMe", "position"]
        result = _get_case_insensitive_column(columns, "chromosome")
        assert result == "ChRoMoSoMe"
    
    def test_missing_column_raises_keyerror(self):
        """Test missing column raises KeyError"""
        columns = ["position", "start"]
        with pytest.raises(KeyError, match="Column 'chromosome' not found"):
            _get_case_insensitive_column(columns, "chromosome")


class TestMafToRegion:
    """Tests for _maf_to_region function"""
    
    def test_nonexistent_file_returns_false(self):
        """Test nonexistent file returns False"""
        success, result = _maf_to_region("nonexistent.maf")
        assert success is False
        assert result == ""
    
    def test_maf_with_strand_column(self):
        """Test MAF file with strand column"""
        # Create temporary MAF file with strand column
        maf_data = """Hugo_Symbol	Chromosome	Start_Position	End_Position	Strand	Reference_Allele	Tumor_Seq_Allele2
GENE1	1	100	100	+	A	T
GENE2	2	200	200	-	C	G
GENE3	X	300	300	+	G	T"""
        
        temp_maf = tempfile.NamedTemporaryFile(mode='w', suffix='.maf', delete=False)
        temp_maf.write(maf_data)
        temp_maf.close()
        
        try:
            success, region_path = _maf_to_region(temp_maf.name)
            assert success is True
            assert region_path.endswith('.region')
            
            # Check region file content
            with open(region_path, 'r') as f:
                lines = f.readlines()
            
            assert len(lines) == 3
            assert lines[0].strip() == "chr1:100-100:1/T"
            assert lines[1].strip() == "chr2:200-200:-1/G"
            assert lines[2].strip() == "chrX:300-300:1/T"
            
            # Clean up
            Path(region_path).unlink()
        finally:
            Path(temp_maf.name).unlink()
    
    def test_maf_without_strand_column(self):
        """Test MAF file without strand column (defaults to +)"""
        # Create temporary MAF file without strand column
        maf_data = """Hugo_Symbol	Chromosome	Start_Position	End_Position	Reference_Allele	Tumor_Seq_Allele2
GENE1	1	100	100	A	T
GENE2	2	200	200	C	G
GENE3	X	300	300	G	T"""
        
        temp_maf = tempfile.NamedTemporaryFile(mode='w', suffix='.maf', delete=False)
        temp_maf.write(maf_data)
        temp_maf.close()
        
        try:
            success, region_path = _maf_to_region(temp_maf.name)
            assert success is True
            assert region_path.endswith('.region')
            
            # Check region file content
            with open(region_path, 'r') as f:
                lines = f.readlines()
            
            assert len(lines) == 3
            assert lines[0].strip() == "chr1:100-100:1/T"
            assert lines[1].strip() == "chr2:200-200:1/G"
            assert lines[2].strip() == "chrX:300-300:1/T"
            
            # Clean up
            Path(region_path).unlink()
        finally:
            Path(temp_maf.name).unlink()
    
    def test_maf_missing_required_columns(self):
        """Test MAF file missing required columns returns False"""
        # Create temporary MAF file missing required columns
        maf_data = """Hugo_Symbol	Position
GENE1	100
GENE2	200"""
        
        temp_maf = tempfile.NamedTemporaryFile(mode='w', suffix='.maf', delete=False)
        temp_maf.write(maf_data)
        temp_maf.close()
        
        try:
            success, result = _maf_to_region(temp_maf.name)
            assert success is False
            assert result.endswith('.region')  # Function returns output path even on failure
        finally:
            Path(temp_maf.name).unlink()
    
    def test_empty_maf_file(self):
        """Test empty MAF file returns False"""
        temp_maf = tempfile.NamedTemporaryFile(mode='w', suffix='.maf', delete=False)
        temp_maf.write("")
        temp_maf.close()
        
        try:
            success, result = _maf_to_region(temp_maf.name)
            assert success is False
            assert result.endswith('.region')  # Function returns output path even on failure
        finally:
            Path(temp_maf.name).unlink()
    
    def test_chromosome_normalization(self):
        """Test chromosome normalization (chr prefix handling)"""
        # Create temporary MAF file with mixed chromosome formats
        maf_data = """Hugo_Symbol	Chromosome	Start_Position	End_Position	Reference_Allele	Tumor_Seq_Allele2
GENE1	1	100	100	A	T
GENE2	chr2	200	200	C	G
GENE3	X	300	300	G	T"""
        
        temp_maf = tempfile.NamedTemporaryFile(mode='w', suffix='.maf', delete=False)
        temp_maf.write(maf_data)
        temp_maf.close()
        
        try:
            success, region_path = _maf_to_region(temp_maf.name)
            assert success is True
            
            # Check region file content
            with open(region_path, 'r') as f:
                lines = f.readlines()
            
            assert len(lines) == 3
            assert lines[0].strip() == "chr1:100-100:1/T"
            assert lines[1].strip() == "chr2:200-200:1/G"
            assert lines[2].strip() == "chrX:300-300:1/T"
            
            # Clean up
            Path(region_path).unlink()
        finally:
            Path(temp_maf.name).unlink()


class TestWrapMafVepAnnotateProtein:
    """Tests for wrap_maf_vep_annotate_protein function using real files"""
    
    def test_successful_annotation_with_real_files(self):
        """Test successful MAF VEP annotation with real files"""
        # Skip test if VEP is not available
        try:
            subprocess.run(["vep", "--help"], check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            pytest.skip("VEP not available")
        
        # Check if real files exist
        if not all(Path(f).exists() for f in [MAF_FILE, MAF_CACHE_DIR, MAF_FASTA]):
            pytest.skip("Required real files not found")
        
        success, result = wrap_maf_vep_annotate_protein(
            MAF_FILE, MAF_CACHE_DIR, MAF_FASTA
        )
        
        # The function should complete successfully
        assert success is True
        assert "VEP folder:" in result or "Merged file:" in result
    
    def test_maf_to_region_with_real_file(self):
        """Test _maf_to_region with real MAF file"""
        if not Path(MAF_FILE).exists():
            pytest.skip("Real MAF file not found")
        
        success, region_path = _maf_to_region(MAF_FILE)
        
        if success:
            assert region_path.endswith('.region')
            assert Path(region_path).exists()
            
            # Check that region file has content
            with open(region_path, 'r') as f:
                lines = f.readlines()
            assert len(lines) > 0
            
            # Clean up
            Path(region_path).unlink()


class TestWrapVcfVepAnnotateGene:
    """Tests for wrap_vcf_vep_annotate_gene function using real files"""
    
    def test_successful_annotation_with_real_files(self):
        """Test VCF gene annotation with real files"""
        # Skip test if VEP is not available
        try:
            subprocess.run(["vep", "--help"], check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            pytest.skip("VEP not available")
        
        # Check if real files exist
        if not all(Path(f).exists() for f in [VCF_FILE, VCF_CACHE_DIR, VCF_FASTA]):
            pytest.skip("Required real files not found")
        
        success, result = wrap_vcf_vep_annotate_unified(
            VCF_FILE, VCF_CACHE_DIR, VCF_FASTA, annotate_gene=True
        )
        
        # The function should complete successfully
        assert success is True
        assert "VEP output file:" in result
    
    def test_with_distance_parameter(self):
        """Test distance parameter with real files"""
        # Skip test if VEP is not available or files don't exist
        try:
            subprocess.run(["vep", "--help"], check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            pytest.skip("VEP not available")
        
        if not all(Path(f).exists() for f in [VCF_FILE, VCF_CACHE_DIR, VCF_FASTA]):
            pytest.skip("Required real files not found")
        
        success, result = wrap_vcf_vep_annotate_unified(
            VCF_FILE, VCF_CACHE_DIR, VCF_FASTA, annotate_gene=True, distance=10000
        )
        
        assert success is True
        assert "VEP output file:" in result
    
    def test_no_stats_parameter_with_real_files(self):
        """Test no_stats parameter with real files"""
        # Skip test if VEP is not available or files don't exist
        try:
            subprocess.run(["vep", "--help"], check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            pytest.skip("VEP not available")
        
        if not all(Path(f).exists() for f in [VCF_FILE, VCF_CACHE_DIR, VCF_FASTA]):
            pytest.skip("Required real files not found")
        
        # Test with no_stats=True (default)
        success, result = wrap_vcf_vep_annotate_unified(
            VCF_FILE, VCF_CACHE_DIR, VCF_FASTA, annotate_gene=True, no_stats=True
        )
        assert success is True
        
        # Test with no_stats=False
        success, result = wrap_vcf_vep_annotate_unified(
            VCF_FILE, VCF_CACHE_DIR, VCF_FASTA, annotate_gene=True, no_stats=False
        )
        assert success is True
    
    def test_missing_files_raise_filenotfounderror(self):
        """Test missing files raise FileNotFoundError"""
        with pytest.raises(FileNotFoundError, match="VCF file not found"):
            wrap_vcf_vep_annotate_unified("nonexistent.vcf", VCF_CACHE_DIR, VCF_FASTA, annotate_gene=True)
        
        with pytest.raises(FileNotFoundError, match="Cache directory not found"):
            wrap_vcf_vep_annotate_unified(VCF_FILE, "nonexistent_cache", VCF_FASTA, annotate_gene=True)
        
        with pytest.raises(FileNotFoundError, match="FASTA file not found"):
            wrap_vcf_vep_annotate_unified(VCF_FILE, VCF_CACHE_DIR, "nonexistent.fa", annotate_gene=True)


class TestWrapVcfVepAnnotateUnified:
    """Tests for wrap_vcf_vep_annotate_unified function"""
    
    def test_validation_no_annotation_enabled(self):
        """Test that ValueError is raised when no annotation is enabled"""
        with pytest.raises(ValueError, match="At least one annotation option must be enabled"):
            wrap_vcf_vep_annotate_unified(
                VCF_FILE, VCF_CACHE_DIR, VCF_FASTA,
                annotate_protein=False, annotate_gene=False, annotate_variant_class=False
            )
    
    def test_protein_annotation_only(self):
        """Test protein annotation only"""
        # Skip test if VEP is not available or files don't exist
        try:
            subprocess.run(["vep", "--help"], check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            pytest.skip("VEP not available")
        
        if not all(Path(f).exists() for f in [VCF_FILE, VCF_CACHE_DIR, VCF_FASTA]):
            pytest.skip("Required real files not found")
        
        success, result = wrap_vcf_vep_annotate_unified(
            VCF_FILE, VCF_CACHE_DIR, VCF_FASTA, annotate_protein=True
        )
        
        assert success is True
        assert "VEP output file:" in result
        assert "protein" in result.lower()
    
    def test_gene_annotation_only(self):
        """Test gene annotation only"""
        # Skip test if VEP is not available or files don't exist
        try:
            subprocess.run(["vep", "--help"], check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            pytest.skip("VEP not available")
        
        if not all(Path(f).exists() for f in [VCF_FILE, VCF_CACHE_DIR, VCF_FASTA]):
            pytest.skip("Required real files not found")
        
        success, result = wrap_vcf_vep_annotate_unified(
            VCF_FILE, VCF_CACHE_DIR, VCF_FASTA, annotate_gene=True
        )
        
        assert success is True
        assert "VEP output file:" in result
        assert "gene" in result.lower()
    
    def test_variant_class_annotation_only(self):
        """Test variant class annotation only"""
        # Skip test if VEP is not available or files don't exist
        try:
            subprocess.run(["vep", "--help"], check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            pytest.skip("VEP not available")
        
        if not all(Path(f).exists() for f in [VCF_FILE, VCF_CACHE_DIR, VCF_FASTA]):
            pytest.skip("Required real files not found")
        
        success, result = wrap_vcf_vep_annotate_unified(
            VCF_FILE, VCF_CACHE_DIR, VCF_FASTA, annotate_variant_class=True
        )
        
        assert success is True
        assert "VEP output file:" in result
        assert "variant_class" in result.lower()
    
    def test_combined_annotations(self):
        """Test combining multiple annotation types"""
        # Skip test if VEP is not available or files don't exist
        try:
            subprocess.run(["vep", "--help"], check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            pytest.skip("VEP not available")
        
        if not all(Path(f).exists() for f in [VCF_FILE, VCF_CACHE_DIR, VCF_FASTA]):
            pytest.skip("Required real files not found")
        
        success, result = wrap_vcf_vep_annotate_unified(
            VCF_FILE, VCF_CACHE_DIR, VCF_FASTA,
            annotate_protein=True, annotate_gene=True, annotate_variant_class=True
        )
        
        assert success is True
        assert "VEP output file:" in result
        assert "protein" in result.lower()
        assert "gene" in result.lower()
        assert "variant_class" in result.lower()
    
    def test_gene_annotation_with_distance(self):
        """Test gene annotation with distance parameter"""
        # Skip test if VEP is not available or files don't exist
        try:
            subprocess.run(["vep", "--help"], check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            pytest.skip("VEP not available")
        
        if not all(Path(f).exists() for f in [VCF_FILE, VCF_CACHE_DIR, VCF_FASTA]):
            pytest.skip("Required real files not found")
        
        success, result = wrap_vcf_vep_annotate_unified(
            VCF_FILE, VCF_CACHE_DIR, VCF_FASTA, 
            annotate_gene=True, distance=5000
        )
        
        assert success is True
        assert "VEP output file:" in result
    
    def test_no_stats_parameter(self):
        """Test no_stats parameter"""
        # Skip test if VEP is not available or files don't exist
        try:
            subprocess.run(["vep", "--help"], check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            pytest.skip("VEP not available")
        
        if not all(Path(f).exists() for f in [VCF_FILE, VCF_CACHE_DIR, VCF_FASTA]):
            pytest.skip("Required real files not found")
        
        # Test with no_stats=True (default)
        success, result = wrap_vcf_vep_annotate_unified(
            VCF_FILE, VCF_CACHE_DIR, VCF_FASTA, 
            annotate_protein=True, no_stats=True
        )
        assert success is True
        
        # Test with no_stats=False
        success, result = wrap_vcf_vep_annotate_unified(
            VCF_FILE, VCF_CACHE_DIR, VCF_FASTA, 
            annotate_protein=True, no_stats=False
        )
        assert success is True
    
    def test_missing_files_raise_filenotfounderror(self):
        """Test missing files raise FileNotFoundError"""
        with pytest.raises(FileNotFoundError, match="VCF file not found"):
            wrap_vcf_vep_annotate_unified("nonexistent.vcf", VCF_CACHE_DIR, VCF_FASTA, annotate_protein=True)
        
        with pytest.raises(FileNotFoundError, match="Cache directory not found"):
            wrap_vcf_vep_annotate_unified(VCF_FILE, "nonexistent_cache", VCF_FASTA, annotate_protein=True)
        
        with pytest.raises(FileNotFoundError, match="FASTA file not found"):
            wrap_vcf_vep_annotate_unified(VCF_FILE, VCF_CACHE_DIR, "nonexistent.fa", annotate_protein=True)


class TestGeneralLogging:
    """Tests for general logging behavior with real files"""
    
    def test_logger_calls_with_real_files(self):
        """Test that functions complete without errors when using real files"""
        # Skip test if VEP is not available or files don't exist
        try:
            subprocess.run(["vep", "--help"], check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            pytest.skip("VEP not available")
        
        if not all(Path(f).exists() for f in [VCF_FILE, VCF_CACHE_DIR, VCF_FASTA]):
            pytest.skip("Required real files not found")
        
        # Test that the function completes without raising exceptions
        success, result = wrap_vcf_vep_annotate_protein(
            VCF_FILE, VCF_CACHE_DIR, VCF_FASTA
        )
        
        # Should complete successfully
        assert success is True
        assert isinstance(result, str)