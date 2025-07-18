import pytest
import subprocess
import tempfile
import gzip
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock, call
import pandas as pd

from src.pyMut.annotate.vep_annotate import (
    _extract_assembly_and_version_from_cache,
    _get_case_insensitive_column,
    _maf_to_region,
    wrap_maf_vep_annotate_protein,
    wrap_vcf_vep_annotate_protein,
    wrap_vcf_vep_annotate_gene
)


class TestExtractAssemblyAndVersionFromCache:
    """Tests for _extract_assembly_and_version_from_cache function"""
    
    def test_valid_cache_name_grch38(self):
        """Test valid cache name returns correct assembly and version"""
        cache_dir = Path("homo_sapiens_vep_114_GRCh38")
        assembly, version = _extract_assembly_and_version_from_cache(cache_dir)
        assert assembly == "GRCh38"
        assert version == "114"
    
    def test_valid_cache_name_with_dot(self):
        """Test cache name with dot in assembly (GRCh38.p13) extracts correctly"""
        cache_dir = Path("homo_sapiens_vep_114_GRCh38.p13")
        assembly, version = _extract_assembly_and_version_from_cache(cache_dir)
        assert assembly == "GRCh38.p13"
        assert version == "114"
    
    def test_invalid_cache_name_raises_valueerror(self):
        """Test invalid cache name raises ValueError"""
        cache_dir = Path("invalid_cache_name")
        with pytest.raises(ValueError, match="Cache directory name 'invalid_cache_name' doesn't match expected format"):
            _extract_assembly_and_version_from_cache(cache_dir)


class TestGetCaseInsensitiveColumn:
    """Tests for _get_case_insensitive_column function"""
    
    def test_case_insensitive_match(self):
        """Test case-insensitive column matching"""
        columns = ["Chromosome", "Start_Position", "End_Position"]
        result = _get_case_insensitive_column(columns, "chromosome")
        assert result == "Chromosome"
        
        result = _get_case_insensitive_column(columns, "START_POSITION")
        assert result == "Start_Position"
    
    def test_missing_column_raises_keyerror(self):
        """Test missing column raises KeyError"""
        columns = ["Chromosome", "Start_Position"]
        with pytest.raises(KeyError, match="Column 'Missing_Column' not found in MAF file"):
            _get_case_insensitive_column(columns, "Missing_Column")


class TestMafToRegion:
    """Tests for _maf_to_region function"""
    
    def test_nonexistent_file_returns_false(self):
        """Test nonexistent file returns (False, '')"""
        success, path = _maf_to_region("nonexistent_file.maf")
        assert success is False
        assert path == ""
    
    def test_maf_with_strand_column(self):
        """Test MAF file with Strand column produces correct region format"""
        # Create temporary MAF file with Strand column
        with tempfile.NamedTemporaryFile(mode='w', suffix='.maf', delete=False) as f:
            f.write("Chromosome\tStart_Position\tEnd_Position\tTumor_Seq_Allele2\tStrand\n")
            f.write("chr1\t100\t100\tA\t+\n")
            f.write("2\t200\t200\tG\t-\n")
            temp_maf = f.name
        
        try:
            success, region_path = _maf_to_region(temp_maf)
            assert success is True
            
            # Check region file content
            with open(region_path, 'r') as f:
                lines = f.readlines()
            
            assert len(lines) == 2
            assert lines[0].strip() == "chr1:100-100:1/A"
            assert lines[1].strip() == "chr2:200-200:-1/G"
            
            # Clean up
            Path(region_path).unlink()
        finally:
            Path(temp_maf).unlink()
    
    def test_maf_without_strand_column(self):
        """Test MAF file without Strand column uses default '+' strand"""
        # Create temporary MAF file without Strand column
        with tempfile.NamedTemporaryFile(mode='w', suffix='.maf', delete=False) as f:
            f.write("Chromosome\tStart_Position\tEnd_Position\tTumor_Seq_Allele2\n")
            f.write("chr1\t100\t100\tA\n")
            temp_maf = f.name
        
        try:
            with patch('src.pyMut.annotate.vep_annotate.logger') as mock_logger:
                success, region_path = _maf_to_region(temp_maf)
                assert success is True
                
                # Check warning was logged
                mock_logger.warning.assert_called_with("Strand column not found in MAF file, using default value '+'")
                
                # Check region file content
                with open(region_path, 'r') as f:
                    lines = f.readlines()
                
                assert len(lines) == 1
                assert lines[0].strip() == "chr1:100-100:1/A"
                
                # Clean up
                Path(region_path).unlink()
        finally:
            Path(temp_maf).unlink()
    
    def test_maf_missing_required_columns(self):
        """Test MAF file missing required columns returns (False, path)"""
        # Create temporary MAF file with missing columns
        with tempfile.NamedTemporaryFile(mode='w', suffix='.maf', delete=False) as f:
            f.write("OnlyColumn\n")
            f.write("value\n")
            temp_maf = f.name
        
        try:
            success, region_path = _maf_to_region(temp_maf)
            assert success is False
            assert region_path == temp_maf.replace('.maf', '.region')
        finally:
            Path(temp_maf).unlink()
    
    def test_empty_maf_file(self):
        """Test empty MAF file returns (False, path)"""
        # Create empty MAF file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.maf', delete=False) as f:
            pass  # Empty file
        temp_maf = f.name
        
        try:
            success, region_path = _maf_to_region(temp_maf)
            assert success is False
            assert region_path == temp_maf.replace('.maf', '.region')
        finally:
            Path(temp_maf).unlink()
    
    def test_chromosome_normalization(self):
        """Test chromosome names are normalized correctly"""
        # Create temporary MAF file with various chromosome formats
        with tempfile.NamedTemporaryFile(mode='w', suffix='.maf', delete=False) as f:
            f.write("Chromosome\tStart_Position\tEnd_Position\tTumor_Seq_Allele2\n")
            f.write("1\t100\t100\tA\n")
            f.write("chr2\t200\t200\tG\n")
            f.write("X\t300\t300\tT\n")
            temp_maf = f.name
        
        try:
            success, region_path = _maf_to_region(temp_maf)
            assert success is True
            
            # Check region file content
            with open(region_path, 'r') as f:
                lines = f.readlines()
            
            assert len(lines) == 3
            assert lines[0].strip() == "chr1:100-100:1/A"
            assert lines[1].strip() == "chr2:200-200:1/G"
            assert lines[2].strip() == "chrX:300-300:1/T"
            
            # Clean up
            Path(region_path).unlink()
        finally:
            Path(temp_maf).unlink()


class TestWrapMafVepAnnotateProtein:
    """Tests for wrap_maf_vep_annotate_protein function"""
    
    @patch('src.pyMut.annotate.vep_annotate.merge_maf_with_vep_annotations')
    @patch('src.pyMut.annotate.vep_annotate.subprocess.run')
    @patch('src.pyMut.annotate.vep_annotate._maf_to_region')
    @patch('src.pyMut.annotate.vep_annotate.logger')
    def test_successful_annotation_with_mocks(self, mock_logger, mock_maf_to_region, mock_subprocess, mock_merge):
        """Test successful MAF VEP annotation with mocks"""
        # Setup mocks
        mock_maf_to_region.return_value = (True, "/tmp/test.region")
        mock_subprocess.return_value = Mock(stderr="", returncode=0)
        mock_merge.return_value = (Mock(), "/tmp/merged.maf")
        
        # Create temporary files for testing
        with tempfile.NamedTemporaryFile(suffix='.maf') as maf_file, \
             tempfile.TemporaryDirectory() as cache_dir, \
             tempfile.NamedTemporaryFile(suffix='.fa') as fasta_file:
            
            # Create cache directory structure
            cache_path = Path(cache_dir) / "homo_sapiens_vep_114_GRCh38"
            cache_path.mkdir(parents=True)
            
            success, result = wrap_maf_vep_annotate_protein(
                maf_file.name, cache_path, fasta_file.name
            )
            
            assert success is True
            assert "VEP folder:" in result
            assert "Merged file:" in result
            
            # Verify subprocess was called with correct arguments
            mock_subprocess.assert_called_once()
            call_args = mock_subprocess.call_args[0][0]
            assert "--format" in call_args
            assert "region" in call_args
            assert "--protein" in call_args
            assert "--uniprot" in call_args
            assert "--domains" in call_args
            assert "--symbol" in call_args
            assert "--pick" in call_args
            assert "--no_stats" not in call_args  # no_stats=True by default
    
    @patch('src.pyMut.annotate.vep_annotate._maf_to_region')
    @patch('src.pyMut.annotate.vep_annotate.logger')
    def test_maf_to_region_failure(self, mock_logger, mock_maf_to_region):
        """Test failure in _maf_to_region returns False"""
        mock_maf_to_region.return_value = (False, "")
        
        with tempfile.NamedTemporaryFile(suffix='.maf') as maf_file, \
             tempfile.TemporaryDirectory() as cache_dir, \
             tempfile.NamedTemporaryFile(suffix='.fa') as fasta_file:
            
            cache_path = Path(cache_dir) / "homo_sapiens_vep_114_GRCh38"
            cache_path.mkdir(parents=True)
            
            success, result = wrap_maf_vep_annotate_protein(
                maf_file.name, cache_path, fasta_file.name
            )
            
            assert success is False
            assert result == ""
    
    @patch('src.pyMut.annotate.vep_annotate.merge_maf_with_vep_annotations')
    @patch('src.pyMut.annotate.vep_annotate.subprocess.run')
    @patch('src.pyMut.annotate.vep_annotate._maf_to_region')
    @patch('src.pyMut.annotate.vep_annotate.logger')
    def test_subprocess_error(self, mock_logger, mock_maf_to_region, mock_subprocess, mock_merge):
        """Test subprocess error returns False"""
        mock_maf_to_region.return_value = (True, "/tmp/test.region")
        mock_subprocess.side_effect = subprocess.CalledProcessError(1, "vep")
        
        with tempfile.NamedTemporaryFile(suffix='.maf') as maf_file, \
             tempfile.TemporaryDirectory() as cache_dir, \
             tempfile.NamedTemporaryFile(suffix='.fa') as fasta_file:
            
            cache_path = Path(cache_dir) / "homo_sapiens_vep_114_GRCh38"
            cache_path.mkdir(parents=True)
            
            success, result = wrap_maf_vep_annotate_protein(
                maf_file.name, cache_path, fasta_file.name
            )
            
            assert success is False
    
    @patch('src.pyMut.annotate.vep_annotate.merge_maf_with_vep_annotations')
    @patch('src.pyMut.annotate.vep_annotate.subprocess.run')
    @patch('src.pyMut.annotate.vep_annotate._maf_to_region')
    @patch('src.pyMut.annotate.vep_annotate.logger')
    def test_merge_failure_returns_true_with_error_message(self, mock_logger, mock_maf_to_region, mock_subprocess, mock_merge):
        """Test merge failure returns True with error message"""
        mock_maf_to_region.return_value = (True, "/tmp/test.region")
        mock_subprocess.return_value = Mock(stderr="", returncode=0)
        mock_merge.side_effect = Exception("Merge failed")
        
        with tempfile.NamedTemporaryFile(suffix='.maf') as maf_file, \
             tempfile.TemporaryDirectory() as cache_dir, \
             tempfile.NamedTemporaryFile(suffix='.fa') as fasta_file:
            
            cache_path = Path(cache_dir) / "homo_sapiens_vep_114_GRCh38"
            cache_path.mkdir(parents=True)
            
            success, result = wrap_maf_vep_annotate_protein(
                maf_file.name, cache_path, fasta_file.name
            )
            
            assert success is True
            assert "Merge failed" in result


class TestWrapVcfVepAnnotateProtein:
    """Tests for wrap_vcf_vep_annotate_protein function"""
    
    @patch('src.pyMut.annotate.vep_annotate.subprocess.run')
    @patch('src.pyMut.annotate.vep_annotate.logger')
    def test_successful_annotation_with_vcf_flag(self, mock_logger, mock_subprocess):
        """Test VCF annotation includes --vcf flag and protein flags"""
        mock_subprocess.return_value = Mock(stderr="", returncode=0)
        
        with tempfile.NamedTemporaryFile(suffix='.vcf') as vcf_file, \
             tempfile.TemporaryDirectory() as cache_dir, \
             tempfile.NamedTemporaryFile(suffix='.fa') as fasta_file:
            
            cache_path = Path(cache_dir) / "homo_sapiens_vep_114_GRCh38"
            cache_path.mkdir(parents=True)
            
            success, result = wrap_vcf_vep_annotate_protein(
                vcf_file.name, cache_path, fasta_file.name
            )
            
            assert success is True
            assert "VEP output file:" in result
            
            # Verify subprocess was called with correct arguments
            mock_subprocess.assert_called_once()
            call_args = mock_subprocess.call_args[0][0]
            assert "--vcf" in call_args
            assert "--protein" in call_args
            assert "--uniprot" in call_args
            assert "--domains" in call_args
            assert "--symbol" in call_args
            assert "--pick" in call_args
    
    @patch('src.pyMut.annotate.vep_annotate.subprocess.run')
    @patch('src.pyMut.annotate.vep_annotate.logger')
    def test_no_stats_true_omits_flag(self, mock_logger, mock_subprocess):
        """Test no_stats=True omits --no_stats flag"""
        mock_subprocess.return_value = Mock(stderr="", returncode=0)
        
        with tempfile.NamedTemporaryFile(suffix='.vcf') as vcf_file, \
             tempfile.TemporaryDirectory() as cache_dir, \
             tempfile.NamedTemporaryFile(suffix='.fa') as fasta_file:
            
            cache_path = Path(cache_dir) / "homo_sapiens_vep_114_GRCh38"
            cache_path.mkdir(parents=True)
            
            wrap_vcf_vep_annotate_protein(
                vcf_file.name, cache_path, fasta_file.name, no_stats=True
            )
            
            call_args = mock_subprocess.call_args[0][0]
            assert "--no_stats" not in call_args
    
    @patch('src.pyMut.annotate.vep_annotate.subprocess.run')
    @patch('src.pyMut.annotate.vep_annotate.logger')
    def test_no_stats_false_includes_flag(self, mock_logger, mock_subprocess):
        """Test no_stats=False includes --no_stats flag"""
        mock_subprocess.return_value = Mock(stderr="", returncode=0)
        
        with tempfile.NamedTemporaryFile(suffix='.vcf') as vcf_file, \
             tempfile.TemporaryDirectory() as cache_dir, \
             tempfile.NamedTemporaryFile(suffix='.fa') as fasta_file:
            
            cache_path = Path(cache_dir) / "homo_sapiens_vep_114_GRCh38"
            cache_path.mkdir(parents=True)
            
            wrap_vcf_vep_annotate_protein(
                vcf_file.name, cache_path, fasta_file.name, no_stats=False
            )
            
            call_args = mock_subprocess.call_args[0][0]
            assert "--no_stats" in call_args
    
    def test_missing_files_raise_filenotfounderror(self):
        """Test missing files raise FileNotFoundError"""
        with pytest.raises(FileNotFoundError, match="VCF file not found"):
            wrap_vcf_vep_annotate_protein("nonexistent.vcf", "cache", "fasta.fa")
        
        with tempfile.NamedTemporaryFile(suffix='.vcf') as vcf_file:
            with pytest.raises(FileNotFoundError, match="Cache directory not found"):
                wrap_vcf_vep_annotate_protein(vcf_file.name, "nonexistent_cache", "fasta.fa")
            
            with tempfile.TemporaryDirectory() as cache_dir:
                with pytest.raises(FileNotFoundError, match="FASTA file not found"):
                    wrap_vcf_vep_annotate_protein(vcf_file.name, cache_dir, "nonexistent.fa")


class TestWrapVcfVepAnnotateGene:
    """Tests for wrap_vcf_vep_annotate_gene function"""
    
    @patch('src.pyMut.annotate.vep_annotate.subprocess.run')
    @patch('src.pyMut.annotate.vep_annotate.logger')
    def test_with_distance_includes_nearest_flags(self, mock_logger, mock_subprocess):
        """Test distance parameter includes --nearest symbol --distance flags"""
        mock_subprocess.return_value = Mock(stderr="", returncode=0)
        
        with tempfile.NamedTemporaryFile(suffix='.vcf') as vcf_file, \
             tempfile.TemporaryDirectory() as cache_dir, \
             tempfile.NamedTemporaryFile(suffix='.fa') as fasta_file:
            
            cache_path = Path(cache_dir) / "homo_sapiens_vep_114_GRCh38"
            cache_path.mkdir(parents=True)
            
            wrap_vcf_vep_annotate_gene(
                vcf_file.name, cache_path, fasta_file.name, distance=10000
            )
            
            call_args = mock_subprocess.call_args[0][0]
            assert "--nearest" in call_args
            assert "symbol" in call_args
            assert "--distance" in call_args
            assert "10000" in call_args
    
    @patch('src.pyMut.annotate.vep_annotate.subprocess.run')
    @patch('src.pyMut.annotate.vep_annotate.logger')
    def test_without_distance_omits_nearest_flags(self, mock_logger, mock_subprocess):
        """Test without distance parameter omits --nearest and --distance flags"""
        mock_subprocess.return_value = Mock(stderr="", returncode=0)
        
        with tempfile.NamedTemporaryFile(suffix='.vcf') as vcf_file, \
             tempfile.TemporaryDirectory() as cache_dir, \
             tempfile.NamedTemporaryFile(suffix='.fa') as fasta_file:
            
            cache_path = Path(cache_dir) / "homo_sapiens_vep_114_GRCh38"
            cache_path.mkdir(parents=True)
            
            wrap_vcf_vep_annotate_gene(
                vcf_file.name, cache_path, fasta_file.name
            )
            
            call_args = mock_subprocess.call_args[0][0]
            assert "--nearest" not in call_args
            assert "--distance" not in call_args
    
    @patch('src.pyMut.annotate.vep_annotate.subprocess.run')
    @patch('src.pyMut.annotate.vep_annotate.logger')
    def test_no_stats_control(self, mock_logger, mock_subprocess):
        """Test no_stats parameter control"""
        mock_subprocess.return_value = Mock(stderr="", returncode=0)
        
        with tempfile.NamedTemporaryFile(suffix='.vcf') as vcf_file, \
             tempfile.TemporaryDirectory() as cache_dir, \
             tempfile.NamedTemporaryFile(suffix='.fa') as fasta_file:
            
            cache_path = Path(cache_dir) / "homo_sapiens_vep_114_GRCh38"
            cache_path.mkdir(parents=True)
            
            # Test no_stats=True (default)
            wrap_vcf_vep_annotate_gene(
                vcf_file.name, cache_path, fasta_file.name, no_stats=True
            )
            call_args = mock_subprocess.call_args[0][0]
            assert "--no_stats" not in call_args
            
            # Test no_stats=False
            wrap_vcf_vep_annotate_gene(
                vcf_file.name, cache_path, fasta_file.name, no_stats=False
            )
            call_args = mock_subprocess.call_args[0][0]
            assert "--no_stats" in call_args
    
    def test_missing_files_raise_filenotfounderror(self):
        """Test missing files raise FileNotFoundError"""
        with pytest.raises(FileNotFoundError, match="VCF file not found"):
            wrap_vcf_vep_annotate_gene("nonexistent.vcf", "cache", "fasta.fa")
        
        with tempfile.NamedTemporaryFile(suffix='.vcf') as vcf_file:
            with pytest.raises(FileNotFoundError, match="Cache directory not found"):
                wrap_vcf_vep_annotate_gene(vcf_file.name, "nonexistent_cache", "fasta.fa")
            
            with tempfile.TemporaryDirectory() as cache_dir:
                with pytest.raises(FileNotFoundError, match="FASTA file not found"):
                    wrap_vcf_vep_annotate_gene(vcf_file.name, cache_dir, "nonexistent.fa")


class TestGeneralLogging:
    """Tests for general logging behavior"""
    
    @patch('src.pyMut.annotate.vep_annotate.subprocess.run')
    @patch('src.pyMut.annotate.vep_annotate.logger')
    def test_logger_info_calls(self, mock_logger, mock_subprocess):
        """Test that logger.info is called appropriately"""
        mock_subprocess.return_value = Mock(stderr="", returncode=0)
        
        with tempfile.NamedTemporaryFile(suffix='.vcf') as vcf_file, \
             tempfile.TemporaryDirectory() as cache_dir, \
             tempfile.NamedTemporaryFile(suffix='.fa') as fasta_file:
            
            cache_path = Path(cache_dir) / "homo_sapiens_vep_114_GRCh38"
            cache_path.mkdir(parents=True)
            
            wrap_vcf_vep_annotate_protein(
                vcf_file.name, cache_path, fasta_file.name
            )
            
            # Verify logger.info was called
            assert mock_logger.info.called
    
    @patch('src.pyMut.annotate.vep_annotate.subprocess.run')
    @patch('src.pyMut.annotate.vep_annotate.logger')
    def test_logger_error_on_subprocess_failure(self, mock_logger, mock_subprocess):
        """Test that logger.error is called on subprocess failure"""
        mock_subprocess.side_effect = subprocess.CalledProcessError(1, "vep")
        
        with tempfile.NamedTemporaryFile(suffix='.vcf') as vcf_file, \
             tempfile.TemporaryDirectory() as cache_dir, \
             tempfile.NamedTemporaryFile(suffix='.fa') as fasta_file:
            
            cache_path = Path(cache_dir) / "homo_sapiens_vep_114_GRCh38"
            cache_path.mkdir(parents=True)
            
            wrap_vcf_vep_annotate_protein(
                vcf_file.name, cache_path, fasta_file.name
            )
            
            # Verify logger.error was called
            assert mock_logger.error.called
    
    @patch('src.pyMut.annotate.vep_annotate.subprocess.run')
    @patch('src.pyMut.annotate.vep_annotate.logger')
    def test_logger_warning_on_stderr(self, mock_logger, mock_subprocess):
        """Test that logger.warning is called when subprocess has stderr output"""
        mock_subprocess.return_value = Mock(stderr="Some warning message", returncode=0)
        
        with tempfile.NamedTemporaryFile(suffix='.vcf') as vcf_file, \
             tempfile.TemporaryDirectory() as cache_dir, \
             tempfile.NamedTemporaryFile(suffix='.fa') as fasta_file:
            
            cache_path = Path(cache_dir) / "homo_sapiens_vep_114_GRCh38"
            cache_path.mkdir(parents=True)
            
            wrap_vcf_vep_annotate_protein(
                vcf_file.name, cache_path, fasta_file.name
            )
            
            # Verify logger.warning was called
            assert mock_logger.warning.called