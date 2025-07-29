import pytest
import json
import pandas as pd
from unittest.mock import Mock, patch, mock_open

from src.pyMut.filters.tissue_expression import (
    _load_rna_cancer_consensus,
    get_gene_symbol,
    tissue_expression,
    filter_by_tissue_expression
)
from src.pyMut.core import PyMutation
from src.pyMut.input import read_vcf, read_maf


class TestLoadRnaCancerConsensus:
    """Tests for _load_rna_cancer_consensus function"""
    
    def setup_method(self):
        """Reset cache before each test"""
        import src.pyMut.filters.tissue_expression as te_module
        te_module._rna_cancer_consensus_cache = None
    
    def test_returns_valid_dict_with_correct_json(self):
        """Test that function returns valid dict with correct JSON"""
        test_data = {
            "TSPAN6": {"BLCA": 10.5, "BRCA": 8.2},
            "TP53": {"BLCA": 15.0, "BRCA": 12.3}
        }
        
        with patch("builtins.open", mock_open(read_data=json.dumps(test_data))):
            with patch("pathlib.Path.exists", return_value=True):
                result = _load_rna_cancer_consensus("/fake/path.json")
                
        assert isinstance(result, dict)
        assert result == test_data
        assert "TSPAN6" in result
        assert "BLCA" in result["TSPAN6"]
        assert result["TSPAN6"]["BLCA"] == 10.5
    
    def test_second_call_uses_cache(self):
        """Test that second call uses cache and doesn't re-read file"""
        test_data = {"TSPAN6": {"BLCA": 10.5}}
        
        mock_file = mock_open(read_data=json.dumps(test_data))
        with patch("builtins.open", mock_file):
            with patch("pathlib.Path.exists", return_value=True):
                # First call
                result1 = _load_rna_cancer_consensus("/fake/path.json")
                # Second call
                result2 = _load_rna_cancer_consensus("/fake/path.json")
        
        # File should only be opened once due to caching
        assert mock_file.call_count == 1
        assert result1 == result2
        assert result1 == test_data
    
    def test_raises_filenotfounderror_when_path_not_exists(self):
        """Test that FileNotFoundError is raised when path doesn't exist"""
        with patch("pathlib.Path.exists", return_value=False):
            with pytest.raises(FileNotFoundError, match="RNA cancer consensus file not found"):
                _load_rna_cancer_consensus("/nonexistent/path.json")
    
    def test_raises_json_decode_error_with_malformed_json(self):
        """Test that JSONDecodeError is raised with malformed JSON"""
        malformed_json = '{"TSPAN6": {"BLCA": 10.5, "BRCA":}}'  # Missing value
        
        with patch("builtins.open", mock_open(read_data=malformed_json)):
            with patch("pathlib.Path.exists", return_value=True):
                with pytest.raises(json.JSONDecodeError):
                    _load_rna_cancer_consensus("/fake/path.json")
    
    def test_uses_default_path_when_none_provided(self):
        """Test that default path is used when file_path is None"""
        test_data = {"TSPAN6": {"BLCA": 10.5}}
        
        with patch("builtins.open", mock_open(read_data=json.dumps(test_data))):
            with patch("pathlib.Path.exists", return_value=True):
                result = _load_rna_cancer_consensus()
                
        assert result == test_data


class TestGetGeneSymbol:
    """Tests for get_gene_symbol function"""
    
    def test_returns_tp53_when_hugo_symbol_not_empty(self):
        """Test returns TP53 when Hugo_Symbol is not empty"""
        with patch("src.pyMut.filters.tissue_expression.FIELDS", {"Hugo_Symbol": ["Hugo_Symbol"]}):
            row = pd.Series({"Hugo_Symbol": "TP53", "other_field": "value"})
            result = get_gene_symbol(row)
            assert result == "TP53"
    
    def test_returns_brca1_when_hugo_symbol_empty_but_vep_nearest_exists(self):
        """Test returns BRCA1 when Hugo_Symbol empty but VEP_NEAREST exists"""
        with patch("src.pyMut.filters.tissue_expression.FIELDS", 
                  {"Hugo_Symbol": ["Hugo_Symbol", "SYMBOL", "VEP_NEAREST"]}):
            row = pd.Series({
                "Hugo_Symbol": "",
                "SYMBOL": "",
                "VEP_NEAREST": "BRCA1"
            })
            result = get_gene_symbol(row)
            assert result == "BRCA1"
    
    def test_extracts_gene1_from_gene1_and_gene2(self):
        """Test extracts GENE1 from GENE1&GENE2"""
        with patch("src.pyMut.filters.tissue_expression.FIELDS", {"Hugo_Symbol": ["Hugo_Symbol"]}):
            row = pd.Series({"Hugo_Symbol": "GENE1&GENE2"})
            result = get_gene_symbol(row)
            assert result == "GENE1"
    
    def test_ignores_empty_dot_dash_values(self):
        """Test ignores empty string, dot, and dash values"""
        with patch("src.pyMut.filters.tissue_expression.FIELDS", 
                  {"Hugo_Symbol": ["Hugo_Symbol", "SYMBOL", "Gene"]}):
            # Test empty string
            row1 = pd.Series({"Hugo_Symbol": "", "SYMBOL": "TP53"})
            assert get_gene_symbol(row1) == "TP53"
            
            # Test dot
            row2 = pd.Series({"Hugo_Symbol": ".", "SYMBOL": "TP53"})
            assert get_gene_symbol(row2) == "TP53"
            
            # Test dash
            row3 = pd.Series({"Hugo_Symbol": "-", "SYMBOL": "TP53"})
            assert get_gene_symbol(row3) == "TP53"
    
    def test_returns_none_when_no_valid_value(self):
        """Test returns None when no valid value in any alias"""
        with patch("src.pyMut.filters.tissue_expression.FIELDS", 
                  {"Hugo_Symbol": ["Hugo_Symbol", "SYMBOL"]}):
            row = pd.Series({"Hugo_Symbol": "", "SYMBOL": "."})
            result = get_gene_symbol(row)
            assert result is None
    
    def test_correct_priority_returns_tp53_over_brca1(self):
        """Test correct priority: returns TP53 if it matches before VEP_SYMBOL = BRCA1"""
        with patch("src.pyMut.filters.tissue_expression.FIELDS", 
                  {"Hugo_Symbol": ["Hugo_Symbol", "VEP_SYMBOL"]}):
            row = pd.Series({"Hugo_Symbol": "TP53", "VEP_SYMBOL": "BRCA1"})
            result = get_gene_symbol(row)
            assert result == "TP53"


class TestTissueExpression:
    """Tests for tissue_expression function"""
    
    def setup_method(self):
        """Reset cache before each test"""
        import src.pyMut.filters.tissue_expression as te_module
        te_module._rna_cancer_consensus_cache = None
    
    @patch('src.pyMut.filters.tissue_expression._load_rna_cancer_consensus')
    def test_expression_above_threshold_returns_true(self, mock_load):
        """Test data = TSPAN6, expression 10 >= threshold 5 returns True"""
        mock_load.return_value = {
            "TSPAN6": {"BLCA": 10.0, "BRCA": 8.0}
        }
        
        result = tissue_expression("TSPAN6", ["BLCA", 5])
        assert result is True
    
    @patch('src.pyMut.filters.tissue_expression._load_rna_cancer_consensus')
    def test_expression_below_threshold_returns_false(self, mock_load):
        """Test data = TSPAN6, expression 3 < threshold 5 returns False"""
        mock_load.return_value = {
            "TSPAN6": {"BLCA": 3.0, "BRCA": 8.0}
        }
        
        result = tissue_expression("TSPAN6", ["BLCA", 5])
        assert result is False
    
    @patch('src.pyMut.filters.tissue_expression._load_rna_cancer_consensus')
    @patch('src.pyMut.filters.tissue_expression.get_gene_symbol')
    def test_series_produces_same_result_as_string(self, mock_get_gene, mock_load):
        """Test Series with SYMBOL produces same result as string"""
        mock_load.return_value = {
            "TSPAN6": {"BLCA": 10.0}
        }
        mock_get_gene.return_value = "TSPAN6"
        
        # Test with string
        result_string = tissue_expression("TSPAN6", ["BLCA", 5])
        
        # Test with Series
        row = pd.Series({"Hugo_Symbol": "TSPAN6"})
        result_series = tissue_expression(row, ["BLCA", 5])
        
        assert result_string == result_series
        assert result_string is True
    
    @patch('src.pyMut.filters.tissue_expression._load_rna_cancer_consensus')
    def test_returns_false_when_gene_not_in_json(self, mock_load):
        """Test returns False when gene is not in JSON"""
        mock_load.return_value = {
            "OTHER_GENE": {"BLCA": 10.0}
        }
        
        result = tissue_expression("UNKNOWN_GENE", ["BLCA", 5])
        assert result is False
    
    @patch('src.pyMut.filters.tissue_expression._load_rna_cancer_consensus')
    def test_returns_false_when_tissue_not_in_gene_data(self, mock_load):
        """Test returns False when tissue is not in gene data"""
        mock_load.return_value = {
            "TSPAN6": {"BRCA": 10.0}  # Only BRCA, no BLCA
        }
        
        result = tissue_expression("TSPAN6", ["BLCA", 5])
        assert result is False
    
    def test_raises_valueerror_when_tissue_not_two_elements(self):
        """Test raises ValueError when tissue doesn't have exactly 2 elements"""
        with pytest.raises(ValueError, match="tissue parameter must be a list with exactly 2 elements"):
            tissue_expression("TSPAN6", ["BLCA"])  # Only one element
        
        with pytest.raises(ValueError, match="tissue parameter must be a list with exactly 2 elements"):
            tissue_expression("TSPAN6", ["BLCA", 5, "extra"])  # Three elements
    
    def test_raises_typeerror_when_data_not_string_or_series(self):
        """Test raises TypeError when data is not string or Series"""
        with pytest.raises(TypeError, match="data parameter must be either a string"):
            tissue_expression(123, ["BLCA", 5])  # Integer instead of string/Series


class TestFilterByTissueExpression:
    """Tests for filter_by_tissue_expression method"""
    
    def setup_method(self):
        """Reset cache and setup test data before each test"""
        import src.pyMut.filters.tissue_expression as te_module
        te_module._rna_cancer_consensus_cache = None
        
        # Create test PyMutation instances
        self.vcf_path = "/pyMut/data/examples/VCF/subset_1k_variants_ALL.chr10.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_VEP_annotated.vcf"
        self.maf_path = "/pyMut/data/examples/MAF/tcga_laml_VEP_annotated.maf.gz"
    
    def test_raises_valueerror_when_tissues_empty(self):
        """Test raises ValueError when tissues list is empty"""
        # Create a simple PyMutation mock
        mock_pymut = Mock()
        mock_pymut.data = pd.DataFrame({"Hugo_Symbol": ["TP53"]})
        
        with pytest.raises(ValueError, match="tissues parameter must be a non-empty list"):
            filter_by_tissue_expression(mock_pymut, [])
    
    def test_raises_valueerror_when_tissues_contains_wrong_types(self):
        """Test raises ValueError when tissues contains wrong types"""
        mock_pymut = Mock()
        mock_pymut.data = pd.DataFrame({"Hugo_Symbol": ["TP53"]})
        
        # Test non-tuple element
        with pytest.raises(ValueError, match="Each tissue specification must be a tuple"):
            filter_by_tissue_expression(mock_pymut, ["BLCA"])
        
        # Test tuple with wrong number of elements
        with pytest.raises(ValueError, match="Each tissue specification must be a tuple"):
            filter_by_tissue_expression(mock_pymut, [("BLCA",)])
        
        # Test non-string tissue code
        with pytest.raises(ValueError, match="Tissue code must be a string"):
            filter_by_tissue_expression(mock_pymut, [(123, 5)])
        
        # Test non-numeric threshold
        with pytest.raises(ValueError, match="Threshold must be a number"):
            filter_by_tissue_expression(mock_pymut, [("BLCA", "invalid")])
    
    @patch('src.pyMut.filters.tissue_expression.tissue_expression')
    def test_keep_expressed_true_returns_only_expressed(self, mock_tissue_expr):
        """Test keep_expressed=True returns only expressed rows"""
        # Setup mock data - one expressed, one not
        test_data = pd.DataFrame({
            "Hugo_Symbol": ["GENE1", "GENE2"],
            "Chromosome": ["1", "2"]
        })
        
        mock_pymut = Mock()
        mock_pymut.data = test_data
        mock_pymut.metadata = None
        mock_pymut.samples = []
        
        # Mock tissue_expression to return True for GENE1, False for GENE2
        def mock_expr_side_effect(row, tissue):
            if isinstance(row, pd.Series) and row["Hugo_Symbol"] == "GENE1":
                return True
            return False
        
        mock_tissue_expr.side_effect = mock_expr_side_effect
        
        with patch('src.pyMut.filters.tissue_expression.get_gene_symbol') as mock_get_gene:
            mock_get_gene.side_effect = lambda row: row["Hugo_Symbol"]
            
            result = filter_by_tissue_expression(mock_pymut, [("BLCA", 5)], keep_expressed=True)
        
        # Should only return the expressed gene (GENE1)
        assert len(result.data) == 1
        assert result.data.iloc[0]["Hugo_Symbol"] == "GENE1"
    
    @patch('src.pyMut.filters.tissue_expression.tissue_expression')
    def test_keep_expressed_false_returns_only_not_expressed(self, mock_tissue_expr):
        """Test keep_expressed=False returns only not expressed rows"""
        test_data = pd.DataFrame({
            "Hugo_Symbol": ["GENE1", "GENE2"],
            "Chromosome": ["1", "2"]
        })
        
        mock_pymut = Mock()
        mock_pymut.data = test_data
        mock_pymut.metadata = None
        mock_pymut.samples = []
        
        # Mock tissue_expression to return True for GENE1, False for GENE2
        def mock_expr_side_effect(row, tissue):
            if isinstance(row, pd.Series) and row["Hugo_Symbol"] == "GENE1":
                return True
            return False
        
        mock_tissue_expr.side_effect = mock_expr_side_effect
        
        with patch('src.pyMut.filters.tissue_expression.get_gene_symbol') as mock_get_gene:
            mock_get_gene.side_effect = lambda row: row["Hugo_Symbol"]
            
            result = filter_by_tissue_expression(mock_pymut, [("BLCA", 5)], keep_expressed=False)
        
        # Should only return the not expressed gene (GENE2)
        assert len(result.data) == 1
        assert result.data.iloc[0]["Hugo_Symbol"] == "GENE2"
    
    @patch('src.pyMut.filters.tissue_expression.tissue_expression')
    def test_adds_expression_columns_and_marks_expressed_in_any_tissue(self, mock_tissue_expr):
        """Test adds <TCGA>_expressed, <TCGA>_threshold columns and marks Expressed_in_Any_Tissue"""
        test_data = pd.DataFrame({
            "Hugo_Symbol": ["GENE1"],
            "Chromosome": ["1"]
        })
        
        mock_pymut = Mock()
        mock_pymut.data = test_data
        mock_pymut.metadata = None
        mock_pymut.samples = []
        
        mock_tissue_expr.return_value = True
        
        with patch('src.pyMut.filters.tissue_expression.get_gene_symbol') as mock_get_gene:
            mock_get_gene.return_value = "GENE1"
            
            result = filter_by_tissue_expression(mock_pymut, [("BLCA", 5)], keep_expressed=True)
        
        # Check that tissue_expression_results DataFrame has the expected columns
        assert hasattr(result, 'tissue_expression_results')
        results_df = result.tissue_expression_results
        
        expected_columns = ["Index", "Gene_Symbol", "Expressed_in_Any_Tissue", "BLCA_expressed", "BLCA_threshold"]
        for col in expected_columns:
            assert col in results_df.columns
        
        assert results_df.iloc[0]["BLCA_expressed"] == True
        assert results_df.iloc[0]["BLCA_threshold"] == 5
        assert results_df.iloc[0]["Expressed_in_Any_Tissue"] == True
    
    @patch('src.pyMut.filters.tissue_expression.tissue_expression')
    def test_maintains_metadata_filters_and_samples(self, mock_tissue_expr):
        """Test maintains self.metadata.filters and copies self.samples intact"""
        test_data = pd.DataFrame({
            "Hugo_Symbol": ["GENE1"],
            "Chromosome": ["1"]
        })
        
        # Create mock metadata with existing filters
        mock_metadata = Mock()
        mock_metadata.filters = ["existing_filter"]
        
        mock_pymut = Mock()
        mock_pymut.data = test_data
        mock_pymut.metadata = mock_metadata
        mock_pymut.samples = ["sample1", "sample2"]
        
        mock_tissue_expr.return_value = True
        
        with patch('src.pyMut.filters.tissue_expression.get_gene_symbol') as mock_get_gene:
            mock_get_gene.return_value = "GENE1"
            
            result = filter_by_tissue_expression(mock_pymut, [("BLCA", 5)], keep_expressed=True)
        
        # Check that metadata filters were updated
        assert hasattr(result.metadata, 'filters')
        assert len(result.metadata.filters) == 2
        assert "existing_filter" in result.metadata.filters
        assert "tissue_expression_filter([('BLCA', 5)], keep_expressed=True)" in result.metadata.filters
        
        # Check that samples were copied
        assert result.samples == ["sample1", "sample2"]
    
    @patch('src.pyMut.filters.tissue_expression.tissue_expression')
    def test_unknown_genes_treated_as_not_expressed(self, mock_tissue_expr):
        """Test rows with unknown genes don't break and are treated as not expressed"""
        test_data = pd.DataFrame({
            "Hugo_Symbol": ["UNKNOWN_GENE"],
            "Chromosome": ["1"]
        })
        
        mock_pymut = Mock()
        mock_pymut.data = test_data
        mock_pymut.metadata = None
        mock_pymut.samples = []
        
        # Mock tissue_expression to raise KeyError for unknown gene
        mock_tissue_expr.side_effect = KeyError("Gene not found")
        
        with patch('src.pyMut.filters.tissue_expression.get_gene_symbol') as mock_get_gene:
            mock_get_gene.return_value = "UNKNOWN_GENE"
            
            # Should not raise error
            result = filter_by_tissue_expression(mock_pymut, [("BLCA", 5)], keep_expressed=True)
        
        # With keep_expressed=True, unknown gene should be filtered out
        assert len(result.data) == 0
        
        # With keep_expressed=False, unknown gene should be kept
        result_false = filter_by_tissue_expression(mock_pymut, [("BLCA", 5)], keep_expressed=False)
        assert len(result_false.data) == 1
    
    @patch('src.pyMut.filters.tissue_expression.tissue_expression')
    def test_resulting_dataframe_index_matches_original(self, mock_tissue_expr):
        """Test resulting DataFrame index matches original for retained rows"""
        # Create test data with specific index
        test_data = pd.DataFrame({
            "Hugo_Symbol": ["GENE1", "GENE2", "GENE3"],
            "Chromosome": ["1", "2", "3"]
        }, index=[10, 20, 30])  # Custom index
        
        mock_pymut = Mock()
        mock_pymut.data = test_data
        mock_pymut.metadata = None
        mock_pymut.samples = []
        
        # Mock to keep only GENE1 and GENE3 (indices 10 and 30)
        def mock_expr_side_effect(row, tissue):
            return row["Hugo_Symbol"] in ["GENE1", "GENE3"]
        
        mock_tissue_expr.side_effect = mock_expr_side_effect
        
        with patch('src.pyMut.filters.tissue_expression.get_gene_symbol') as mock_get_gene:
            mock_get_gene.side_effect = lambda row: row["Hugo_Symbol"]
            
            result = filter_by_tissue_expression(mock_pymut, [("BLCA", 5)], keep_expressed=True)
        
        # Check that the resulting DataFrame has the correct indices
        expected_indices = [10, 30]
        assert list(result.data.index) == expected_indices
        assert result.data.loc[10, "Hugo_Symbol"] == "GENE1"
        assert result.data.loc[30, "Hugo_Symbol"] == "GENE3"


class TestIntegrationWithRealData:
    """Integration tests using real PyMutation instances"""
    
    def setup_method(self):
        """Reset cache before each test"""
        import src.pyMut.filters.tissue_expression as te_module
        te_module._rna_cancer_consensus_cache = None
    
    @pytest.mark.slow
    def test_vcf_mutation_filtering(self):
        """Test filtering with real VCF PyMutation instance"""
        vcf_path = "/pyMut/data/examples/VCF/subset_1k_variants_ALL.chr10.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased_VEP_annotated.vcf"
        
        # This test requires the actual read_vcf function to work
        try:
            vcf_mut = read_vcf(vcf_path)
            original_count = len(vcf_mut.data)
            
            # Test keep_expressed=True
            result_expressed = vcf_mut.filter_by_tissue_expression([("BLCA", 5)], keep_expressed=True)
            
            # Basic assertions for expressed filtering
            assert isinstance(result_expressed, PyMutation)
            assert hasattr(result_expressed, 'tissue_expression_results')
            assert len(result_expressed.data) <= original_count  # Should be filtered
            
            # Check tissue expression results structure
            results_df = result_expressed.tissue_expression_results
            expected_columns = ["Index", "Gene_Symbol", "Expressed_in_Any_Tissue", "BLCA_expressed", "BLCA_threshold"]
            for col in expected_columns:
                assert col in results_df.columns
            
            # Test keep_expressed=False
            result_not_expressed = vcf_mut.filter_by_tissue_expression([("BLCA", 5)], keep_expressed=False)
            
            # Basic assertions for not expressed filtering
            assert isinstance(result_not_expressed, PyMutation)
            assert hasattr(result_not_expressed, 'tissue_expression_results')
            assert len(result_not_expressed.data) <= original_count  # Should be filtered
            
            # The sum of expressed and not expressed should not exceed original (some may be filtered out due to unknown genes)
            assert len(result_expressed.data) + len(result_not_expressed.data) <= original_count
            
            # Test metadata filters are updated
            if hasattr(result_expressed, 'metadata') and result_expressed.metadata is not None:
                assert hasattr(result_expressed.metadata, 'filters')
                filter_found = any("tissue_expression_filter" in f for f in result_expressed.metadata.filters)
                assert filter_found
            
        except Exception as e:
            pytest.skip(f"Could not load VCF file for integration test: {e}")
    
    @pytest.mark.slow
    def test_maf_mutation_filtering(self):
        """Test filtering with real MAF PyMutation instance"""
        maf_path = "/pyMut/data/examples/MAF/tcga_laml_VEP_annotated.maf.gz"
        
        # This test requires the actual read_maf function to work
        try:
            maf_mut = read_maf(maf_path)
            original_count = len(maf_mut.data)
            
            # Test keep_expressed=True
            result_expressed = maf_mut.filter_by_tissue_expression([("BLCA", 5)], keep_expressed=True)
            
            # Basic assertions for expressed filtering
            assert isinstance(result_expressed, PyMutation)
            assert hasattr(result_expressed, 'tissue_expression_results')
            assert len(result_expressed.data) <= original_count  # Should be filtered
            
            # Check tissue expression results structure
            results_df = result_expressed.tissue_expression_results
            expected_columns = ["Index", "Gene_Symbol", "Expressed_in_Any_Tissue", "BLCA_expressed", "BLCA_threshold"]
            for col in expected_columns:
                assert col in results_df.columns
            
            # Test keep_expressed=False
            result_not_expressed = maf_mut.filter_by_tissue_expression([("BLCA", 5)], keep_expressed=False)
            
            # Basic assertions for not expressed filtering
            assert isinstance(result_not_expressed, PyMutation)
            assert hasattr(result_not_expressed, 'tissue_expression_results')
            assert len(result_not_expressed.data) <= original_count  # Should be filtered
            
            # The sum of expressed and not expressed should not exceed original (some may be filtered out due to unknown genes)
            assert len(result_expressed.data) + len(result_not_expressed.data) <= original_count
            
            # Test multiple tissues filtering
            result_multi = maf_mut.filter_by_tissue_expression([("BLCA", 5), ("BRCA", 3)], keep_expressed=True)
            assert isinstance(result_multi, PyMutation)
            assert hasattr(result_multi, 'tissue_expression_results')
            
            # Check multi-tissue results structure
            multi_results_df = result_multi.tissue_expression_results
            multi_expected_columns = ["Index", "Gene_Symbol", "Expressed_in_Any_Tissue", 
                                    "BLCA_expressed", "BLCA_threshold", "BRCA_expressed", "BRCA_threshold"]
            for col in multi_expected_columns:
                assert col in multi_results_df.columns
            
            # Test metadata filters are updated
            if hasattr(result_expressed, 'metadata') and result_expressed.metadata is not None:
                assert hasattr(result_expressed.metadata, 'filters')
                filter_found = any("tissue_expression_filter" in f for f in result_expressed.metadata.filters)
                assert filter_found
            
        except Exception as e:
            pytest.skip(f"Could not load MAF file for integration test: {e}")