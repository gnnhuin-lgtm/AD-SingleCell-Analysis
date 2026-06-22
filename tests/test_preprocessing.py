"""Tests for preprocessing module."""

import pytest
import numpy as np
import anndata as ad
from src.preprocessing import QualityControl, Normalize, VariableGeneSelection


@pytest.fixture
def sample_adata():
    """Create sample AnnData object for testing."""
    X = np.random.randint(0, 10, size=(100, 200))
    obs = {"cell_type": np.random.choice(["A", "B", "C"], 100)}
    var_names = [f"Gene_{i}" for i in range(200)]
    var_names[0:10] = [f"MT-{i}" for i in range(10)]
    
    adata = ad.AnnData(X=X, obs=obs, var={"gene_names": var_names})
    adata.var_names = var_names
    
    return adata


def test_quality_control_initialization(sample_adata):
    """Test QualityControl initialization."""
    qc = QualityControl(sample_adata)
    assert qc.adata.shape == (100, 200)


def test_calculate_qc_metrics(sample_adata):
    """Test QC metrics calculation."""
    qc = QualityControl(sample_adata)
    qc.calculate_qc_metrics()
    
    assert "n_genes_by_counts" in sample_adata.obs.columns
    assert "total_counts" in sample_adata.obs.columns
    assert "pct_counts_mt" in sample_adata.obs.columns


def test_normalize_initialization(sample_adata):
    """Test Normalize initialization."""
    normalizer = Normalize(sample_adata)
    assert normalizer.adata.shape == (100, 200)


def test_variable_gene_selection(sample_adata):
    """Test HVG selection initialization."""
    hvg = VariableGeneSelection(sample_adata)
    assert hvg.adata.shape == (100, 200)


if __name__ == "__main__":
    pytest.main([__file__])
