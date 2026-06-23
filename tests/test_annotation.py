"""Tests for annotation module."""

import pytest
import numpy as np
import anndata as ad
from src.annotation import CellTypeAnnotation, ADGeneSetAnalyzer, BRAIN_CELL_MARKERS


@pytest.fixture
def sample_adata():
    """Create sample AnnData for annotation testing."""
    n_cells, n_genes = 50, 100
    rng = np.random.default_rng(42)

    X = rng.poisson(lam=2, size=(n_cells, n_genes)).astype(float)

    var_names = list(BRAIN_CELL_MARKERS.keys())[:5]
    # Ensure some marker genes exist
    marker_subset = []
    for _, markers in BRAIN_CELL_MARKERS.items():
        marker_subset.extend(markers[:3])
    var_names = var_names + marker_subset

    # Pad with random gene names
    extra = n_genes - len(var_names)
    if extra > 0:
        var_names += [f"Gene_{i}" for i in range(extra)]

    adata = ad.AnnData(
        X=X,
        obs={"cell_barcode": [f"Cell_{i}" for i in range(n_cells)]},
        var={"gene_symbol": var_names},
    )
    adata.var_names = var_names

    return adata


def test_cell_type_annotation_initialization(sample_adata):
    """Test CellTypeAnnotation initialization."""
    annotator = CellTypeAnnotation(sample_adata)
    assert annotator.adata.shape == sample_adata.shape
    assert annotator.marker_dict is not None


def test_scoring_cell_types(sample_adata):
    """Test cell type scoring."""
    annotator = CellTypeAnnotation(sample_adata, min_score=0.0)
    scores = annotator.score_cell_types()
    assert isinstance(scores, type(pd.DataFrame()))
    assert "cell_type" in sample_adata.obs.columns
    assert "cell_type_score" in sample_adata.obs.columns


def test_ad_gene_set_analyzer_initialization(sample_adata):
    """Test ADGeneSetAnalyzer initialization."""
    analyzer = ADGeneSetAnalyzer(sample_adata)
    assert analyzer.adata is sample_adata
    assert len(analyzer.gene_sets) > 0


def test_ad_gene_set_scoring(sample_adata):
    """Test AD gene set scoring."""
    analyzer = ADGeneSetAnalyzer(sample_adata)
    scores = analyzer.score_gene_sets()
    assert isinstance(scores, type(pd.DataFrame()))


def test_ad_gene_set_summary(sample_adata):
    """Test gene set summary."""
    analyzer = ADGeneSetAnalyzer(sample_adata)
    analyzer.score_gene_sets()
    summary = analyzer.get_gene_set_summary()
    assert isinstance(summary, type(pd.DataFrame()))
