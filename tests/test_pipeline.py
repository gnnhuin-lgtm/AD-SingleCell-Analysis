"""Tests for the pipeline module."""

import os
import tempfile
import pytest
import numpy as np
import pandas as pd
import anndata as ad
from pathlib import Path

from src.pipeline import ADSingleCellPipeline


@pytest.fixture
def sample_h5ad():
    """Create a small test h5ad file."""
    X = np.random.default_rng(42).poisson(lam=2, size=(30, 50)).astype(float)
    var_names = [f"Gene_{i}" for i in range(50)]
    obs = pd.DataFrame({"condition": ["AD"] * 15 + ["Control"] * 15})
    adata = ad.AnnData(X=X, obs=obs, var={"gene_names": var_names})
    adata.var_names = var_names

    tmpdir = tempfile.mkdtemp()
    filepath = os.path.join(tmpdir, "test_data.h5ad")
    adata.write(filepath)

    return filepath, tmpdir


def test_pipeline_initialization():
    """Test pipeline initialization."""
    with tempfile.TemporaryDirectory() as tmpdir:
        pipeline = ADSingleCellPipeline(output_dir=tmpdir)
        assert pipeline.adata is None
        assert pipeline.output_dir == Path(tmpdir)


def test_pipeline_load_data(sample_h5ad):
    """Test data loading."""
    filepath, tmpdir = sample_h5ad
    pipeline = ADSingleCellPipeline(output_dir=tmpdir)
    pipeline.load_data(filepath)
    assert pipeline.adata is not None
    assert isinstance(pipeline.adata, ad.AnnData)
