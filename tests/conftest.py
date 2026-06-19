import numpy as np
import anndata as ad
import pytest


@pytest.fixture
def adata():
    """100-obs AnnData with three groups; seeded for reproducibility."""
    rng = np.random.default_rng(42)
    X = rng.standard_normal((100, 5))
    a = ad.AnnData(X=X)
    a.obs["Group"] = rng.choice(["test", "control", "test_2"], size=100)
    return a


@pytest.fixture
def linspace_adata():
    """Deterministic AnnData with linearly spaced features and two alternating groups."""
    data = ad.AnnData(np.linspace(0, 100, 1000).reshape(100, 10))
    data.obs["Group"] = ["test", "control"] * 50
    return data
