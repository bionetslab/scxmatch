import pytest
import scanpy as sc

from scxmatch import test, estimate_peak_RAM_GB


# ── input validation ──────────────────────────────────────────────────────────

def test_invalid_group_by(adata):
    with pytest.raises(KeyError):
        test(adata, group_by="foo", test_group="test", reference="control", k=10)


def test_single_test_group(adata):
    test(adata, group_by="Group", test_group="test", reference="control", k=10)


def test_list_test_group(adata):
    test(adata, group_by="Group", test_group=["test", "test_2"], reference="control", k=10)


def test_invalid_test_group(adata):
    with pytest.raises(ValueError):
        test(adata, group_by="Group", test_group="foo", reference="control", k=10)


def test_single_reference(adata):
    test(adata, group_by="Group", test_group="test", reference="control", k=10)


def test_list_reference(adata):
    test(adata, group_by="Group", test_group="test", reference=["control", "test_2"], k=10)


def test_invalid_reference(adata):
    with pytest.raises(ValueError):
        test(adata, group_by="Group", test_group="test", reference="foo", k=10)


@pytest.mark.parametrize("metric", ["euclidean", "sqeuclidean", "jaccard"])
def test_valid_metric(adata, metric):
    test(adata, group_by="Group", test_group="test", reference="control", metric=metric, k=10)


def test_invalid_metric(adata):
    with pytest.raises(ValueError):
        test(adata, group_by="Group", test_group="test", reference="control", metric="foo", k=10)


@pytest.mark.parametrize("rank", [False, True])
def test_rank_parameter(adata, rank):
    test(adata, group_by="Group", test_group="test", reference="control", rank=rank, k=10)


@pytest.mark.parametrize("k", [10, None])
def test_k_parameter(adata, k):
    test(adata, group_by="Group", test_group="test", reference="control",
         k=k, total_RAM_available_gb=15)


def test_invalid_k_type(adata):
    with pytest.raises(ValueError):
        test(adata, group_by="Group", test_group="test", reference="control", k="foo")


def test_invalid_k_negative(adata):
    with pytest.raises(ValueError):
        test(adata, group_by="Group", test_group="test", reference="control", k=-1.5)


# ── output structure ──────────────────────────────────────────────────────────

def test_partner_column_added(adata):
    test(adata, group_by="Group", test_group="test", reference="control", k=10)
    assert "XMatch_partner_test_vs_control" in adata.obs.columns


def test_partner_column_groups(adata):
    test(adata, group_by="Group", test_group="test", reference="control", k=10)
    matched = adata[adata.obs["XMatch_partner_test_vs_control"].notna()]
    assert set(matched.obs["Group"].unique()) == {"test", "control"}


# ── numerical results ─────────────────────────────────────────────────────────

def test_results_simulated(linspace_adata):
    result = test(linspace_adata, group_by="Group", test_group="test", reference="control",
                  metric="euclidean", rank=False, k=10)
    assert result["p_value"] == pytest.approx(0.9999999999999821, abs=1e-12)
    assert result["z_score"] == pytest.approx(6.96419413859206, rel=1e-6)
    assert result["coverage"] == pytest.approx(1.0)


@pytest.mark.slow
def test_results_krumsiek11_paired():
    adata = sc.datasets.krumsiek11()
    result = test(adata, group_by="cell_type", test_group="Mo", reference="Ery",
                  k=100, metric="sqeuclidean", rank=False)
    assert result["p_value"] == pytest.approx(1.1679837230153187e-24, rel=1e-6)
    assert result["z_score"] == pytest.approx(-8.972174757807267, rel=1e-6)
    assert result["coverage"] == pytest.approx(1.0)


@pytest.mark.slow
def test_results_krumsiek11_vs_rest():
    adata = sc.datasets.krumsiek11()
    result = test(adata, group_by="cell_type", test_group="Mo", reference=None,
                  k=100, metric="sqeuclidean", rank=False)
    assert result["p_value"] == pytest.approx(6.1476610170941865e-53, rel=1e-6)
    assert result["z_score"] == pytest.approx(-17.975222537080874, rel=1e-6)
    assert result["coverage"] == pytest.approx(1.0)


@pytest.mark.slow
def test_results_krumsiek11_multigroup():
    adata = sc.datasets.krumsiek11()
    result = test(adata, group_by="cell_type", test_group=["Mo", "Ery"],
                  reference=["Mk", "Neu"], k=100, metric="sqeuclidean", rank=False)
    assert result["p_value"] == pytest.approx(9.668885179334899e-49, rel=1e-6)
    assert result["z_score"] == pytest.approx(-12.6688586787564, rel=1e-6)
    assert result["coverage"] == pytest.approx(1.0)


# ── estimate_peak_RAM_GB ──────────────────────────────────────────────────────

def test_ram_estimate_value():
    assert estimate_peak_RAM_GB(N=1000, k=100) == pytest.approx(1.974865094992254)


@pytest.mark.parametrize("kwargs,exc", [
    ({"N": 1000, "k": 1000}, ValueError),  # k >= N
    ({"N": -1,   "k": 10},   ValueError),  # negative N
    ({"N": 0,    "k": 10},   ValueError),  # zero N
    ({"N": 1000, "k": -5},   ValueError),  # negative k
    ({"N": 1000, "k": 0},    ValueError),  # zero k
])
def test_ram_estimate_invalid_inputs(kwargs, exc):
    with pytest.raises(exc):
        estimate_peak_RAM_GB(**kwargs)
