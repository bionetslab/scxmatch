import numpy as np
import anndata as ad
import sys
from scipy.spatial.distance import cdist
from scxmatch import *
import scanpy as sc


np.random.seed(42)

def simulate_data(n_obs, n_var):
    samples = [np.random.normal(0, 1, n_var) for _ in range(n_obs)]
    adata = ad.AnnData(np.array(samples))
    return adata


def test_group_by(adata, test_group, reference, metric, rank, k):
    try:
        test(adata, group_by="foo", test_group=test_group, reference=reference,
             metric=metric, rank=rank, k=k)
    except KeyError:
        print(" +++++++++++++++ Successfully threw KeyError for invalid group_by. +++++++++++++++ ")
        return
    raise AssertionError("KeyError not raised for invalid group_by.")


def test_test_group(adata, group_by, test_group, test_group_2, reference, metric, rank, k):
    test(adata, group_by=group_by, test_group=test_group, reference=reference,
        metric=metric, rank=rank, k=k)
    print(" +++++++++++++++ Successfully worked with single group as test_group. +++++++++++++++ ")
    
    test(adata, group_by=group_by, test_group=[test_group, test_group_2], reference=reference,
        metric=metric, rank=rank, k=k)
    print(" +++++++++++++++ Successfully worked with list as test_group. +++++++++++++++ ")

    try:
        test(adata, group_by=group_by, test_group="foo", reference=reference,
             metric=metric, rank=rank, k=k)
    except ValueError:
        print(" +++++++++++++++ Successfully threw ValueError for invalid test_group. +++++++++++++++ ")
        return
    raise AssertionError("KeyError not raised for invalid test_group.")


def test_reference(adata, group_by, test_group, test_group_2, reference, metric, rank, k):
    test(adata, group_by=group_by, test_group=test_group, reference=reference,
        metric=metric, rank=rank, k=k)
    print(" +++++++++++++++ Successfully worked with single group as reference. +++++++++++++++ ")
    
    test(adata, group_by=group_by, test_group=test_group, reference=[reference, test_group_2],
        metric=metric, rank=rank, k=k)
    print(" +++++++++++++++ Successfully worked with list as test_group. +++++++++++++++ ")

    try:
        test(adata, group_by=group_by, test_group=test_group, reference="foo",
             metric=metric, rank=rank, k=k)
    except ValueError:
        print(" +++++++++++++++ Successfully threw ValueError for invalid reference. +++++++++++++++ ")
        return
    raise AssertionError("KeyError not raised for invalid test_group.")


def test_metrics(adata, group_by, test_group, reference, k, rank):
    metrics = ["euclidean", "sqeuclidean", "jaccard"]
    for metric in metrics:
        print(f"Testing metric: {metric}")
        test(adata, group_by=group_by, test_group=test_group, reference=reference,
             metric=metric, rank=rank, k=k)
    try:
        test(adata, group_by=group_by, test_group=test_group, reference=reference,
             metric="foo", rank=rank, k=k)
    except ValueError:
        print(" +++++++++++++++ Successfully threw ValueError for invalid metric. +++++++++++++++ ")
        return
    raise AssertionError("ValueError not raised for invalid metric.")


def test_rank(adata, group_by, test_group, reference, k, metric):
    for rank in [False, True]:
        test(adata, group_by=group_by, test_group=test_group, reference=reference,
                 metric=metric, rank=rank, k=k)
    print(" +++++++++++++++ Successfully worked for ranked and original values. +++++++++++++++ ")


def test_auto_k(adata, group_by, test_group, reference, rank):        
    test(adata, group_by=group_by, test_group=test_group, reference=reference, rank=rank, k="auto", total_RAM_available_gb=16)
    
    try:
        test(adata, group_by=group_by, test_group=test_group, reference=reference, rank=rank, k="auto", total_RAM_available_gb=0.1)
    except ValueError:
        print(" +++++++++++++++ Successfully threw ValueError for RAM under 0.21. +++++++++++++++ ")
        try:
            test(adata, group_by=group_by, test_group=test_group, reference=reference, rank=rank, k="auto", total_RAM_available_gb=-10)
        except ValueError:
            print(" +++++++++++++++ Successfully threw ValueError for negative RAM. +++++++++++++++ ")
            return
    raise AssertionError("ValueError not raised for RAM under 0.21 or negative RAM.")
    

def test_k(adata, group_by, test_group, reference, metric, rank):
    for k in [10, "auto", "full"]:
        test(adata, group_by=group_by, test_group=test_group, reference=reference,
            metric=metric, rank=rank, k=k, total_RAM_available_gb=15)
    try:
        test(adata, group_by=group_by, test_group=test_group, reference=reference,
            metric=metric, rank=rank, k="foo")    
    except ValueError:
        print(" +++++++++++++++ Successfully threw ValueError for invalid k. +++++++++++++++ ")
        try:
            test(adata, group_by=group_by, test_group=test_group, reference=reference,
                 metric=metric, rank=rank, k=-1.5)
        except ValueError:
            print(" +++++++++++++++ Successfully threw ValueError for negative k. +++++++++++++++ ")
            return
        return      
    raise AssertionError("ValueError not raised for invalid k.")
        
        
def test_added_columns(adata, group_by, test_group, reference, metric, rank, k):
    test(adata, group_by=group_by, test_group=test_group, reference=reference,
            metric=metric, rank=rank, k=k)
    if "XMatch_partner_test_vs_control" not in adata.obs.columns:
        raise AssertionError("XMatch_partner_test_vs_control column not added to adata.obs.")
    if set(adata[adata.obs["XMatch_partner_test_vs_control"].notna()].obs[group_by].unique()) != set([reference, test_group]):
        raise AssertionError("XMatch_partner_test_vs_control column does not contain expected values.")   
    

def test_results_simulated():
    data = ad.AnnData(np.linspace(0, 100, 1000).reshape(100, 10))
    data.obs["Group"] = [0, 1] * 50
    p, z, s = test(data, group_by="Group", test_group="test", reference="control",
                   metric="euclidean", rank=False, k=10)

    assert np.isclose(p, 0.9999999999999821, atol=1e-16), f"p value mismatch on simulated data"
    assert np.isclose(z, 6.96419413859206, atol=1e-16), f"z score mismatch on simulated data"
    assert np.isclose(s, 1.0, atol=1e-16), f"support mismatch on simulated data"
    print(" +++++++++++++++ Successfully tested simulated data. +++++++++++++++ ")


def test_results_krumsiek11():
    adata = sc.datasets.krumsiek11()
    p, z, s = test(adata, group_by="cell_type", test_group="Mo", reference="Ery", k=100, metric="sqeuclidean", rank=False)
    assert np.isclose(p, 1.1679837230153187e-24, atol=1e-16), f"p value mismatch on krumsiek11 data"
    assert np.isclose(z, -8.972174757807267, atol=1e-16), f"z score mismatch on krumsiek11 data"
    assert np.isclose(s, 1.0, atol=1e-16), f"support mismatch on krumsiek11 data"

    p, z, s = test(adata, group_by="cell_type", test_group="Mo", reference=None, k=100, metric="sqeuclidean", rank=False)
    assert np.isclose(p, 6.1476610170941865e-53, atol=1e-16), f"p value mismatch on krumsiek11 data"
    assert np.isclose(z, -17.975222537080874, atol=1e-16), f"z score mismatch on krumsiek11 data"
    assert np.isclose(s, 1.0, atol=1e-16), f"support mismatch on krumsiek11 data"    
    
    p, z, s = test(adata, group_by="cell_type", test_group=["Mo", "Ery"], reference=["Mk", "Neu"], k=100, metric="sqeuclidean", rank=False)
    assert np.isclose(p, 9.668885179334899e-49, atol=1e-16), f"p value mismatch on krumsiek11 data"
    assert np.isclose(z, -12.6688586787564, atol=1e-16), f"z score mismatch on krumsiek11 data"
    assert np.isclose(s, 1.0, atol=1e-16), f"support mismatch on krumsiek11 data"    
    print(" +++++++++++++++ Successfully tested krumsiek11 data. +++++++++++++++ ")
    print(adata)

    
def main():
    test_results_simulated()
    test_results_krumsiek11()
     
    n_obs = 100
    n_var = 5
    k = 10
    group_by = "Group"
    reference = "control"
    test_group = "test"
    test_group_2 = "test_2"

    metric = "euclidean"
    rank = False

    adata = simulate_data(n_obs, n_var)
    adata.obs[group_by] = np.random.choice([reference, test_group, test_group_2], size=n_obs)

    # Individual tests
    test_group_by(adata, test_group, reference, metric, rank, k)
    test_test_group(adata, group_by, test_group, test_group_2, reference, metric, rank, k)
    test_reference(adata, group_by, test_group, test_group_2, reference, metric, rank, k)
    test_metrics(adata, group_by, test_group, reference, k, rank)
    test_rank(adata, group_by, test_group, reference, k, metric)
    test_auto_k(adata, group_by, test_group, reference, rank)    
    test_k(adata, group_by, test_group, reference, metric, rank)
    test_added_columns(adata, group_by, test_group, reference, metric, rank, k)
    print("All tests passed successfully!")
    
    
if __name__ == "__main__":
    main()
