import numpy as np
import anndata as ad
import sys
from scipy.spatial.distance import cdist

sys.path.append("..")
from scxmatch import *

np.random.seed(42)

def simulate_data(n_obs, n_var):
    samples = [np.random.normal(0, 1, n_var) for _ in range(n_obs)]
    adata = ad.AnnData(np.array(samples))
    return adata


def test_group_by(adata, test_group, reference, metric, rank, k, return_matching):
    try:
        test(adata, group_by="foo", test_group=test_group, reference=reference,
             metric=metric, rank=rank, k=k, return_matching=return_matching)
    except KeyError:
        print(" +++++++++++++++ Successfully threw KeyError for invalid group_by. +++++++++++++++ ")
        return
    raise AssertionError("KeyError not raised for invalid group_by.")


def test_test_group(adata, group_by, test_group, test_group_2, reference, metric, rank, k, return_matching):
    test(adata, group_by=group_by, test_group=test_group, reference=reference,
        metric=metric, rank=rank, k=k, return_matching=return_matching)
    print(" +++++++++++++++ Successfully worked with single group as test_group. +++++++++++++++ ")
    
    test(adata, group_by=group_by, test_group=[test_group, test_group_2], reference=reference,
        metric=metric, rank=rank, k=k, return_matching=return_matching)
    print(" +++++++++++++++ Successfully worked with list as test_group. +++++++++++++++ ")

    try:
        test(adata, group_by=group_by, test_group="foo", reference=reference,
             metric=metric, rank=rank, k=k, return_matching=return_matching)
    except ValueError:
        print(" +++++++++++++++ Successfully threw ValueError for invalid test_group. +++++++++++++++ ")
        return
    raise AssertionError("KeyError not raised for invalid test_group.")


def test_reference(adata, group_by, test_group, test_group_2, reference, metric, rank, k, return_matching):
    test(adata, group_by=group_by, test_group=test_group, reference=reference,
        metric=metric, rank=rank, k=k, return_matching=return_matching)
    print(" +++++++++++++++ Successfully worked with single group as reference. +++++++++++++++ ")
    
    test(adata, group_by=group_by, test_group=test_group, reference=[reference, test_group_2],
        metric=metric, rank=rank, k=k, return_matching=return_matching)
    print(" +++++++++++++++ Successfully worked with list as test_group. +++++++++++++++ ")

    #try:
    #    test(adata, group_by=group_by, test_group=test_group, reference="foo",
    #         metric=metric, rank=rank, k=k, return_matching=return_matching)
    #except ValueError:
    #    print(" +++++++++++++++ Successfully threw ValueError for invalid test_group. +++++++++++++++ ")
    #    return
    #raise AssertionError("KeyError not raised for invalid test_group.")


def test_metrics(adata, group_by, test_group, reference, k, rank, return_matching):
    metrics = ["euclidean", "sqeuclidean", "jaccard"]
    for metric in metrics:
        print(f"Testing metric: {metric}")
        test(adata, group_by=group_by, test_group=test_group, reference=reference,
             metric=metric, rank=rank, k=k, return_matching=return_matching)
    try:
        test(adata, group_by=group_by, test_group=test_group, reference=reference,
             metric="foo", rank=rank, k=k, return_matching=return_matching)
    except ValueError:
        print(" +++++++++++++++ Successfully threw ValueError for invalid metric. +++++++++++++++ ")
        return
    raise AssertionError("ValueError not raised for invalid metric.")


def test_rank(adata, group_by, test_group, reference, k, metric, return_matching):
    for rank in [False, True]:
        test(adata, group_by=group_by, test_group=test_group, reference=reference,
                 metric=metric, rank=rank, k=k, return_matching=return_matching)
    print(" +++++++++++++++ Successfully worked for ranked and original values. +++++++++++++++ ")


def test_return_matching(adata, group_by, test_group, reference, k, metric, rank):
    _, _, _ = test(adata, group_by=group_by, test_group=test_group, reference=reference,
                metric=metric, rank=rank, k=k, return_matching=False)
    (_, _, _), _, _ = test(adata, group_by=group_by, test_group=test_group, reference=reference,
            metric=metric, rank=rank, k=k, return_matching=True)
    print(" +++++++++++++++ Successfully worked with and without returning matching. +++++++++++++++ ")



def test_approximate_k(num_samples):
    assert 671914 == approximate_k(16, num_samples)
    assert 425 == approximate_k(0.22, num_samples)
    print(" +++++++++++++++ Successfully calculated k. +++++++++++++++ ")

    ## assert 99 == approximate_k(total_RAM_available_gb, num_samples)
    ##try:
    ##    approximate_k(0.00001, num_samples)
    ##except ValueError:
    ##    print(" +++++++++++++++ Successfully threw ValueError for RAM under 0.21. +++++++++++++++ ")
    ##   try:
    ##        approximate_k(total_RAM_available_gb, 0)
    ##    except ValueError:
    ##        print(" +++++++++++++++ Successfully threw ValueError for len(adata) == 0. +++++++++++++++ ")
    ##       return
    ##    raise AssertionError("ValueError not raised for len(adata) == 0.")
    ##raise AssertionError("ValueError not raised for RAM under 0.21.")
    

def main():
    n_obs = 100
    n_var = 5
    k = 10
    group_by = "Group"
    reference = "control"
    test_group = "test"
    test_group_2 = "test_2"

    metric = "euclidean"
    rank = False
    return_matching = False

    adata = simulate_data(n_obs, n_var)
    adata.obs[group_by] = np.random.choice([reference, test_group, test_group_2], size=n_obs)

    # Individual tests
    test_group_by(adata, test_group, reference, metric, rank, k, return_matching)
    test_test_group(adata, group_by, test_group, test_group_2, reference, metric, rank, k, return_matching)
    test_reference(adata, group_by, test_group, test_group_2, reference, metric, rank, k, return_matching)
    test_metrics(adata, group_by, test_group, reference, k, rank, return_matching)
    test_rank(adata, group_by, test_group, reference, k, metric, return_matching)
    test_return_matching(adata, group_by, test_group, reference, k, metric, rank)
    
    test_approximate_k(len(adata))


if __name__ == "__main__":
    main()
