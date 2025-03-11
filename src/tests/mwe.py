import numpy as np
import anndata as ad
import sys
sys.path.append("../src")
from scxmatch import *

def simulate_data(n_obs, n_var):
    samples = [np.random.normal(0, 1, n_var) for _ in range(n_obs)]
    adata = ad.AnnData(np.array(samples))
    return adata


def main():
    n_obs = 100
    n_var = 2
    k = 10
    group_by = "Group"
    reference = "control"
    test_group="test"

    adata = simulate_data(n_obs, n_var)
    adata.obs[group_by] = np.random.choice([reference, test_group], size=n_obs)
    p, z, s = rosenbaum(adata, group_by=group_by, test_group=test_group, reference=reference, metric="sqeuclidean", rank=False, k=k, return_matching=False)
    p, z, s = rosenbaum(adata, group_by=group_by, test_group=test_group, reference=reference, metric="sqeuclidean", rank=False, k=None, return_matching=False)


if __name__ == "__main__":
    main()