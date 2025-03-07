import anndata as ad
from itertools import chain
from math import comb, factorial, log, exp
import numpy as np
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import rankdata

from scxmatch.matching import *


def cross_match_count(Z, matching, test_group):
    print("counting cross matches.")
    pairs = [(Z.iloc[i], Z.iloc[j]) for (i, j) in matching]
    filtered_pairs = [pair for pair in pairs if (pair[0] == test_group) ^ (pair[1] == test_group)] # cross-match pairs contain test group exactly once
    a1 = len(filtered_pairs)
    return a1


def get_p_value(a1, n, N, I):
    p_value = 0
    for A1 in range(a1 + 1):  # For all A1 <= a1
        A2 = (n - A1) / 2 
        A0 = I - (n + A1) / 2 

        if int(A0) != A0:
            continue
        if int(A2) != A2:
            continue 
        if A0 < 0 or A2 < 0:
            continue  

        A0 = int(A0)
        A2 = int(A2)
        
        log_numerator = A1 * log(2) + log(factorial(I))
        log_denominator = log(comb(N, n)) + log(factorial(A0)) + log(factorial(A1)) + log(factorial(A2))
        p_value += exp(log_numerator - log_denominator)

    return p_value

    
def get_z_score(a1, n, N):
    m = N - n
    E = n * m / (N - 1) # Eq. 3 in Rosenbaum paper
    var = 2 * n * (n - 1) * m * (m - 1) / ((N - 3) * (N - 1)**2)
    z = (a1 - E) / np.sqrt(var)
    return z


def get_relative_support(N, Z):
    max_support = len(Z) - (len(Z) % 2)
    return N / max_support
    

def rosenbaum_test(Z, matching, test_group):
    used_elements = list(chain.from_iterable(matching))
    n = sum(1 for el in used_elements if Z.iloc[el] == test_group)
    N = len(matching) * 2
    I = len(matching)

    a1 = cross_match_count(Z, matching, test_group)
    
    p_value = get_p_value(a1, n, N, I)
    z_score = get_z_score(a1, n, N)
    relative_support = get_relative_support(N, Z)
    return p_value, z_score, relative_support


def kNN(adata, k, metric):
    print("calculating kNN graph.")
    if sp.issparse(adata.X):
        adata.X = adata.X.toarray()  # Convert only if it's sparse
    sc.pp.neighbors(adata, n_neighbors=k, metric=metric, n_pcs=0, transformer='pynndescent')


def rosenbaum(adata, group_by, test_group, reference=None, metric="sqeuclidean", rank=False, k=None, return_matching=False):
    """
    Perform Rosenbaum's matching-based test for checking the association between two groups 
    using a distance-based matching approach.

    Parameters:
    -----------
    data : anndata.AnnData or pd.DataFrame
        The input data containing the samples and their respective features. If the input is an
        `AnnData` object, the samples and their corresponding features should be stored in `data.X` and the
        group labels in `data.obs[group_by]`. If using a `pandas.DataFrame`, the group labels should be in the
        column specified by `group_by`, and the feature matrix should be the remaining columns.

    group_by : str
        The column in `data.obs` or `data` (in case of a `pandas.DataFrame`) containing the group labels.
        The values of this column should include the `test_group` and potentially the `reference` group.

    test_group : str
        The group of interest that is being tested for association. This group will be compared against the `reference` group.

    reference : str, optional, default="rest"
        The group used as a comparison to the `test_group`. If set to "rest", all groups other than `test_group`
        are treated as the reference group.

    metric : str, optional, default="mahalanobis"
        The distance metric used for calculating distances between the samples during the matching process. 
        It can be any valid metric recognized by `scipy.spatial.distance.cdist`.

    rank : bool, optional, default=True
        If `True`, ranks the features in the data matrix before performing the matching. This can help reduce
        the impact of varying scales of the features on the distance computation.

    Returns:
    --------
    p_value : float
        The p-value from Rosenbaum's test, indicating the statistical significance of the observed matching.

    a1 : int
        The count of cross-matched pairs that contain `test_group` exactly once. This is used to compute the p-value.

    Raises:
    -------
    TypeError : If the input `data` is neither an `AnnData` object nor a `pandas.DataFrame`.
    ValueError : If the input `test_group` is not in the data.

    Notes:
    ------
    Rosenbaum's test compares how likely it is to observe a matching between the `test_group` and the `reference`
    group, using a matching algorithm based on distance metrics (such as "mahalanobis"). The test computes a p-value
    based on the number of cross-matched pairs between the two groups.

    The function internally uses the `match_samples` function to compute a matching of the samples based on the chosen
    distance metric. The resulting matching is then used in the `rosenbaum_test` to calculate the p-value.
    """

    if not isinstance(adata, ad.AnnData):
        raise TypeError("the input must be an AnnData object or a pandas DataFrame.")
        
    if not isinstance(test_group, list): 
        test_group = [test_group]
    
    if reference != None:
        if not isinstance(reference, list): 
            reference = [reference]

    for t in test_group:
        if t not in adata.obs[group_by].values:
            raise ValueError(f"the test group {t} is not contained in your data.")
        
    if reference != None:       
        adata = adata[adata.obs[group_by].isin(test_group + reference), :]
        
    if rank:
        print("computing variable-wise ranks.")
        adata.X = np.apply_along_axis(rankdata, axis=0, arr=adata.X)
    

    adata.obs["XMatch_group"] = np.where(adata.obs[group_by].isin(test_group), "test", "reference")
    
    group_by = "XMatch_group"
    test_group = "test"
    print(adata.obs[group_by].value_counts())
       
    if k:
        kNN(adata, k, metric)
    
    num_samples = len(adata)
    if k:
        G = construct_graph_via_kNN(adata)

    else:
        distances = calculate_distances(adata.X, metric)
        G = construct_graph_from_distances(distances)
    matching = match(G, num_samples)

    if return_matching:
        return rosenbaum_test(Z=adata.obs[group_by], matching=matching, test_group=test_group), G, matching
    return rosenbaum_test(Z=adata.obs[group_by], matching=matching, test_group=test_group)