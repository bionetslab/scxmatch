import anndata as ad
import numpy as np
from scipy.stats import rankdata
import warnings

from ._match import _kNN, _calculate_distances, _extract_matching, _construct_graph_via_kNN, _construct_graph_from_distances, _match, _add_partners_to_adata
from ._count import _cross_match_count, _get_p_value, _get_z_score, _get_relative_support, _rosenbaum_test

def test(adata, group_by, test_group, reference=None, metric="sqeuclidean", rank=False, k=100):
    """
    Perform Rosenbaum's matching-based test for checking the association between two groups 
    using a distance-based matching approach.

    Parameters:
    -----------
    adata : anndata.AnnData or pd.DataFrame
        The input `AnnData` object containing the samples and their respective features.
        The samples and their corresponding features should be stored in `adata.X` and the
        group labels in `adata.obs[group_by]`. 
        
    group_by : str
        The column in `data.obs` containing the group labels.
        The values of this column should include the `test_group` and potentially the `reference` group.

    test_group : str
        The group of interest that is being tested. This group will be compared against the `reference` group.

    reference : str, optional, default=None
        The group used as a comparison to the `test_group`. If set to None, all groups other than `test_group`
        are treated as the reference group.

    metric : str, optional, default="sqeuclidean"
        The distance metric used for calculating distances between the samples during the matching process. 
        It can be any valid metric recognized by `scipy.spatial.distance.cdist`.

    rank : bool, optional, default=True
        If `True`, ranks the features in the data matrix before performing the matching. This can help reduce
        the impact of varying scales of the features on the distance computation.
    
    k : int, optional, default=100
        The number of nearest neighbors to consider for each sample. This parameter is used to limit the
        number of samples considered during the matching process, which can help reduce computational
        complexity and memory usage. If set to None, all samples are considered.

    Returns:
    --------
    p_value : float
        The p-value, indicating the statistical significance of the observed matching.
        
    z_score : float
        The z-score, indicating the standardized difference between the observed and expected.

    relative_support : float   
        The relative support of the matching, calculated as the ratio of the number of
        samples covered by the matching and the total number of samples.
        
    Raises:
    -------
    TypeError : If the input `data` is not an `AnnData` object.
    ValueError : If the input `test_group` is not in the data.

    Notes:
    ------
    scXMatch describes how likely it is to observe cross matches of the `test_group` and the `reference`
    group in a Maximum-Cardinality-Maximum-Weight-Matching computed based on distances. 
    """
    if not isinstance(adata, ad.AnnData):
        raise TypeError("the input must be an AnnData object.")
    
    if adata.is_view:
        adata = adata.copy()
        adata.obs.reset_index(inplace=True)
        adata.obs.index = adata.obs.index.astype(str)
        warnings.warn("The input AnnData object is a view, matching results will be written to a copy of it.")
        
    if not isinstance(test_group, list): 
        test_group = [test_group]
    
    if reference != None:
        if not isinstance(reference, list): 
            reference = [reference]

    for t in test_group:
        if t not in adata.obs[group_by].values:
            raise ValueError(f"the test group {t} is not contained in your data.")
        

    if reference != None:       
        for r in reference:
            if r not in adata.obs[group_by].values:
                raise ValueError(f"the reference group {r} is not contained in your data.")
        subset = adata[adata.obs[group_by].isin(test_group + reference), :].copy()
    else:
        subset = adata.copy()
    
    if rank:
        print("computing variable-wise ranks.")
        subset.X = np.apply_along_axis(rankdata, axis=0, arr=subset.X)
    
    if not isinstance(k, int):
        if k is not None:
            raise ValueError("k must be an integer or None.")
 
    subset.obs["XMatch_group"] = np.where(subset.obs[group_by].isin(test_group), "test", "reference")
    
    if isinstance(k, int):
        _kNN(subset, k, metric)
    
    num_samples = len(subset)
    if isinstance(k, int):
        G = _construct_graph_via_kNN(subset)

    else:
        distances = _calculate_distances(subset.X, metric)
        G = _construct_graph_from_distances(distances)
    matching = _match(G, num_samples)
    
    _add_partners_to_adata(adata, subset, matching, reference, test_group)
        
    group_by = "XMatch_group"
    test_group = "test"
    
    p_value, a1, E, var, z_score, relative_support = _rosenbaum_test(Z=subset.obs[group_by], matching=matching, test_group=test_group)
    effect_strength_ratio = 1 - np.min([a1, E]) / E if E != 0 else np.nan
    
    return {
        "p_value": p_value,
        "a1": a1,
        "E": E,
        "var": var,
        "z_score": z_score,
        "coverage": relative_support,
        "effect_strength_ratio": effect_strength_ratio
    }

def estimate_peak_RAM_GB(N, k):  
    if N <= 0:
        raise ValueError("N must be a positive integer.")
    if k <= 0:
        raise ValueError("k must be a positive integer.")
    if k >= N:
        raise ValueError("k must be less than N.")

    beta_0 = 1.9611825714677176
    beta_1 = 1.368252352453625e-07
    return beta_0 + beta_1 * N * k