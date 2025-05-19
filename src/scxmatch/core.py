import anndata as ad
import numpy as np
from scipy.stats import rankdata

from .match import *
from .count import *


def test(adata, group_by, test_group, reference=None, metric="sqeuclidean", rank=False, k=100, return_matching=False):
    """
    Perform Rosenbaum's matching-based test for checking the association between two groups 
    using a distance-based matching approach.

    Parameters:
    -----------
    adata : anndata.AnnData or pd.DataFrame
        The input `AnnData` object containing the samples and their respective features.
        The samples and their corresponding features should be stored in `data.X` and the
        group labels in `data.obs[group_by]`. 
        
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

    return_matching : bool, optional, default=False
        If `True`, returns the matching object along with the p-value, the z-score, and the relative support.
        If `False`, only the p-value, z-score, and relative support are returned.

    Returns:
    --------
    p_value : float
        The p-value, indicating the statistical significance of the observed matching.
        
    z_score : float
        The z-score, indicating the standardized difference between the observed and expected.

    relative_support : float   
        The relative support of the matching, calculated as the ratio of the number of
        samples covered by the matching and the total number of samples.
    
    G : graph-tool `Graph` object
        The graph representation of the samples and their distances, constructed from the input data.
        This is returned only if `return_matching` is set to `True`.
        
    matching : list of tuples
        The matching pairs of samples, where each pair is represented as a tuple of indices.
        This is returned only if `return_matching` is set to `True`.
        
    
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
        _kNN(adata, k, metric)
    
    num_samples = len(adata)
    if k:
        G = _construct_graph_via_kNN(adata)

    else:
        distances = _calculate_distances(adata.X, metric)
        G = _construct_graph_from_distances(distances)
    matching = _match(G, num_samples)

    if return_matching:
        return _rosenbaum_test(Z=adata.obs[group_by], matching=matching, test_group=test_group), G, matching
    return _rosenbaum_test(Z=adata.obs[group_by], matching=matching, test_group=test_group)



def approximate_k(total_RAM_available_gb, num_samples):
    """
    Calculate the maximum number of neighbors (k) that can be used for kNN matching based on the available RAM and the number of samples.
    The function uses the parameters from the paper.
    Parameters:
    ----------          
    total_RAM_available_gb : float
        The total available RAM in gigabytes.
    num_samples : int

        The number of samples in the dataset.
        
    Returns:
    -------
    k : int
        An approximation of the maximum number of neighbors (k) that can be used for kNN matching.
    """
    k = 10**7 * (total_RAM_available_gb - 0.21) / (2.35 * num_samples)
    return int(np.floor(k))
    
    