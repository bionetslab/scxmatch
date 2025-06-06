import numpy as np
from graph_tool.topology import max_cardinality_matching
import graph_tool.all as gt
from scipy.spatial.distance import cdist as cpu_cdist
from scipy.sparse import csr_matrix, tril, triu, issparse
import scanpy as sc


try:
    from cupyx.scipy.spatial.distance import cdist as gpu_cdist
    import cupy as cp
    GPU = True
    print("found cupy installation, will try use the GPU to calculate the full distance matrix.")
except:
    from scipy.spatial.distance import cdist as cpu_cdist
    GPU = False
    print("will use the CPU to calculate the distance matrix.")
    pass



def _kNN(adata, k, metric):
    print("calculating kNN graph.")
    if issparse(adata.X):
        adata.X = adata.X.toarray()  # Convert only if it's sparse
    sc.pp.neighbors(adata, n_neighbors=k, metric=metric, n_pcs=0, transformer='pynndescent')




def _calculate_distances(samples, metric):
    if not isinstance(samples, np.ndarray):  # Check if it's a scipy sparse matrix
        samples = samples.toarray()

    try:
        if GPU:
            print("trying to use GPU to calculate distance matrix.")
        distances = cp.asnumpy(gpu_cdist(cp.array(samples), cp.array(samples), metric=metric)) 

    except:
        if GPU:
            print("using CPU to calculate distance matrix due to chosen metric.")
        else:
            print("using CPU to calculate distance matrix.")
        distances = cpu_cdist(samples, samples, metric=metric)

    num_samples = len(samples)

    if num_samples % 2 != 0: # with an uneven number of samples, a minimal-distance column is added 
        distances = np.pad(distances, [(0, 1), (0, 1)], mode='constant', constant_values=0)
        num_samples += 1

    max_distance = np.max(distances)
    distances = max_distance + 1 - distances
    return distances


def _extract_matching(matching_map):
    matching_list = []
    matched = set()  # To keep track of processed vertices

    for v in matching_map.get_array().nonzero()[0]:  # Only consider vertices with matches
        partner = matching_map[v]
        if partner != -1 and partner not in matched:
            try:
                v = v.item()
            except:
                pass
            try: 
                partner = partner.item()
            except:
                pass
            matching_list.append((v, partner))
            matched.add(v)
            matched.add(partner)
    return matching_list



def _construct_graph_from_distances(distances):
    num_samples = distances.shape[0]
    print("creating distance graph with", num_samples, "samples")
    transposed_distances = distances.transpose()
    combined_distances = np.maximum(distances, transposed_distances)
    sparse_weights = csr_matrix(combined_distances)
    G = gt.Graph(sparse_weights, directed=False)
    return G


def _construct_graph_via_kNN(adata):
    distances = adata.obsp["distances"]
    max_dist = distances.max() 
    # only transform non-zero entries
    distances.data = max_dist + 1 - distances.data # transform so that weight minimization ~ weight maximization
    # the following seems a little cumbersome, however, if you pass a 
    # csr matrix to the graph-tool Graph constructor with directed=False, 
    # it will automatically ignore the lower diagonal. The distances matrix 
    # is technically directed, because a can be b's closest neighbor while b is not a's.
    # For the max cardinality min weight matching, we are generous and make all
    # directed edges undirected edges. 

    lower_tri = tril(distances)
    upper_tri = triu(distances)

    transposed_lower_tri = lower_tri.transpose()
    sparse_distances = transposed_lower_tri.maximum(upper_tri)

    G = gt.Graph(sparse_distances, directed=False) 
    return G


def _match(G, num_samples):
    matching = max_cardinality_matching(G, weight=G.edge_properties["weight"], minimize=False) # "minimize=True" only works with a heuristic, therefore we use (max_distance + 1 - distance_ij) and maximize 
    matching_list = _extract_matching(matching)
    matching_list = [p for p in matching_list if ((p[0] < num_samples) and (p[1] < num_samples))]
    return matching_list
            

def _add_partners_to_adata(adata, subset, matching_list, reference, test_group):
    if reference is None:
        ref_string = "rest"
    else:
        ref_string = "_".join(reference)
    test_string = "_".join(test_group)
    
    adata.obs[f"XMatch_partner_{test_string}_vs_{ref_string}"] = np.nan
    adata.obs[f"XMatch_partner_{test_string}_vs_{ref_string}"] = adata.obs[f"XMatch_partner_{test_string}_vs_{ref_string}"].astype("object")

    for sample_a, sample_b in matching_list:
        sampla_a_index = subset.obs.iloc[sample_a].name
        sampla_b_index = subset.obs.iloc[sample_b].name
        adata.obs.loc[sampla_a_index, f"XMatch_partner_{test_string}_vs_{ref_string}"] = sampla_b_index
        adata.obs.loc[sampla_b_index, f"XMatch_partner_{test_string}_vs_{ref_string}"] = sampla_a_index