import scanpy as sc

def scanpy_setup(adata):
    if 'counts' in adata.layers:
        adata.X = adata.layers['counts'].copy()
    else:
        adata.layers['counts'] = adata.X.copy()
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var.highly_variable]
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    return adata
