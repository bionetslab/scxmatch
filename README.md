<img src="./logo.svg" alt="scXMatch" width="400"/>




![Python package Conda](https://github.com/bionetslab/scxmatch/actions/workflows/python-package-conda.yml/badge.svg)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/scxmatch/README.html)




**scXMatch** (single-cell cross match) is a Python package that implements Rosenbaum's cross-match test using distance-based matching to assess statistical dependence between two groups of high-dimensional data. This is particularly useful in analyzing multivariate distributions in structured data, such as single-cell RNA-seq.

This package provides a Python implementation inspired by the methodology described in [Rosenbaum (2005)](https://faculty.wharton.upenn.edu/wp-content/uploads/2012/04/Multivariate-distributions.pdf).

---

## Installation
Due to its dependence on graph-tool, this package can only be installed from conda, not from PyPI. The channels need to be specified.
```bash
conda install scxmatch -c conda-forge -c bioconda
```
---

## Requirements

- Python ≥ 3.9
- `anndata`
- `scanpy`
- `scipy`
- `graph-tool` $\geq$ 2.92

---

## API Documentation

### `scxmatch.test`

```python
scxmatch.test(
    adata,
    group_by,
    test_group,
    reference=None,
    metric="sqeuclidean",
    rank=False,
    k=100,
    total_RAM_available_gb=None
)
```

#### Description

Performs Rosenbaum’s matching-based test to determine if there is a statistically significant difference between two groups of samples using a distance-based graph matching approach.

#### Parameters
- `adata` (`anndata.AnnData`): The input data matrix. Features should be in `adata.X`, and group labels in `adata.obs[group_by]`.
- `group_by` (`str`): Column in `adata.obs` indicating group labels.
- `test_group` (`str` or `list of str`): The group(s) to be tested.
- `reference` (`str` or `list of str`, optional): The reference group(s). If `None`, all non-test samples are used as reference.
- `metric` (`str`, default `"sqeuclidean"`): Distance metric for matching. Follows `scipy.spatial.distance.cdist` standards.
- `rank` (`bool`, default `True`): If `True`, features are rank-transformed before distance computation.
- `k` (`int`, `"auto"`, or `"full"`, default `100`): Number of nearest neighbors to use for graph construction. If `full`, a full distance matrix will be calculated.
- `total_RAM_available_gb` (`float`, optional): Required if `k="auto"`.
  
#### Returns
- `p_value` (`float`): P-value from the Rosenbaum crossmatch test.
- `z_score` (`float`): Standardized test statistic.
- `relative_support` (`float`): Proportion of samples included in the matching.

#### Raises
- `TypeError`: If the input `adata` is not an `AnnData` object.
- `ValueError`: If `test_group` or `reference` contains values not present in `adata.obs[group_by]`.
- `ValueError`: If `k="auto"` and `total_RAM_available_gb` is not provided.
- `ValueError`: If `k` is not an integer, `"auto"`, or `"full"`.

#### Modifies:

- Modifies `adata.obs` **in-place** by adding the following columns:
  - `XMatch_partner_<test_group>_vs_<reference>`: The index of each sample’s matched partner in the MWMCM.
---

## Example Usage

```python
import anndata as ad
import scxmatch

# Load your AnnData object
adata = ad.read_h5ad("your_data.h5ad")

# Run test
p_val, z, support = scxmatch.test(
    adata=adata,
    group_by="cell_type",
    test_group="treated",
    reference="control",
    metric="euclidean",
    rank=True,
    k=50
)

print(f"P-value: {p_val:.4f}, Z-score: {z:.2f}, Support: {support:.2%}")
```

---

## Citation

If you use `scXMatch` in your research, please cite the original paper and our publication:

> Rosenbaum, P. R. (2005). An exact distribution-free test comparing two multivariate distributions based on adjacency. *Journal of the Royal Statistical Society: Series B*, 67(4), 515–530.


> TBD

---

## License

MIT License
