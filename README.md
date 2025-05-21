![scXMatch](./logo.svg)

**scXMatch** is a Python package that implements Rosenbaum's cross-match test using distance-based matching to assess statistical dependence between two groups of high-dimensional data. This is particularly useful in analyzing multivariate distributions in structured data, such as single-cell RNA-seq.

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
    return_matching=False
)
```

#### Description

Performs Rosenbaum’s matching-based test to determine if there is a statistically significant difference between two groups of samples using a distance-based graph matching approach.

#### Parameters

- **adata** (`anndata.AnnData` | `pandas.DataFrame`):  
  The dataset containing features in `adata.X` and group labels in `adata.obs[group_by]`.

- **group_by** (`str`):  
  Column name in `adata.obs` indicating group labels.

- **test_group** (`str` or `List[str]`):  
  The group being tested for statistical dependence.

- **reference** (`str` or `List[str]`, optional):  
  Reference group for comparison. If `None`, all non-`test_group` samples are used.

- **metric** (`str`, default: `"sqeuclidean"`):  
  Distance metric used for computing pairwise distances. Follows `scipy.spatial.distance.cdist` standards.

- **rank** (`bool`, default: `False`):  
  If `True`, rank-transform features before computing distances.

- **k** (`int`, optional, default: `100`):  
  Number of nearest neighbors to consider during graph construction. If `None`, all samples are considered.

- **return_matching** (`bool`, default: `False`):  
  If `True`, returns the matching graph and matched pairs along with test results.

#### Returns

- **p_value** (`float`):  
  p-value of the statistical test.

- **z_score** (`float`):  
  Standardized z-score of the test statistic.

- **relative_support** (`float`):  
  Proportion of samples involved in the matching.

- **G** (`graph-tool.Graph`, optional):  
  Constructed graph (only returned if `return_matching=True`).

- **matching** (`List[Tuple[int, int]]`, optional):  
  List of index pairs representing the matching (only returned if `return_matching=True`).

#### Raises

- `TypeError` if input is not an `AnnData` object.  
- `ValueError` if `test_group` is not found in the group labels.

---

### `scxmatch.approximate_k`

```python
scxmatch.approximate_k(total_RAM_available_gb, num_samples)
```

#### Description

Estimates the largest feasible number of neighbors `k` for kNN graph construction based on available RAM and number of samples.

#### Parameters

- **total_RAM_available_gb** (`float`):  
  Available RAM in gigabytes.

- **num_samples** (`int`):  
  Number of samples in your dataset.

#### Returns

- **k** (`int`):  
  Approximate value of `k` that can be used without exceeding memory limits.

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
