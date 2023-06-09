import numpy as np
import anndata as ad
import scanpy as sc
from params import Params
import matplotlib
matplotlib.use('Agg')

p = Params()
adata = ad.read_h5ad(p.file_path)
adata.obs["sample"] = adata.obs["day"].astype(str) + "_" + adata.obs["treatment"].astype(str)
for prefix in ["harmony_", "old_pca_"]: # "new_pca_",
    print(prefix)
    rep = "X_pca"
    if prefix == "harmony_":
        rep = "X_pca_harmony"

    adata.obsm[f'X_{prefix}pool_sub_umap'] = np.empty((len(adata), 2), dtype=np.float64)
    unique_pools = np.unique(adata.obs['pool'])
    for pool in unique_pools:
        print(pool)
        adata_p = adata[adata.obs["pool"] == pool].copy()
        if prefix == "new_pca_":
            sc.tl.pca(adata_p)
            print("pca")
        sc.pp.neighbors(adata_p, n_neighbors=10, n_pcs=50, use_rep=rep)
        print("neighbors")
        sc.tl.umap(adata_p)
        print("umap")
        adata.obsm[f'X_{prefix}pool_sub_umap'][adata.obs['pool'] == pool] = adata_p.obsm['X_umap']

    adata.obsm[f'X_{prefix}sample_sub_umap'] = np.empty((len(adata), 2), dtype=np.float64)
    unique_samples = np.unique(adata.obs['sample'])
    for sample in unique_samples:
        print(sample)
        adata_p = adata[adata.obs["sample"] == sample].copy()
        if prefix == "new_pca_":
            sc.tl.pca(adata_p)
            print("pca")
        sc.pp.neighbors(adata_p, n_neighbors=10, n_pcs=50, use_rep=rep)
        print("neighbors")
        sc.tl.umap(adata_p)
        print("umap")
        adata.obsm[f'X_{prefix}sample_sub_umap'][adata.obs['sample'] == sample] = adata_p.obsm['X_umap']

adata.write(p.file_path)


