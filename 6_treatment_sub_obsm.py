import numpy as np
import anndata as ad
import scanpy as sc
from params import Params
import matplotlib
matplotlib.use('Agg')

p = Params()
adata = ad.read_h5ad(p.file_path)
adata.obs["day_pool"] = adata.obs["day"].astype(str) + "_" + adata.obs["pool"].astype(str)
for prefix in ["new_pca_"]: # , "harmony_", "old_pca_"
    print(prefix)
    rep = "X_pca"
    if prefix == "harmony_":
        rep = "X_pca_harmony"

    adata.obsm[f'X_{prefix}day_pool_sub_umap'] = np.empty((len(adata), 2), dtype=np.float64)
    unique_day_pool = np.unique(adata.obs['day_pool'])
    for day_pool in unique_day_pool:
        print(day_pool)
        adata_p = adata[adata.obs["day_pool"] == day_pool].copy()
        if prefix == "new_pca_":
            sc.tl.pca(adata_p)
            print("pca")
        sc.pp.neighbors(adata_p, n_neighbors=10, n_pcs=50, use_rep=rep)
        print("neighbors")
        sc.tl.umap(adata_p)
        print("umap")
        adata.obsm[f'X_{prefix}day_pool_sub_umap'][adata.obs['day_pool'] == day_pool] = adata_p.obsm['X_umap']

adata.write(p.file_path)


