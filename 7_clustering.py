import scanpy as sc
import anndata as ad
import numpy as np
from params import Params

p = Params()
adata = ad.read_h5ad(p.file_path)
# adata.obs = adata.obs.drop([col for col in adata.obs.columns if 'louvain' in col], axis=1)
# adata.write(f'{folder}als_n1.h5ad')
for prefix in ["new_pca_"]: # "new_pca_",
    rep = "X_pca"
    if prefix == "harmony_":
        rep = "X_pca_harmony"
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50, use_rep=rep)
    print("done neighbors")
    unique_pools = np.unique(adata.obs['pool'])

    for resolution in [0.6]: # , 1.0, 2.0
        print(f"Resolution {resolution}")
        sc.tl.louvain(adata, key_added=prefix + "louvain_" + str(resolution), resolution=resolution)
        adata.obs[prefix + 'sub_louvain_' + str(resolution)] = ''
        for pool in unique_pools:
            print(pool)
            adata_p = adata[adata.obs["pool"] == pool].copy()
            if prefix == "new_pca_":
                sc.tl.pca(adata_p)
            sc.pp.neighbors(adata_p, n_neighbors=10, n_pcs=50, use_rep=rep)
            sc.tl.louvain(adata_p, key_added=prefix + "sub_louvain_" + str(resolution), resolution=resolution)
            adata.obs.loc[adata.obs['pool'] == pool,
                          prefix + 'sub_louvain_' + str(resolution)] = adata_p.obs[prefix
                                                                                   + 'sub_louvain_' + str(resolution)]

adata.write(p.file_path)