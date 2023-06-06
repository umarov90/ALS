import numpy as np
import anndata as ad
import scanpy as sc
from params import Params
import matplotlib
matplotlib.use('Agg')

p = Params()
adata = ad.read_h5ad(p.file_path)
for prefix in ["harmony_", "old_pca_"]: # "new_pca_",
    print(prefix)
    rep = "X_pca"
    if prefix == "harmony_":
        rep = "X_pca_harmony"
    adata.obsm[f'X_{prefix}sub_umap'] = np.empty((adata.n_obs, 2), dtype=np.float64)

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

    unique_samples = np.unique(adata.obs['sample_id'])
    for sample in unique_samples:
        print(sample)
        adata_p = adata[adata.obs["sample_id"] == sample].copy()
        if prefix == "new_pca_":
            sc.tl.pca(adata_p)
            print("pca")
        sc.pp.neighbors(adata_p, n_neighbors=10, n_pcs=50, use_rep=rep)
        print("neighbors")
        sc.tl.umap(adata_p)
        print("umap")
        adata.obsm[f'X_{prefix}sample_sub_umap'][adata.obs['sample_id'] == sample] = adata_p.obsm['X_umap']

del adata.obsm["X_umap"]
adata.write(p.file_path)

# import plotnine as p
# from plotnine.save import ggsave
#     p.options.figure_size = (6, 5)
#     fig = (
#         p.ggplot(p.aes(x = 'umap_1', y = 'umap_2', color = 'individual'), adata_p.obs)
#         + p.geom_point(shape='.', size=0.1)
#         + p.theme_minimal()
#     )
#
#     ggsave(fig, filename=f'figs/samples/sample_{sample}.png', dpi=300, facecolor="white")
