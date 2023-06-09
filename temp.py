import numpy as np
import anndata as ad
import scanpy as sc
from params import Params
import matplotlib
matplotlib.use('Agg')
import plotnine as p
from plotnine import *


p = Params()
adata = ad.read_h5ad(p.file_path)

rep = "X_pca"
# rep = "X_pca_harmony"

adata_p = adata[adata.obs["day"] == "D14"].copy()
if prefix == "new_pca_":
    sc.tl.pca(adata_p)
    print("pca")
sc.pp.neighbors(adata_p, n_neighbors=10, n_pcs=50, use_rep=rep)
print("neighbors")
sc.tl.umap(adata_p)
print("umap")
adata.obsm[f'X_{prefix}pool_sub_umap'][adata.obs['pool'] == pool] = adata_p.obsm['X_umap']

p.options.figure_size = (6, 5)
fig = (
    p.ggplot(p.aes(x = 'umap_1', y = 'umap_2', color = 'individual'), adata_p.obs)
    + p.geom_point(shape='.', size=0.1)
    + p.theme_minimal()
)

ggsave(fig, filename=f'figs/samples/sample_{sample}.png', dpi=300, facecolor="white")