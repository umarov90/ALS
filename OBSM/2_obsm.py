import anndata as ad
import numpy as np
import scanpy as sc
import scanpy.external as sce
from params import Params
import utils.common as cm
import matplotlib
matplotlib.use('macosx')

p = Params()
rn = 0
adata = ad.read_h5ad(p.file_path)

sc.tl.pca(adata, use_highly_variable=True, random_state=rn)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50, random_state=rn)

sc.tl.umap(adata, random_state=rn)
sc.pl.umap(adata, color="manual_anno_L0")
adata.obsm['X_umap_original'] = adata.obsm['X_umap']

sce.pp.harmony_integrate(adata, 'pool', max_iter_harmony=10, max_iter_kmeans=20,
                         random_state=rn, block_size=0.01, sigma=0.2,
                         epsilon_cluster=float('-inf'), epsilon_harmony=float('-inf'))

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50, use_rep="X_pca_harmony", random_state=rn)
sc.tl.umap(adata, random_state=rn)
sc.pl.umap(adata, color="manual_anno_L0")
adata.obsm['X_umap_harmony'] = adata.obsm['X_umap']

sc.tl.tsne(adata, random_state=rn)
sc.pl.tsne(adata, color="MYC", use_raw=False)
adata.obsm['X_tsne_original'] = adata.obsm['X_tsne']

sc.tl.tsne(adata, use_rep="X_pca_harmony", random_state=rn)
adata.obsm['X_tsne_harmony'] = adata.obsm['X_tsne']
sc.pl.tsne(adata, color="manual_anno_L0")

del adata.obsm["X_umap"]
del adata.obsm["X_tsne"]
cm.safe_save(adata, p.file_path)