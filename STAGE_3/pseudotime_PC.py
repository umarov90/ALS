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
sc.pl.embedding(basis="X_sc_umap", adata=adata, color="pseudotime")
sc.tl.pca(adata, use_highly_variable=True, random_state=rn)
sce.pp.harmony_integrate(adata, 'pool', max_iter_harmony=20, max_iter_kmeans=40,
                         random_state=rn, block_size=0.01, sigma=0.2,
                         epsilon_cluster=float('-inf'), epsilon_harmony=float('-inf'))
sc.pl.pca(adata, color="day", dimensions=[(0, 1), (1, 2)])
adata.obs["pseudotime"] = adata.obsm["X_pca_harmony"][:,0]
sc.pl.embedding(basis="X_sc_umap", adata=adata, color="pseudotime")
# cm.safe_save(adata, p.file_path)