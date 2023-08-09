import joblib
import scvelo as scv
from params import Params
import scanpy as sc
import anndata as ad
import scanpy.external as sce
import utils.common as cm
import numpy as np
import matplotlib

matplotlib.use('macosx')
scv.set_figure_params()

if __name__ == '__main__':
    p = Params()
    adata = ad.read_h5ad(p.file_path)
    adata.layers = joblib.load(p.folder + "layers.p")

    # scv.pl.proportions(adata)
    sub_info = "manual_anno_L0"
    adata.obsm[f'X_{sub_info}_sub_umap'] = np.empty((len(adata), 2), dtype=np.float64)
    info_vals = adata.obs[sub_info].unique()
    for info_val in info_vals:
        print(f"{sub_info} {info_val}")
        adata_d = adata[adata.obs[sub_info] == info_val].copy()

        sc.tl.pca(adata_d)
        sce.pp.harmony_integrate(adata_d, 'pool', max_iter_harmony=10, max_iter_kmeans=20, block_size=0.01, sigma=0.2,
                                 epsilon_cluster=float('-inf'), epsilon_harmony=float('-inf'))
        sc.pp.neighbors(adata_d, n_neighbors=10, n_pcs=50, use_rep="X_pca_harmony")
        sc.tl.umap(adata_d)
        sc.pl.umap(adata_d)
        adata.obsm[f'X_{sub_info}_sub_umap'][adata.obs[sub_info] == info_val] = adata_d.obsm['X_umap']

    adata.layers.clear()
    cm.safe_save(adata, p.file_path)

