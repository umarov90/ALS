import scvelo as scv
from utils.params import Params
import scanpy as sc
import anndata as ad
import numpy as np
from scipy.sparse import csr_matrix
import matplotlib
matplotlib.use('macosx')
scv.set_figure_params()


if __name__ == '__main__':
    p = Params()
    adata = ad.read_h5ad(p.file_path)

    sc.tl.pca(adata, n_comps=30)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
    scv.tl.louvain(adata, key_added="clusters", resolution=0.6)
    scv.tl.umap(adata)

    # scv.pl.proportions(adata)

    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.tl.velocity(adata)
    scv.tl.velocity_graph(adata, approx=True)

    scv.pl.velocity_embedding_stream(adata, basis='umap')

    print("All done")
