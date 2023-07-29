import multiprocessing
multiprocessing.set_start_method("fork")
import sys
import scvelo as scv
from params import Params
import anndata as ad
import joblib
import scanpy as sc
from utils import common as cm
import numpy as np
from scipy.sparse import csr_matrix
import matplotlib
matplotlib.use('macosx')


if __name__ == '__main__':
    p = Params()
    color_column = "manual_anno_L0"
    custom_palette = joblib.load(p.folder + color_column + "_palette.p")
    adata = ad.read_h5ad(p.file_path)
    # adata = cm.prepare_adata_layers(adata)
    # print("Layers copied")
    # joblib.dump(adata.layers, p.folder + "layers_m.p")
    adata.layers = joblib.load(p.folder + "layers_m.p")
    # scv.pl.proportions(adata)
    # adata.obsm['X_umap'] = adata.obsm['X_umap_harmony']
    # adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']

    sc.tl.pca(adata, use_highly_variable=True, n_comps=30)
    scv.pp.moments(adata)
    scv.tl.velocity(adata)
    scv.tl.velocity_graph(adata)
    scv.pl.velocity_embedding_stream(adata, basis='X_umap_harmony', color=color_column, palette=custom_palette,
                                     save="meta_ALL.png", dpi=300)
    scv.pl.velocity_embedding_stream(adata, basis='X_umap_harmony', color="day",
                                     save="meta_ALL_day.png", dpi=300)
