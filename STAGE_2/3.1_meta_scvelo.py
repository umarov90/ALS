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

    print(len(adata))
    adata = adata[adata.obs["pool"] == "pool34"].copy()
    # adata = adata[adata.obs["gender"] == "Female"].copy()
    print(len(adata))
    sc.tl.pca(adata, use_highly_variable=True, n_comps=30)
    scv.pp.moments(adata)
    scv.tl.recover_dynamics(adata, n_jobs=1, n_top_genes=2000) # , n_top_genes=2000
    scv.tl.velocity(adata, mode='dynamical')
    scv.tl.velocity_graph(adata, n_jobs=1)
    scv.tl.recover_latent_time(adata)
    try:
        scv.pl.velocity_embedding_stream(adata, basis='X_sc_umap', color=color_column, palette=custom_palette,
                                         save="meta_ALL_dynamical.png", dpi=300)
        scv.pl.velocity_embedding_stream(adata, basis='X_sc_umap', color="day",
                                         save="meta_ALL_day_dynamical.png", dpi=300)
        scv.pl.velocity_embedding_stream(adata, basis='X_sc_umap', color="latent_time", color_map="plasma",
                                         save="meta_ALL_latent_time.png", dpi=300)
    except:
        pass

    adata.layers.clear()
    # cm.safe_save(adata, p.file_path)
