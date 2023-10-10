import anndata as ad
import numpy as np
import scanpy as sc
import scanpy.external as sce
import joblib
import scvelo as scv
from params import Params
import utils.common as cm
import matplotlib
matplotlib.use('macosx')


if __name__ == '__main__':
    p = Params()
    rn = 0
    adata = ad.read_h5ad(p.file_path)
    # adata = ad.read_h5ad(p.folder + "dpt.h5ad")
    # print(len(adata))
    # adata = adata[adata.obs["exclude_dpt_meta"] != "yes"]
    print(len(adata))
    adata = adata[adata.obs["manual_anno_L1"].isin(['iPSC.D00', 'NSC.D04', 'NSC.D07', 'NSC.D014', 'pMN.D14', 'mMN.D14'])].copy()
    adata = adata[~adata.obs["manual_anno_L2"].isin(['NSC.D14.1', 'NSC.D14.5'])].copy()
    sc.pl.embedding(basis="X_sc_umap", adata=adata, color="manual_anno_L2")
    adata.uns['iroot'] = np.flatnonzero(adata.obs['manual_anno_L0'] == 'iPSC')[0]
    sc.tl.pca(adata, use_highly_variable=True, n_comps=30)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30, method='gauss')
    sc.tl.diffmap(adata, n_comps=10)
    sc.tl.dpt(adata, n_dcs=10) # n_branchings
    # cm.safe_save(adata, p.folder + "als_sc_dpt.h5ad")
    sc.pl.embedding(basis="X_sc_umap", adata=adata, color="dpt_pseudotime")
    exit(0)
    #
    # adata.layers = joblib.load(p.folder + "layers_m.p")
    # scv.pp.moments(adata)
    # scv.tl.velocity(adata)
    # scv.tl.velocity_graph(adata)
    #
    # # cm.safe_save(adata, p.folder + "paga.h5ad")
    # adata.uns['velocity_graph']
    adata = ad.read_h5ad(p.folder + "paga.h5ad")
    sc.tl.paga(adata, groups='manual_anno_L0', use_rna_velocity=True, copy=False)
    sc.tl.paga(adata, groups='manual_anno_L0', use_rna_velocity=False, copy=False)
    # a = sc.tl.paga(adata, groups='manual_anno_L1', use_rna_velocity=False, copy=True)
    # adata.uns['connectivities'] = a
    # adata.uns['paga']['transitions_confidence']
    # adata.uns['paga']['connectivities']
    sc.pl.paga(adata, transitions="transitions_confidence")
    # cm.safe_save(adata, p.file_path)