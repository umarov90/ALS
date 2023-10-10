import anndata as ad
import numpy as np
import scanpy as sc
import scanpy.external as sce
from params import Params
import pandas as pd
import utils.common as cm
import matplotlib
matplotlib.use('macosx')

p = Params()
rn = 0
adata = ad.read_h5ad(p.file_path)
adata_small = adata[adata.obs["manual_anno_L0"].isin(['NSC'])].copy()
# adata = adata[adata.obs["gender"] == "Female"].copy()
# genes = pd.read_csv(p.folder + "classification_genes2.txt")["gene"].tolist()
genes = pd.read_csv(p.folder + "classification_genes4.txt")["gene"].tolist()
# genes = list(set(genes + genes2))
adata_small = adata_small[:, genes]
# sc.tl.pca(adata, use_highly_variable=False, random_state=rn, n_comps=50)
sc.pp.neighbors(adata_small, n_neighbors=10, use_rep="X")
sc.tl.umap(adata_small, random_state=rn)
sc.pl.umap(adata_small, color="day")
sc.pl.umap(adata_small, color="status")

# adata.obsm['X_umap_status_classification'] = np.zeros((len(adata), 2), dtype=np.float64)
# adata.obsm['X_umap_status_classification'][adata.obs["manual_anno_L0"].isin(['iPSC', 'pMN', 'NSC', 'oMN'])] = adata_small.obsm['X_umap']
# cm.safe_save(adata, p.file_path)