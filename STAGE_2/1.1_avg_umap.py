import SEACells
import anndata as ad
import pandas as pd
from params import Params
import matplotlib
matplotlib.use('macosx')
import scanpy as sc
import numpy as np
import utils.common as cm
params = Params()

adata_s = ad.read_h5ad(params.folder + "ad_files/als_filtered.h5ad")
adata_m = ad.read_h5ad(params.folder + "als_meta.h5ad")
groups = pd.read_csv(params.folder + "meta_grouping.csv")
adata_m.obsm[f'X_sc_umap'] = np.empty((len(adata_m), 2), dtype=np.float64)
name_to_row = {name: i for i, name in enumerate(adata_s.obs_names)}
for i, obs_name in enumerate(adata_m.obs_names):
    sublist = adata_m.obs.loc[obs_name, 'old_obs_names'].split("\t")
    rows = [name_to_row[name] for name in sublist if name in name_to_row]
    adata_m.obsm[f'X_sc_umap'][i, :] = adata_s.obsm[f'X_umap_harmony'][rows, :].mean(axis=0)

cm.safe_save(adata_m, params.file_path)
adata_m.obsm[f'X_umap'] = adata_m.obsm[f'X_sc_umap']
sc.pl.umap(adata_m, color="manual_anno_L0")
