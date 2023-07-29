import scanpy as sc
from params import Params
import anndata as ad
import numpy as np
import utils.common as cm
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('macosx')

p = Params()
adata_old = ad.read_h5ad(p.folder + "ad_files/als.h5ad")
adata = ad.read_h5ad(p.file_path)

avg_scores = []
for i, row in adata.obs.iterrows():
    names = row["old_obs_names"].split("\t")
    avg_scores.append(adata_old.obs.loc[names, 'doublet_score'].mean())

adata.obs["doublet_score"] = avg_scores
cm.safe_save(adata, p)
