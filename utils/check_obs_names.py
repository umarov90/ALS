import anndata as ad
from datetime import datetime
from params import Params
import numpy as np
import utils.common as cm

p = Params()

adata2 = ad.read_h5ad(p.file_path) # als_meta_vis0.6.h5ad
# a = adata2.obs["old_obs_names"]
# a.to_csv(p.folder + "meta_grouping.csv")
adata1 = ad.read_h5ad(p.folder + "als_combined.h5ad")

if set(adata1.obs_names) == set(adata2.obs_names):
    print("obs_names match")

common_genes = np.intersect1d(adata1.var_names, adata2.var_names)[:100]
adata1_common_genes = adata1[:, common_genes].X
adata2_common_genes = adata2.raw[:, common_genes].X
if np.allclose(adata1_common_genes.toarray(), adata2_common_genes.toarray()):
    print("The values of the common genes match")

if set(adata1.obs["old_obs_names"]) == set(adata1.obs["old_obs_names"]):
    print("old_obs_names match")

print("Done")