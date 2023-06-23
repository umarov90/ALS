import sys
import scvelo as scv
from utils.params import Params
import scanpy as sc
import anndata as ad
import numpy as np
from scipy.sparse import csr_matrix
import matplotlib
matplotlib.use('macosx')
scv.set_figure_params()



p = Params()
sample_id = "TFHS000394" # sys.argv[1] #
adata = ad.read_h5ad(p.file_path)
adata = adata[adata.obs["sample_id"] == sample_id].copy()
ldata = scv.read(p.folder + "looms/" + sample_id + ".loom", cache=True)
ldata.obs_names = ldata.obs_names.str.replace(":", "_")
ldata.obs_names = ldata.obs_names.str[:-1]

mats = ldata.layers
m = len(adata)

new_mats = mats.copy()
for key in new_mats:
    new_mats[key] = csr_matrix((m, mats[key].shape[1]), dtype=np.float64)

name_to_row = {name: i for i, name in enumerate(ldata.obs_names)}
for i, obs_name in enumerate(adata.obs_names):
    sublist = adata.obs.loc[obs_name, 'old_obs_names']
    for mat_keys in mats.keys():
        rows = [name_to_row[name] for name in sublist if name in name_to_row]
        if rows:
            for key in new_mats:
                new_mats[key][i, :] = mats[key][rows, :].sum(axis=0)

adata.layers = new_mats
adata.write(p.file_path)