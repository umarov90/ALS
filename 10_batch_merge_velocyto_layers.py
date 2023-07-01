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


if __name__ == '__main__':
    p = Params()
    wd = sys.argv[1]
    sample_id = sys.argv[2]
    adata = ad.read_h5ad(p.file_path)
    adata = adata[adata.obs["sample_id"] == sample_id].copy()
    ldata = scv.read(wd + "/" + sample_id + "/velocyto/" + sample_id + ".loom", cache=True)
    ldata.obs_names = ldata.obs_names.str.replace(":", "_")
    ldata.obs_names = ldata.obs_names.str[:-1]

    mats = ldata.layers
    m = len(adata)

    for key in mats:
        adata.layers[key] = csr_matrix((m, mats[key].shape[1]), dtype=np.float64)

    name_to_row = {name: i for i, name in enumerate(ldata.obs_names)}
    for i, obs_name in enumerate(adata.obs_names):
        if i % 50 == 0:
            print(i, end=" ")
        sublist = adata.obs.loc[obs_name, 'old_obs_names'].split("\t")
        rows = [name_to_row[name] for name in sublist if name in name_to_row]
        for key in mats:
            adata.layers[key][i, :] = mats[key][rows, :].sum(axis=0)

    adata.write(p.folder + "parts/" + sample_id + ".h5ad")
