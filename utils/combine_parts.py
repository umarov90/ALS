import os
import scvelo as scv
from utils.params import Params
import anndata as ad
import matplotlib
matplotlib.use('macosx')
scv.set_figure_params()


p = Params()
directory = p.folder + "parts"
adatas = []
for filename in os.listdir(directory):
    filepath = os.path.join(directory, filename)
    if os.path.isfile(filepath):
        adata = ad.read_h5ad(filepath)
        adatas.append(adata)

adata = ad.concat(adatas)
adata.write(p.file_path)
