import os
from params import Params
import anndata as ad
import utils.common as cm

p = Params()
directory = p.folder + "parts"
adatas = []
for filename in os.listdir(directory):
    filepath = os.path.join(directory, filename)
    if os.path.isfile(filepath):
        try:
            adata = ad.read_h5ad(filepath)
            adatas.append(adata)
        except:
            pass
print(f"Merging {len(adatas)} parts")
adata = ad.concat(adatas)
cm.safe_save(adata, p.folder + "als_combined.h5ad")
