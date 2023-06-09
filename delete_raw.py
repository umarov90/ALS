import anndata as ad
from params import Params

p = Params()
adata = ad.read_h5ad(p.folder + "temp/als_1.h5ad")
del adata.raw
adata.write(p.file_path)