import anndata as ad
from params import Params

p = Params()
# adata1 = ad.read_h5ad(p.folder + p.file_name.replace(".h5ad", "_vis.h5ad"))
adata = ad.read_h5ad(p.file_path)
del adata.raw
adata = adata[:, adata.var.highly_variable]
adata.obs['individual'] = adata.obs['individual'].astype('category')
adata.write(p.folder + p.file_name.replace(".h5ad", "_vis2.h5ad"))