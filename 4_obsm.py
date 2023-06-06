import anndata as ad
import scanpy as sc
import scanpy.external as sce
from params import Params

p = Params()
adata = ad.read_h5ad(p.file_path)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
sc.tl.umap(adata)
adata.obsm['X_umap_original'] = adata.obsm['X_umap']

sce.pp.harmony_integrate(adata, 'pool')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50, use_rep="X_pca_harmony")
sc.tl.umap(adata)
adata.obsm['X_umap_harmony'] = adata.obsm['X_umap']

del adata.obsm["X_umap"]
adata.write(p.file_path)