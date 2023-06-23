import scanpy as sc
import anndata as ad
from utils.params import Params

p = Params()
adata = ad.read_h5ad(p.file_path)
adata.uns['log1p']["base"] = None
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50, use_rep="X_pca_harmony")
sc.tl.louvain(adata, key_added="louvain", resolution=0.6)
# sc.tl.rank_genes_groups(adata, "louvain", method='wilcoxon', key_added="wilcoxon", use_raw=False,
#                         tie_correct=True)
adata.write(p.file_path)