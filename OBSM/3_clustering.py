import scanpy as sc
import anndata as ad
from params import Params
import utils.common as cm

p = Params()
adata = ad.read_h5ad(p.file_path)
adata.uns['log1p']["base"] = None
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50, use_rep="X_pca_harmony")
sc.tl.louvain(adata, key_added="louvain", resolution=1.0)
# sc.tl.rank_genes_groups(adata, "louvain", method='wilcoxon', key_added="wilcoxon", use_raw=False, tie_correct=True)
cm.safe_save(adata, p.file_path)
