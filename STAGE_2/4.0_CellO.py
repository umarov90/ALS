import cello
import scvelo as scv
from params import Params
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use('macosx')
scv.set_figure_params()


if __name__ == '__main__':
    p = Params()
    adata = ad.read_h5ad(p.file_path)
    # sc.pp.normalize_total(adata, target_sum=1e6)
    # sc.pp.log1p(adata)
    cello_resource_loc = "/Users/ramzan"
    model_prefix = "ALS_MODEL"

    cello.scanpy_cello(
        adata,
        'louvain',
        cello_resource_loc,
        out_prefix=model_prefix
    )
    sc.tl.umap(adata)
    adata.write(p.file_path)
    sc.pl.umap(adata, color='louvain')
    sc.pl.umap(adata, color='Most specific cell type')
    cello.cello_probs(adata, '0', cello_resource_loc, 0.2, clust_key='louvain')

