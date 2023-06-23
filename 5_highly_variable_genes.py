import anndata as ad
import scanpy as sc
import numpy as np
from utils.params import Params

p = Params()
adata = ad.read_h5ad(p.file_path)
adata.raw = adata.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.var['mt'] = adata.var_names.str.startswith('MT-')
print("Starting variable genes")
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    subset=False,
    batch_key="pool", span=0.5
)
print("Finished variable genes")
adata.var.loc[adata.var['mt'], 'highly_variable'] = False
adata.var.loc[adata.var_names.str.match('^IG[HIKL]'), 'highly_variable'] = False
adata.write(p.file_path)

