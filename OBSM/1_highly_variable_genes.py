import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
from params import Params
import utils.common as cm

p = Params()
adata = ad.read_h5ad(p.file_path)

exclude_genes = pd.read_csv(p.folder + "chrY_genes.csv")["gene"].tolist()
exclude_genes.append("XIST")
exclude_genes.append("MYC")

adata.var['mt'] = adata.var_names.str.startswith('MT-')
print("Starting variable genes")
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2200,
    subset=False,
    flavor="seurat_v3",
    batch_key="pool", span=1.0
)
print("Finished variable genes")
adata.var.loc[adata.var['mt'], 'highly_variable'] = False
adata.var.loc[adata.var_names.str.match('^IG[HIKL]'), 'highly_variable'] = False
adata.var.loc[adata.var_names.isin(exclude_genes), 'highly_variable'] = False

adata.raw = adata.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# sc.pp.scale(adata, max_value=10)

cm.safe_save(adata, p.file_path)

