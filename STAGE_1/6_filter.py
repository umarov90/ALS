import scanpy as sc
from params import Params
import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('macosx')

p = Params()
adata = ad.read_h5ad(p.file_path)
print(len(adata))
adata = adata[adata.obs["droplet_type"] == "SNG"]
print(len(adata))
adata = adata[adata.obs["wrong_pool"] == False]
print(len(adata))

adata = adata[adata.obs["predicted_doublet"] == False]
print(len(adata))
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

adata = adata[(adata.obs['total_counts'] > 500) &
              (adata.obs['n_genes_by_counts'] > 300) &
              (adata.obs['pct_counts_mt'] < 10), :]
print(len(adata))
adata.write(p.folder + p.file_name.replace(".h5ad", "_filtered.h5ad"))

