import anndata as ad
from datetime import datetime
from params import Params
import utils.common as cm
import numpy as np

p = Params()
adata = ad.read_h5ad(p.file_path)

adata.obs['treatment'] = adata.obs['treatment'].replace('na', 'Ctrl')

adata.obs["sample"] = adata.obs["day"].astype(str) + "_" + adata.obs["treatment"].astype(str)

cm.safe_save(adata, p.file_path)
