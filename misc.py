import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
from params import Params

p = Params()
adata = ad.read_h5ad(p.file_path)
adata.obs["sample"] = adata.obs["day"].astype(str) + "_" + adata.obs["treatment"].astype(str)
adata.write(p.file_path)