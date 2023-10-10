import anndata as ad
from datetime import datetime
import pandas as pd
import utils.common as cm
from params import Params
import numpy as np

# pMN no -
# MALAT without 1
p = Params()
adata = ad.read_h5ad(p.file_path)
df = pd.read_csv(p.folder + "manual_annotations/exclude_dpt_meta.csv", index_col="index", comment="#")
adata.obs.drop(columns=df.columns.tolist(), inplace=True, errors='ignore')
adata.obs = adata.obs.merge(df, left_index=True, right_index=True, how='left')
cm.safe_save(adata, p.folder + "dpt.h5ad")