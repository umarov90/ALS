import anndata as ad
from datetime import datetime
from params import Params
import utils.common as cm
import numpy as np

p = Params()
adata = ad.read_h5ad(p.file_path)

drop_columns = []
for column in adata.obs.columns:
    if adata.obs[column].nunique() == 1:
        drop_columns.append(column)
    if "marker" in column:
        drop_columns.append(column)
    if column.endswith("_x") or column.endswith("_y"):
        drop_columns.append(column)
adata.obs.drop(drop_columns, axis=1, inplace=True)
cm.safe_save(adata, p.file_path)
