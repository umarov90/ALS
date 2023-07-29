import anndata as ad
from datetime import datetime
from params import Params
import numpy as np

p = Params()
adata = ad.read_h5ad(p.file_path)
del adata.raw
adata.layers.clear()
# adata = adata[:, adata.var.highly_variable]
drop_columns = []
for column in adata.obs.columns:
    if adata.obs[column].nunique() == 1:
        drop_columns.append(column)
    if "marker" in column:
        drop_columns.append(column)
adata.obs.drop(drop_columns, axis=1, inplace=True)
# {datetime.now().strftime('%d.%m.%y')}
adata.write(p.folder + p.file_name.replace(".h5ad", "_vis_jul27.h5ad"))
