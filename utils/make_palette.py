import anndata as ad
from datetime import datetime
import seaborn as sns
import joblib
from params import Params
import numpy as np

p = Params()
adata = ad.read_h5ad(p.file_path)

color_column = "manual_anno_L0"

cvals = adata.obs[color_column].unique()

custom_palette = {}
for i, c in enumerate(cvals):
    custom_palette[c] = sns.color_palette('Set1')[i]
joblib.dump(custom_palette, p.folder + color_column + "_palette.p")
