import scanpy as sc
from params import Params
import anndata as ad
import numpy as np
import utils.common as cm
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('macosx')

p = Params()
adata = ad.read_h5ad(p.folder + "ad_files/als_filtered.h5ad")

groups = pd.read_csv(p.folder + "meta_grouping.csv")
dfs = []
for index, row in groups.iterrows():
    obs_names = row['old_obs_names'].split('\t')
    new_df = pd.DataFrame({'obs_names': obs_names})
    new_df["SEACell"] = row["SEACell"]
    dfs.append(new_df)

processed_groups = pd.concat(dfs)
processed_groups.set_index("obs_names", inplace=True)
print("Merge processed_groups")
adata.obs = adata.obs.merge(processed_groups, left_index=True, right_index=True, how='left')

df = pd.read_csv(p.folder + "manual_annotations/chung.csv", index_col="index", comment="#")

print("Merge manual_annotations")
adata.obs = adata.obs.merge(df, left_on="SEACell", right_index=True, how='left')

cm.safe_save(adata, p.folder + "als_filtered_anno.h5ad")