import os
import scanpy as sc
import anndata as ad
import numpy as np
import SEACells
import utils.common as cm
from params import Params
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import traceback


cluster_col = "manual_anno_L0"
params = Params()
adata = ad.read_h5ad(params.file_path)
adata.uns['log1p']["base"] = None
adata = adata[adata.obs["status"] == "healthy"].copy()
cluster_dict = {}
individuals = adata.obs['individual'].unique()
for cluster in adata.obs[cluster_col].unique():
    print(cluster)
    cluster_dict[cluster] = {}
    adata_cluster = adata[adata.obs[cluster_col] == cluster].copy()
    for individual in individuals:
        cluster_dict[cluster][individual] = (adata_cluster.obs['individual'] == individual).sum() / len(adata_cluster)

df = pd.DataFrame.from_dict(cluster_dict)
sns.set(font_scale=0.8)
clustermap = sns.clustermap(df, cmap="rocket", linewidths=0, figsize=(10, 40), yticklabels=True)

# Customize labels and title if needed
clustermap.ax_heatmap.set_xlabel('Individuals')
clustermap.ax_heatmap.set_ylabel('Clusters')
clustermap.fig.suptitle('Cluster Map of Individuals in Clusters', fontsize=16)

plt.savefig("clustermap1.png", bbox_inches='tight', dpi=300)

