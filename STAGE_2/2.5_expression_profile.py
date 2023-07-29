import os
import scanpy as sc
import anndata as ad
import numpy as np
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

# Step 1: Count the number of rows for each individual
individual_counts = adata.obs["individual"].value_counts()
# Step 2: Identify individuals with less than 100 rows
individuals_to_drop = individual_counts[individual_counts < 100].index
# Step 3: Drop the rows for the identified individuals from adata
adata = adata[~adata.obs["individual"].isin(individuals_to_drop)]

groupby = ["treatment", "status"][0]
np.unique(adata.obs[groupby])
if groupby == "treatment":
    vals = ["Ctrl", "Ropi"]
    print(len(adata))
    adata = adata[adata.obs["status"] == "ALS"].copy()
    print(len(adata))
elif groupby == "status":
    vals = ["healthy", "ALS"]
    print(len(adata))
    adata = adata[adata.obs["treatment"] != "Ropi"].copy()
    print(len(adata))

for cluster in adata.obs[cluster_col].unique():
    try:
        adata_cluster = adata[adata.obs[cluster_col] == cluster].copy()
        sc.tl.rank_genes_groups(adata_cluster, groupby, method='wilcoxon', key_added=cluster_col + '_rank',
                                use_raw=False, tie_correct=True)
        gene_rank = sc.get.rank_genes_groups_df(adata, group=cluster, key=cluster_col + '_rank')
        gene_names = cm.process_gene_rank2(gene_rank, adata, cluster, 100)
        gene_expression = adata_cluster[:, gene_names].X
        average_expression = pd.DataFrame(index=adata_cluster.obs["individual"].unique(), columns=gene_names)

        for individual in adata_cluster.obs["individual"].unique():
            individual_indices = np.where(adata_cluster.obs["individual"] == individual)[0]
            average_expression.loc[individual] = np.mean(gene_expression[individual_indices], axis=0)

        for c in gene_names:
            average_expression[c] = average_expression[c] / average_expression[c].mean()
        individual_to_gender = dict(adata_cluster.obs[['individual', 'gender']].values)
        individual_to_status = dict(adata_cluster.obs[['individual', 'status']].values)
        average_expression.index = [f"{ind}_{individual_to_gender[ind]}_{individual_to_status[ind]}"
                                    for ind in average_expression.index]
        heatmap = sns.clustermap(average_expression.astype(float), cmap='viridis') # , yticklabels=True
        plt.setp(heatmap.ax_heatmap.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
        plt.tight_layout()
        heatmap.savefig(f'individual_de_{groupby}_{cluster}.png', dpi=300)
    except Exception as e:
        traceback.print_exc()
        pass
