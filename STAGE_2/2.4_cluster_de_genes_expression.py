import scanpy as sc
import anndata as ad
import utils.common as cm
from params import Params
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
sns.set(font_scale=0.5)

cluster_col = "manual_anno_L0"
params = Params()

adata = ad.read_h5ad(params.file_path)
adata.uns['log1p']["base"] = None
# sc.tl.rank_genes_groups(adata, cluster_col, method='wilcoxon',
#                         key_added="manual_anno_L0_rank", use_raw=False, tie_correct=True)
# cm.safe_save(adata, params.file_path)
adata.uns['log1p']["base"] = None
gene_dict = {}
for cluster in adata.obs[cluster_col].unique():
    gene_rank = sc.get.rank_genes_groups_df(adata, group=cluster, key=cluster_col + '_rank')
    # gl = cm.process_gene_rank2(gene_rank, adata_t, cluster)
    # with open(f"out/{neuron_val}_{cluster}.tsv", 'w') as file:
    #     file.writelines('\n'.join(gl))
    gene_rank = cm.process_gene_rank2(gene_rank, adata, cluster)
    gene_dict[cluster] = gene_rank

all_values = [value for values in gene_dict.values() for value in values.nlargest(20, 'scores')['names'].tolist()]
all_values = list(set(all_values))
# joblib.dump(all_values, "genes.p")
# all_values = joblib.load("genes.p")
exp_list = []
for cluster in adata.obs[cluster_col].unique():
    adata_cluster = adata[adata.obs[cluster_col] == cluster].copy()
    gene_expression = adata_cluster[:, all_values].X.mean(axis=0)
    exp_list.append(gene_expression)
avg_expression = np.concatenate(exp_list, axis=0)
df = pd.DataFrame(avg_expression, index=adata.obs[cluster_col].unique(), columns=all_values)
fig, ax = plt.subplots(figsize=(80,10))
heatmap = sns.clustermap(df, cmap='viridis', xticklabels=1, yticklabels=1, figsize=(20,5))
plt.setp(heatmap.ax_heatmap.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
fig = plt.gcf()
fig.set_size_inches(10, fig.get_size_inches()[1])
heatmap.savefig("DE_expression_clustermap.pdf", dpi=300)


