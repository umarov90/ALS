import os
import scanpy as sc
import anndata as ad
import numpy as np

import utils
from utils.params import Params
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

FDR = 0.01
LOG_FOLD_CHANGE = 1.0


def ensure(folder):
    if not os.path.exists(folder):
        os.makedirs(folder)


def volcano_plot(adata, group_key, title=""):
    result = sc.get.rank_genes_groups_df(adata, group=group_key, key="cluster_de")
    result["-logQ"] = -np.log(result["pvals"].astype("float"))
    lowqval_de = result.loc[abs(result["logfoldchanges"]) > LOG_FOLD_CHANGE]
    other_de = result.loc[abs(result["logfoldchanges"]) <= LOG_FOLD_CHANGE]

    fig, ax = plt.subplots()
    sns.regplot(
        x=other_de["logfoldchanges"],
        y=other_de["-logQ"],
        fit_reg=False,
        scatter_kws={"s": 6},
    )
    sns.regplot(
        x=lowqval_de["logfoldchanges"],
        y=lowqval_de["-logQ"],
        fit_reg=False,
        scatter_kws={"s": 6},
    )
    ax.set_xlabel("log2 FC")
    ax.set_ylabel("-log Q-value")
    ax.set_xlim(-5, 5)
    plt.title(title)
    folder = f"figs/dge/{groupby}/volcano/"
    ensure(folder)
    fig.savefig(f'{folder}{title}.png', dpi=300)
    plt.close(fig)


louvain_cluster_col = "louvain"
chunk_size = 10000
params = Params()
adata = ad.read_h5ad(params.file_path)
adata.uns['log1p']["base"] = None
groupby = ["treatment", "status"][1]
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

for cluster in adata.obs[louvain_cluster_col].unique():
    try:
        adata_cluster = adata[adata.obs[louvain_cluster_col] == cluster].copy()
        if len(adata_cluster[adata_cluster.obs[groupby] == vals[0]]) < 100:
            continue
        if len(adata_cluster[adata_cluster.obs[groupby] == vals[1]]) < 100:
            continue
        # status_values = adata_cluster.obs['status'].values
        # np.random.shuffle(status_values)
        # adata_cluster.obs['status'] = status_values
        sc.tl.rank_genes_groups(adata_cluster, groupby, method='wilcoxon', key_added="cluster_de", use_raw=False, tie_correct=True)

        gene_rank = sc.get.rank_genes_groups_df(adata_cluster, group=vals[1], key='cluster_de')[['names', 'logfoldchanges', "pvals_adj", "pvals"]]
        gene_rank = utils.process_gene_rank(gene_rank, adata_cluster)
        gene_rank = gene_rank.iloc[gene_rank['slogp'].abs().nlargest(200).index]
        gene_names = gene_rank['names']

        gene_expression = adata_cluster[:, gene_names].X
        average_expression = pd.DataFrame(index=adata_cluster.obs["individual"].unique(), columns=gene_names)

        for individual in adata_cluster.obs["individual"].unique():
            individual_indices = np.where(adata_cluster.obs["individual"] == individual)[0]
            average_expression.loc[individual] = np.mean(gene_expression[individual_indices], axis=0)

        heatmap = sns.clustermap(average_expression.astype(float), cmap='viridis')
        plt.setp(heatmap.ax_heatmap.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
        plt.tight_layout()
        folder_out = f"figs/dge/{groupby}/heatmap/"
        ensure(folder_out)
        heatmap.savefig(f'{folder_out}individual_dge_cluster_{cluster}.png', dpi=300)
    except:
        pass
