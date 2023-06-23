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
import gseapy
import traceback
import joblib


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

# for cluster in adata.obs[louvain_cluster_col].unique():
#     adata_cluster = adata[adata.obs[louvain_cluster_col] == cluster].copy()
enr_results = []
for individual in adata.obs["individual"].unique():
    adata_i = adata[adata.obs["individual"] == individual].copy()
    if len(adata_i[adata_i.obs[groupby] == vals[0]]) < 10:
        continue
    if len(adata_i[adata_i.obs[groupby] == vals[1]]) < 10:
        continue
    # status_values = adata_i.obs['status'].values
    # np.random.shuffle(status_values)
    # adata_i.obs['status'] = status_values
    sc.tl.rank_genes_groups(adata_i, groupby, method='wilcoxon', key_added="cluster_de", use_raw=False,
                            tie_correct=True)
    # for key in vals:
    #     volcano_plot(adata_i, key, f"all_pools_{0}_{key}")
    try:
        # GSEA
        gene_rank = sc.get.rank_genes_groups_df(adata_i, group=vals[1], key='cluster_de')[
            ['names', 'logfoldchanges', "pvals_adj", "pvals"]]
        gene_rank = utils.process_gene_rank(gene_rank, adata_i)
        gene_rank['names'] = gene_rank['names'].str.upper()
        res = gseapy.prerank(rnk=gene_rank, gene_sets="c2.cp.kegg.v2023.1.Hs.symbols.gmt")  # 'KEGG_2021_Human'
        # terms = res.res2d.Term
        # folder_out = f"figs/dge/{groupby}/GSEA/gseaplot/"
        # ensure(folder_out)
        # gseapy.gseaplot(figsize=(6, 10), ofname=f"{folder_out}all_pools_{cluster}_{vals[1]}_KEGG.png",
        #                 rank_metric=res.ranking, term=terms[0], **res.results[terms[0]])
        res.res2d['Individual'] = individual
        enr_results.append(res.res2d)
    except Exception as e:
        traceback.print_exc()
        # print(f"{cluster} GSEA not performed")
# Concatenate the GSEA results for each cluster into a single DataFrame
enr_combined = pd.concat(enr_results)

# Pivot the DataFrame to prepare for heatmap plotting
heatmap_data = enr_combined.pivot(index='Term', columns='Individual', values='NES')
heatmap_data = heatmap_data.fillna(0)
top_pathways = heatmap_data.abs().max(axis=1).nlargest(50).index # min p adjust 0.05, NES max at least 1 passing NES 2
heatmap_data = heatmap_data.loc[top_pathways]
joblib.dump(heatmap_data, f"hd_{groupby}.p")
plt.clf()

heatmap = sns.clustermap(
    heatmap_data,
    cmap='RdYlBu_r',
    cbar_kws={'label': 'NES'},
    yticklabels=True,
    figsize=(12, 12)
)

# Set labels and title
heatmap.ax_heatmap.set_xlabel('Cluster')
heatmap.ax_heatmap.set_ylabel('Pathway')
heatmap.ax_heatmap.set_title('GSEA Heatmap ')

# Rotate and align the X-axis labels
plt.setp(heatmap.ax_heatmap.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
plt.tight_layout()
# Save the figure to a PNG file
folder_out = f"figs/dge/{groupby}/GSEA/heatmap/"
ensure(folder_out)
heatmap.savefig(f'{folder_out}cluster_{0}_{vals[1]}_GSEA.png', dpi=300)
