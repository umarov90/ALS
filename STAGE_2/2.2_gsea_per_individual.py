import os
import scanpy as sc
import anndata as ad
import numpy as np
from utils import common as cm
from params import Params
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy
import traceback
import joblib


cluster_col = "manual_anno_L0"
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

enr_results = []
for individual in adata.obs["individual"].unique():
    adata_i = adata[adata.obs["individual"] == individual].copy()
    sc.tl.rank_genes_groups(adata_i, groupby, method='wilcoxon', key_added="cluster_de", use_raw=False,
                            tie_correct=True)
    try:
        # GSEA
        gene_rank = sc.get.rank_genes_groups_df(adata_i, group=vals[1], key='cluster_de')[
            ['names', 'logfoldchanges', "pvals_adj", "pvals"]]
        gene_rank = cm.process_gene_rank(gene_rank, adata_i)
        gene_rank['names'] = gene_rank['names'].str.upper()
        res = gseapy.prerank(rnk=gene_rank, gene_sets="../data/GMTs/c2.cp.kegg.v2023.1.Hs.symbols.gmt")  # 'KEGG_2021_Human'
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
top_pathways = heatmap_data.abs().max(axis=1).nlargest(50).index
heatmap_data = heatmap_data.loc[top_pathways]
joblib.dump(heatmap_data, f"hd_{groupby}.p")
plt.clf()

heatmap = sns.clustermap(
    heatmap_data,
    cmap='RdYlBu_r',
    cbar_kws={'label': 'NES'},
    yticklabels=True,
    figsize=(16, 12)
)

# Set labels and title
heatmap.ax_heatmap.set_xlabel('Cluster')
heatmap.ax_heatmap.set_ylabel('Pathway')
heatmap.ax_heatmap.set_title('GSEA Heatmap ')

# Rotate and align the X-axis labels
plt.setp(heatmap.ax_heatmap.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
plt.tight_layout()
# Save the figure to a PNG file
heatmap.savefig(f'individual_GSEA_heatmap.png', dpi=300)
