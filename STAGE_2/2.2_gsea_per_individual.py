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
# adata.obs["individual"] = adata.obs["individual"].astype(str) + "_" + adata.obs["gender"].astype(str)

groupby = "treatment"
vals = ["Ctrl", "Ropi"]
adata = adata[adata.obs["status"] == "ALS"].copy()

enr_results = []
for individual in adata.obs["individual"].unique():
    try:
        adata_i = adata[adata.obs["individual"] == individual].copy()
        sc.tl.rank_genes_groups(adata_i, groupby, method='wilcoxon', key_added="cluster_de", use_raw=False, tie_correct=True)
        # GSEA
        gene_rank = sc.get.rank_genes_groups_df(adata_i, group=vals[1], key='cluster_de')
        gene_rank = cm.process_gene_rank2(gene_rank, adata_i)
        gene_rank['names'] = gene_rank['names'].str.upper()
        gene_rank = gene_rank.loc[:, ["names", "scores"]]
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
# heatmap_data = joblib.load(f"hd_{groupby}.p")
plt.clf()

heatmap = sns.clustermap(
    heatmap_data,
    cmap='RdYlBu_r',
    cbar_kws={'label': 'NES'},
    yticklabels=True,
    xticklabels=True,
    figsize=(26, 12)
)

xtick_labels = heatmap.ax_heatmap.get_xticklabels()
for label in xtick_labels:
    individual = label.get_text()
    gender = adata.obs[adata.obs["individual"] == individual]["gender"].tolist()[0]
    if gender == "Male":
        label.set_color("blue")
    else:
        label.set_color("green")

# Set labels and title
heatmap.ax_heatmap.set_xlabel('Cluster')
heatmap.ax_heatmap.set_ylabel('Pathway')
heatmap.ax_heatmap.set_title('GSEA Heatmap ')
heatmap.ax_heatmap.set_facecolor('white')
# Rotate and align the X-axis labels
plt.setp(heatmap.ax_heatmap.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')

blue_patch = plt.Line2D([0], [0], marker='o', color='w', label='Male', markerfacecolor='blue', markersize=10)
red_patch = plt.Line2D([0], [0], marker='o', color='w', label='Female', markerfacecolor='green', markersize=10)
plt.legend(handles=[blue_patch, red_patch], loc='lower left')

plt.tight_layout()
# Save the figure to a PNG file
heatmap.savefig(f'individual_{groupby}_GSEA_heatmap.pdf', dpi=300)
