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

# Remove individuals with less than 100 cells
individual_counts = adata.obs["individual"].value_counts()
individuals_to_drop = individual_counts[individual_counts < 100].index
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

# for gender in adata.obs["gender"].unique():
#     adata_g = adata[adata.obs["gender"] == gender].copy()
for cluster in adata.obs[cluster_col].unique():
    try:
        adata_cluster = adata[adata.obs[cluster_col] == cluster].copy()

        sc.tl.rank_genes_groups(adata_cluster, groupby, method='wilcoxon',
                                key_added="rank", use_raw=False, tie_correct=True)
        gene_names = []
        gene_rank = sc.get.rank_genes_groups_df(adata_cluster, group=vals[0], key='rank')
        gene_names.extend(cm.process_gene_rank2(gene_rank, adata_cluster, n=25, filter_low=False))

        gene_rank = sc.get.rank_genes_groups_df(adata_cluster, group=vals[1], key='rank')
        gene_names.extend(cm.process_gene_rank2(gene_rank, adata_cluster, n=25, filter_low=False))

        SEACell_ad = SEACells.core.summarize_by_SEACell(adata_cluster, SEACells_label='individual', summarize_layer='raw')

        SEACell_ad.raw = SEACell_ad
        sc.pp.normalize_total(SEACell_ad, target_sum=1e4)
        sc.pp.log1p(SEACell_ad)

        individual_to_groupby = adata_cluster.obs.groupby('individual')[groupby].first().to_dict()
        SEACell_ad.obs[groupby] = SEACell_ad.obs_names.map(individual_to_groupby)
        # sc.pl.heatmap(SEACell_ad, gene_names, groupby=groupby,  use_raw=False, # standard_scale="var",
        #               save=f'individual_de_{groupby}_{cluster}.pdf',  cmap=sns.color_palette("rocket", as_cmap=True))
        sc.pl.clustermap(SEACell_ad[:, gene_names], obs_keys=groupby, use_raw=False,
                      save=f'individual_de_{groupby}_{cluster}.pdf')
    except Exception as e:
        traceback.print_exc()
        pass
