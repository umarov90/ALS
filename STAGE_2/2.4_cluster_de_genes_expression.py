import scanpy as sc
import anndata as ad
import utils.common as cm
from params import Params
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import SEACells
import pandas as pd
# sns.set(font_scale=0.5)

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
    gene_rank = cm.process_gene_rank2(gene_rank, adata[adata.obs[cluster_col] == cluster])
    gene_dict[cluster] = gene_rank

all_values = [value for values in gene_dict.values() for value in values.nlargest(10, 'scores')['names'].tolist()]
gene_names = list(set(all_values))
# joblib.dump(all_values, "genes.p")
# all_values = joblib.load("genes.p")
SEACell_ad = SEACells.core.summarize_by_SEACell(adata, SEACells_label=cluster_col, summarize_layer='raw')

SEACell_ad.raw = SEACell_ad
sc.pp.normalize_total(SEACell_ad, target_sum=1e4)
sc.pp.log1p(SEACell_ad)
groupby = cluster_col
to_groupby = adata.obs.groupby(groupby)[groupby].first().to_dict()
SEACell_ad.obs[groupby] = SEACell_ad.obs_names.map(to_groupby)
SEACell_ad.obs[groupby] = SEACell_ad.obs[groupby].astype('category')
sc.pl.heatmap(SEACell_ad, gene_names, groupby=groupby, standard_scale="var", use_raw=True,
              save=f'de_{groupby}.pdf',  cmap=sns.color_palette("rocket", as_cmap=True))
sc.pl.clustermap(SEACell_ad[:, gene_names], obs_keys=groupby, use_raw=False,
                      save=f'de_{groupby}.pdf')

