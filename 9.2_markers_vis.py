import os
import scanpy as sc
import anndata as ad
from collections import Counter
from utils import utils
from utils.params import Params
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


FDR = 0.01
LOG_FOLD_CHANGE = 1.0


def ensure(folder):
    if not os.path.exists(folder):
        os.makedirs(folder)


louvain_cluster_col = "louvain"
params = Params()
adata = ad.read_h5ad(params.file_path)
adata.uns['log1p']["base"] = None
gene_dict = {}
for cluster in adata.obs[louvain_cluster_col].unique():
    gene_rank = sc.get.rank_genes_groups_df(adata, group=cluster, key='wilcoxon')[
        ['names', 'logfoldchanges', "pvals_adj", "pvals"]]
    gl = utils.process_gene_rank2(gene_rank, adata, 100)
    with open(f"out/{cluster}.tsv", 'w') as file:
        file.writelines('\n'.join(gl))
    gene_dict[cluster] = utils.process_gene_rank2(gene_rank, adata, 20)

all_values = [value for values in gene_dict.values() for value in values]
common_values = [value for value, count in Counter(all_values).most_common(40)]
common_values.remove("CRABP1")
# common_values.append("CRABP1")
sc.pl.dotplot(adata, common_values, groupby='louvain', dendrogram=True, use_raw=False)
plt.savefig('dotplot.svg')
sc.pl.heatmap(adata, common_values, groupby='louvain', dendrogram=True, use_raw=False)
plt.savefig('heatmap.svg')


