import os
import scanpy.external as sce
import anndata as ad
import matplotlib
import numpy as np
import scanpy as sc
from utils import utils
from utils.params import Params
matplotlib.use('Agg')

FDR = 0.01
LOG_FOLD_CHANGE = 1.0


def ensure(folder):
    if not os.path.exists(folder):
        os.makedirs(folder)


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

adata.obsm[f'X_cluster_dge_sub_umap'] = np.empty((len(adata), 2), dtype=np.float64)

for cluster in adata.obs[louvain_cluster_col].unique():
    print(cluster)
    try:
        adata_cluster = adata[adata.obs[louvain_cluster_col] == cluster].copy()
        if len(adata_cluster[adata_cluster.obs[groupby] == vals[0]]) < 100:
            continue
        if len(adata_cluster[adata_cluster.obs[groupby] == vals[1]]) < 100:
            continue
        # status_values = adata_cluster.obs['status'].values
        # np.random.shuffle(status_values)
        # adata_cluster.obs['status'] = status_values
        sc.tl.rank_genes_groups(adata_cluster, groupby, method='wilcoxon', key_added="cluster_de", use_raw=False,
                                tie_correct=True)

        gene_rank = sc.get.rank_genes_groups_df(adata_cluster, group=vals[1], key='cluster_de')[
            ['names', 'logfoldchanges', "pvals_adj", "pvals"]]
        gene_rank = utils.process_gene_rank(gene_rank, adata_cluster)
        gene_rank = gene_rank.iloc[gene_rank['slogp'].abs().nlargest(200).index]
        gene_names = gene_rank['names']

        adata_cluster = adata_cluster[:, gene_names]

        sc.tl.pca(adata)
        sce.pp.harmony_integrate(adata, 'pool', max_iter_harmony=100, max_iter_kmeans=200)
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50, use_rep="X_pca_harmony")

        sc.tl.umap(adata_cluster)
        adata.obsm[f'X_cluster_dge_sub_umap'][adata.obs[louvain_cluster_col] == cluster] = adata_cluster.obsm['X_umap']
    except:
        pass

adata.write(params.file_path)
