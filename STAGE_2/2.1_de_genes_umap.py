import os
import traceback

import scanpy.external as sce
import anndata as ad
import matplotlib
import numpy as np
import scanpy as sc
import utils.common as cm
import pandas as pd
from params import Params
matplotlib.use('Agg')

params = Params()
FDR = 0.01
LOG_FOLD_CHANGE = 1.0
num_top_genes = 50
cluster_col = "manual_anno_L0"
adata = ad.read_h5ad(params.file_path)
adata.uns['log1p']["base"] = None
groupby = ["treatment", "status"][1]
np.unique(adata.obs[groupby])

adata.obsm[f'X_{groupby}_de_sub_umap'] = np.empty((len(adata), 2), dtype=np.float64)

for cluster in adata.obs[cluster_col].unique():
    print(cluster)
    try:

        adata_cluster = adata[adata.obs[cluster_col] == cluster].copy()
        if groupby == "treatment":
            vals = ["Ctrl", "Ropi"]
            adata_cluster = adata_cluster[adata_cluster.obs["status"] == "ALS"].copy()
        elif groupby == "status":
            vals = ["healthy", "ALS"]
            adata_cluster = adata_cluster[adata_cluster.obs["treatment"] != "Ropi"].copy()
        adata_val0 = adata_cluster[adata_cluster.obs[groupby] == vals[0]]
        adata_val1 = adata_cluster[adata_cluster.obs[groupby] == vals[1]]
        print(len(adata_cluster))
        df = pd.DataFrame(columns=['clusters', 'names', 'scores'])
        sc.tl.rank_genes_groups(adata_cluster, groupby, method='wilcoxon',
                                key_added="rank", use_raw=False, tie_correct=True)

        cm.volcano_plot(adata_cluster, vals[0], title=f"{cluster} {vals[0]}", save=f"figures/{cluster}_{vals[0]}.pdf")
        cm.volcano_plot(adata_cluster, vals[1], title=f"{cluster} {vals[1]}", save=f"figures/{cluster}_{vals[1]}.pdf")

        # gene_rank = sc.get.rank_genes_groups_df(adata_cluster, group=vals[0], key='rank')
        # genes1 = cm.process_gene_rank2(gene_rank, adata_val0,
        #                                vals[0], n=num_top_genes)
        # df = pd.concat([df, genes1])
        # gene_rank = sc.get.rank_genes_groups_df(adata_cluster, group=vals[1], key='rank')
        # genes2 = cm.process_gene_rank2(gene_rank, adata_val1, vals[1], n=num_top_genes)
        # df = pd.concat([df, genes2])
        #
        # adata_cluster = adata_cluster[:, df["names"].unique()]
        #
        # sc.tl.pca(adata_cluster, n_comps=10, use_highly_variable=False)
        # sce.pp.harmony_integrate(adata_cluster, 'pool', max_iter_harmony=10,
        #                          max_iter_kmeans=20, block_size=0.01, sigma=0.2)
        # sc.pp.neighbors(adata_cluster, n_neighbors=10, use_rep="X_pca_harmony")
        #
        # sc.tl.umap(adata_cluster)
        # cm.copy_obsm(adata, adata_cluster, f'X_{groupby}_de_sub_umap', "X_umap")
        # print(f"Finished {cluster}")

    except:
        traceback.print_exc()
        print(f"Failed {cluster}!!!")

# cm.safe_save(adata, params.file_path)
