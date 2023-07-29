import os
import joblib
import scanpy as sc
import anndata as ad
import numpy as np
import utils.common as cm
from params import Params
import pandas as pd
import matplotlib
matplotlib.use('macosx')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
# sns.set(font_scale=0.5)
louvain_cluster_col = "louvain"
params = Params()

adata = ad.read_h5ad(params.file_path)
adata.uns['log1p']["base"] = None
cluster_list = list(adata.obs[louvain_cluster_col].unique())
# Groupping to avoid picking genes that separates obvious groups
groups = {}
groups["IPS"] = [9, 11]
groups["Neurons"] = [15, 0, 1, 2, 5, 18, 16, 14, 17, 12, 22, 20, 4, 13, 21, 7]
named_clusters = []
for group_key in groups.keys():
    groups[group_key] = [str(num) for num in groups[group_key]]
    named_clusters.extend(groups[group_key])
groups["Rest"] = list(set(cluster_list) - set(named_clusters))

num_top_genes = 5
for group_key in groups.keys():
    tag = f"wilcoxon_{group_key}_{num_top_genes}"
    # adata_g = adata
    adata_g = adata[adata.obs["louvain"].isin(groups[group_key])].copy()
    cluster_list = list(adata_g.obs[louvain_cluster_col].unique())
    if os.path.exists(f"temp/{tag}_df.p"):
        df = pd.read_csv(f"temp/{tag}_df.p", dtype={'clusters': str})
    else:
        df = pd.DataFrame(columns=['clusters', 'names', 'scores'])
        for i in range(len(cluster_list)):
            for j in range(i+1, len(cluster_list)):
                if cluster_list[i] != cluster_list[j]:
                    adata_p = adata_g[adata_g.obs[louvain_cluster_col].isin([cluster_list[i], cluster_list[j]])]
                    sc.tl.rank_genes_groups(adata_p, "louvain", method='wilcoxon', key_added="rank", use_raw=False)
                    # t-test_overestim_var
                    gene_rank = sc.get.rank_genes_groups_df(adata_p, group=cluster_list[i], key='rank')
                    genes1 = cm.process_gene_rank2(gene_rank, adata_p, cluster_list[i], n=num_top_genes)
                    df = pd.concat([df, genes1])
                    gene_rank = sc.get.rank_genes_groups_df(adata_p, group=cluster_list[j], key='rank')
                    genes2 = cm.process_gene_rank2(gene_rank, adata_p, cluster_list[j], n=num_top_genes)
                    df = pd.concat([df, genes2])
                    print(f"{cluster_list[i]} {cluster_list[j]}")

        df.to_csv(cm.ensure_file(f"temp/{tag}_df.p"), index=False)

    df = df.sort_values(['names', 'expression'], ascending=False)
    df = df.drop_duplicates('names', keep='first')

    df2 = df.sort_values(['clusters', 'scores'], ascending=False)
    cluster_genes = df2.groupby('clusters')['names'].apply(lambda x: x.tolist()).to_dict()
    cluster_values = df2.groupby('clusters')['scores'].apply(lambda x: x.tolist()).to_dict()
    fig, axs = plt.subplots(len(cluster_genes), 1, figsize=(10, len(cluster_genes) * 5), sharex=True)

    # Iterate over each cluster and subplot
    for i, (cluster, genes) in enumerate(cluster_genes.items()):
        values = cluster_values[cluster]

        genes = genes[:10]
        values = values[:10]

        ax = axs[i]  # Select the current subplot
        # Plot the genes and values for the cluster
        x = np.arange(len(genes))
        ax.plot(x, values, marker='o')
        # Add cluster names as text on the plot
        for j, (gene, value) in enumerate(zip(genes, values)):
            ax.text(j, value, gene, ha='center', va='bottom')
        # Set the x-tick labels as gene names
        ax.set_xticks(x)
        ax.set_xticklabels(genes)
        # Set the subplot title
        ax.set_title(cluster)

    plt.tight_layout()
    plt.savefig(f'figures/{tag}_gene_rank.pdf', bbox_inches='tight')

    cluster_dict = df.groupby('clusters')['names'].apply(lambda x: x.tolist()).to_dict()
    all_values = df['names'].unique()

    for cluster in cluster_dict.keys():
        cluster_genes = list(set(cluster_dict[cluster]))
        with open(cm.ensure_file(f"DGE/{cluster}.txt"), 'w') as file:
            file.writelines('\n'.join(cluster_genes))
        cluster_signature = ",".join(cluster_genes)
        print(f"{cluster}\t{cluster_signature}")

    sc.pl.clustermap(adata_g[:, all_values], obs_keys='louvain', use_raw=False, save=f'_{tag}.png')
    sc.pl.dotplot(adata_g, cluster_dict, groupby='louvain', dendrogram=True, use_raw=False,
                  swap_axes=True, save=f'{tag}.pdf')
    sc.pl.heatmap(adata_g, cluster_dict, groupby='louvain', dendrogram=True, use_raw=False,
                  show_gene_labels=True, swap_axes=True, save=f'_{tag}.pdf')

    corr_matrix = adata.X.T.corr()
    df = pd.DataFrame(corr_matrix, index=adata.var_names, columns=adata.var_names)
    clust_map = sns.clustermap(df, cmap='RdBu_r', figsize=(10, 10))

    row_linkage = clust_map.row_linkage
    column_linkage = clust_map.col_linkage

    # from scipy.spatial import distance
    # from scipy.cluster import hierarchy
    #
    # correlations_array = np.asarray(df.corr())
    # row_linkage = hierarchy.linkage(
    #     distance.pdist(correlations_array), method='average')
    # col_linkage = hierarchy.linkage(
    #     distance.pdist(correlations_array.T), method='average')

    # Use the reordered row and column indices to plot the clustermap
    sc.pl.clustermap(adata, row_linkage=row_linkage, col_linkage=column_linkage, cmap='RdBu_r')

