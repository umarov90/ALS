import os
import scanpy as sc
import anndata as ad
import numpy as np
from params import Params
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy
import gseapy.plot as gplt
from gseapy import gseaplot, heatmap


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
    plt.title(title)
    folder = f"figs/dge/{column}/volcano/"
    ensure(folder)
    fig.savefig(f'{folder}{title}.png', dpi=300)
    plt.close(fig)


p = Params()
adata = ad.read_h5ad(p.file_path)
adata = adata[adata.obs["day"] == "D14"].copy()
column = ["treatment", "status"][1]
vals = np.unique(adata.obs[column])

old_adata = ad.read_h5ad(p.folder + "als_1.h5ad").raw.to_adata()
filtered_raw = old_adata[adata.obs_names].copy()
adata.raw = filtered_raw

adata = adata.raw.to_adata()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

louvain_cluster_col = 'new_pca_sub_louvain_0.5'
unique_pools = np.unique(adata.obs['pool'])
for pool in unique_pools:
    adata_p = adata[adata.obs["pool"] == pool].copy()
    enr_results = []
    cluster_labels = []
    for cluster in adata_p.obs[louvain_cluster_col].unique():
        adata_cluster = adata_p[adata_p.obs[louvain_cluster_col] == cluster].copy()
        if len(adata_cluster[adata_cluster.obs[column] == vals[0]]) < 500:
            continue
        if len(adata_cluster[adata_cluster.obs[column] == vals[1]]) < 500:
            continue
        # status_values = adata_cluster.obs['status'].values
        # np.random.shuffle(status_values)
        # adata_cluster.obs['status'] = status_values
        sc.tl.rank_genes_groups(adata_cluster, column, method='wilcoxon', key_added="cluster_de", use_raw=False)
        for key in vals:
            volcano_plot(adata_cluster, key, f"{pool}_{cluster}_{key}")
        try:
            # Gene Set Analysis
            gene_set_names = gseapy.get_library_name(organism='Human')
            glist = sc.get.rank_genes_groups_df(adata_cluster, group='ALS',
                                                key='cluster_de', log2fc_min=0.25,
                                                pval_cutoff=0.05)['names'].squeeze().str.strip().tolist()
            enr_res = gseapy.enrichr(gene_list=glist,
                                     organism='Human',
                                     gene_sets='GO_Biological_Process_2018',
                                     cutoff=0.5)
            folder_GO = f"figs/dge/{column}/GO/"
            ensure(folder_GO)
            gseapy.barplot(enr_res.res2d, title='GO_Biological_Process_2018',
                           figsize=(6, 10), ofname=f"{folder_GO}{pool}_{cluster}_ALS_GO.png")
        except:
            print(f"{pool}_{cluster} Gene Set Analysis not performed")
        try:
            # GSEA
            gene_rank = sc.get.rank_genes_groups_df(adata_cluster, group='ALS', key='cluster_de')[['names', 'logfoldchanges']]
            gene_rank.sort_values(by=['logfoldchanges'], inplace=True, ascending=False)
            # calculate_qc_metrics will calculate number of cells per gene
            sc.pp.calculate_qc_metrics(adata_cluster, percent_top=None, log1p=False, inplace=True)
            # filter for genes expressed in at least 30 cells.
            gene_rank = gene_rank[gene_rank['names'].isin(adata_cluster.var_names[adata_cluster.var.n_cells_by_counts > 30])]
            res = gseapy.prerank(rnk=gene_rank, gene_sets='KEGG_2021_Human')
            terms = res.res2d.Term
            folder_out = f"figs/dge/{column}/GSEA/gseaplot/"
            ensure(folder_out)
            gseapy.gseaplot(figsize=(6, 10), ofname=f"{folder_out}{pool}_{cluster}_ALS_KEGG.png",
                            rank_metric=res.ranking, term=terms[0], **res.results[terms[0]])
            res.res2d['Cluster'] = cluster
            enr_results.append(res.res2d)
            cluster_labels.append(cluster)
        except:
            print(f"{pool}_{cluster} GSEA not performed")
    # Concatenate the GSEA results for each cluster into a single DataFrame
    enr_combined = pd.concat(enr_results)

    # Pivot the DataFrame to prepare for heatmap plotting
    heatmap_data = enr_combined.pivot(index='Term', columns='Cluster', values='NES')
    heatmap_data = heatmap_data.fillna(0)
    top_pathways = heatmap_data.abs().mean(axis=1).nlargest(30).index
    heatmap_data = heatmap_data.loc[top_pathways]

    plt.clf()
    # Create the heatmap
    fig, ax = plt.subplots(figsize=(12, 10))

    # Set the color map
    cmap = 'RdYlBu_r'

    # Plot the heatmap using seaborn
    heatmap = sns.heatmap(
        heatmap_data,
        cmap=cmap,
        ax=ax,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={'label': 'NES'}  # Add color bar label
    )

    # Set labels and title
    ax.set_xlabel('Cluster')
    ax.set_ylabel('Pathway')
    ax.set_title('GSEA Heatmap ' + pool)

    # Rotate and align the X-axis labels
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
    plt.tight_layout()
    # Save the figure to a PNG file
    folder_out = f"figs/dge/{column}/GSEA/heatmap/"
    ensure(folder_out)
    fig.savefig(f'{folder_out}{pool}_ALS_GSEA.png', dpi=300)
