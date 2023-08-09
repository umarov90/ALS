import scanpy as sc
import anndata as ad
import utils.common as cm
from params import Params
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import numpy as np
# sns.set(font_scale=0.5)
params = Params()

adata = ad.read_h5ad(params.file_path)
tag = "manual"
for lvl in range(3):
    df = pd.DataFrame(columns=['clusters', 'names', 'scores'])
    column_name_level = "manual_anno_L" + str(lvl)
    cluster_list = adata.obs[column_name_level].unique().categories.tolist()
    adata_lvl = adata[adata.obs[column_name_level].isin(cluster_list)]
    lvl_folder = f"{params.folder}markers/L{lvl}/"
    for cluster in cluster_list:
        # all_markers = []
        # for nl in range(3):
        markers = adata_lvl.obs.loc[adata_lvl.obs[column_name_level] == cluster,
                                    "manual_anno_L" + str(lvl) + "_marker"].iloc[0]
        markers = markers[markers.index("-") + 1:]
        markers = markers.split(".")
        # all_markers.extend(markers)
        # markers = all_markers
        genes = []
        for m in markers:
            gene = m.split("_")[1]
            if gene in adata_lvl.var_names:
                genes.append(gene)
        with open(cm.ensure_file(f"{lvl_folder}tsv/{cluster}.tsv"), 'w') as file:
            file.write("\n".join(genes))
        try:
            scores = np.ravel(adata_lvl[:, genes].X.mean(axis=0))
            df_sub = pd.DataFrame({'names': genes, 'scores': scores})
            df_sub['clusters'] = cluster
            df = pd.concat([df, df_sub])
        except:
            print(cluster)
            print(genes)
    df.to_csv(cm.ensure_file(f"temp/{tag}_{lvl}_df.p"), index=False)
    df = df.sort_values(['names', 'scores'], ascending=False)
    df = df.drop_duplicates('names', keep='first')
    cluster_dict = df.groupby('clusters')['names'].apply(lambda x: x.tolist()).to_dict()
    all_genes = df['names'].unique()

    sc.pl.clustermap(adata_lvl[:, all_genes], obs_keys=column_name_level, use_raw=False,
                     save=f'_{tag}_L{lvl}.png')
    sc.pl.dotplot(adata_lvl, cluster_dict, groupby=column_name_level, dendrogram=True, use_raw=False,
                  swap_axes=True, save=f'_{tag}_L{lvl}.pdf')
    sc.pl.heatmap(adata_lvl, cluster_dict, groupby=column_name_level, dendrogram=True, use_raw=False,
                  show_gene_labels=True, swap_axes=True, save=f'_{tag}_L{lvl}.pdf')
    sc.pl.stacked_violin(adata_lvl, cluster_dict, groupby=column_name_level, swap_axes=True, dendrogram=True,
                         show_gene_labels=True, save=f'_{tag}_L{lvl}.pdf')

