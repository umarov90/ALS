import joblib
import scvelo as scv
from params import Params
import scanpy as sc
import pandas as pd
import anndata as ad
import scanpy.external as sce
import matplotlib
# matplotlib.use('macosx')
scv.set_figure_params()

if __name__ == '__main__':
    p = Params()
    adata = ad.read_h5ad(p.folder + "als_filtered.h5ad")
    print("loading layers")
    adata.layers = joblib.load(p.folder + "layers_s.p")
    print("loaded")
    # adata = ad.read_h5ad(p.file_path)
    # adata.layers = joblib.load(p.folder + "layers_m.p")

    groups = pd.read_csv(p.folder + "meta_grouping.csv")
    dfs = []
    for index, row in groups.iterrows():
        obs_names = row['old_obs_names'].split('\t')
        new_df = pd.DataFrame({'obs_names': obs_names})
        new_df["SEACell"] = row["SEACell"]
        dfs.append(new_df)

    processed_groups = pd.concat(dfs)
    processed_groups.set_index("obs_names", inplace=True)
    print("Merge processed_groups")
    adata.obs = adata.obs.merge(processed_groups, left_index=True, right_index=True, how='left')

    df = pd.read_csv(p.folder + "manual_annotations/chung.csv", index_col="index", comment="#")

    print("Merge manual_annotations")
    adata.obs = adata.obs.merge(df, left_on="SEACell", right_index=True, how='left')

    print("Drop without SEACell")
    adata = adata[~adata.obs["SEACell"].isna()]

    color_column = "manual_anno_L0"
    custom_palette = joblib.load(p.folder + color_column + "_palette.p")

    sub_info = "day"
    info_vals = adata.obs[sub_info].unique()
    for info_val in info_vals:
        print(f"{sub_info} {info_val}")
        adata_d = adata[adata.obs[sub_info] == info_val].copy()
        sc.tl.pca(adata_d)
        # sce.pp.harmony_integrate(adata_d, 'pool', max_iter_harmony=10,
        #                          max_iter_kmeans=20, block_size=0.01, sigma=0.2)
        sc.pp.neighbors(adata_d, n_neighbors=30, n_pcs=30) # , use_rep="X_pca_harmony"
        sc.tl.umap(adata_d)
        # sc.pl.umap(adata_d, color="manual_anno_L0", save=f"{sub_info}_{info_val}.png")

        scv.pp.moments(adata_d, n_neighbors=30) # , use_rep="X_pca_harmony"
        scv.tl.velocity(adata_d)
        scv.tl.velocity_graph(adata_d)
        scv.pl.velocity_embedding_stream(adata_d, basis='X_umap',
                                         save=f"{sub_info}_{info_val}.png",
                                         color=color_column, palette=custom_palette)


