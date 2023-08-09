# import multiprocessing
# multiprocessing.set_start_method("fork")
import joblib
import scvelo as scv
from params import Params
import anndata as ad
import scanpy.external as sce
import scanpy as sc
import pandas as pd
from utils import common as cm
# matplotlib.use('macosx')
scv.set_figure_params()


if __name__ == '__main__':
    p = Params()
    color_column = "manual_anno_L0"
    custom_palette = joblib.load(p.folder + color_column + "_palette.p")
    adata = ad.read_h5ad(p.folder + "ad_files/als_filtered.h5ad")
    adata = adata[adata.obs["day"] != "D0"]
    del adata.raw
    adata.obsm["X_pca"] = adata.obsm["X_pca_harmony"].copy()
    del adata.obsm["X_pca_harmony"]
    # adata = cm.prepare_adata_layers(adata)
    # print("Layers copied")
    # joblib.dump(adata.layers, p.folder + "layers_s.p")
    # print("Layers dumped")
    print("Layers loading")
    adata.layers = joblib.load(p.folder + "layers_s.p")
    print("Layers loaded")

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

    # scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    #
    # sc.tl.pca(adata, use_highly_variable=True)
    #
    # sce.pp.harmony_integrate(adata, 'pool', max_iter_harmony=10,
    #                          max_iter_kmeans=20, block_size=0.01, sigma=0.2)
    # sc.pp.neighbors(adata, n_neighbors=30, n_pcs=50, use_rep="X_pca_harmony")
    # sc.tl.umap(adata)
    print("moments")
    # sc.tl.pca(adata, use_highly_variable=True, n_comps=30)
    scv.pp.moments(adata)
    scv.tl.velocity(adata, use_highly_variable=True)
    scv.tl.velocity_graph(adata, approx=True)
    scv.pl.velocity_embedding_stream(adata, basis='X_umap_harmony', color=color_column, palette=custom_palette,
                                     save="sc_ALL.png", dpi=300)
    adata.write(p.folder + "als_sc_velo.h5ad")
