import anndata as ad
import os
from datetime import datetime
from params import Params
import numpy as np
import utils.common as cm
import scvelo as scv
import scanpy as sc
import pandas as pd
import seaborn as sns
import scanpy.external as sce
import joblib
import matplotlib
# matplotlib.use('macosx')


if __name__ == '__main__':
    p = Params()
    adata_m = ad.read_h5ad(p.file_path)
    adata_m.layers = joblib.load(p.folder + "layers_m.p")

    d7 = adata_m[adata_m.obs["day"] == "D7"]
    dd = d7.obs["sample_id"].unique()

    color_column = "manual_anno_L0"
    custom_palette = joblib.load(p.folder + color_column + "_palette.p")

    ids = []

    for sid in dd:
        ids.append(sid)
        print(sid)

        # if os.path.isfile(f"./figures/scvelo_meta_{sid}.png"):
        #     continue
        adata_d7 = adata_m[adata_m.obs["sample_id"] == sid]

        # adata_loom = scv.read(p.folder + f"/velocyto/{sid}.loom")
        # adata_loom.obs_names = adata_loom.obs_names.str.replace(":", "_").str.rstrip("x")
        #
        # groups = pd.read_csv(p.folder + "meta_grouping.csv")
        # dfs = []
        # for index, row in groups.iterrows():
        #     obs_names = row['old_obs_names'].split('\t')
        #     new_df = pd.DataFrame({'obs_names': obs_names})
        #     new_df["SEACell"] = row["SEACell"]
        #     dfs.append(new_df)
        #
        # processed_groups = pd.concat(dfs)
        # processed_groups.set_index("obs_names", inplace=True)
        # adata_loom.obs = adata_loom.obs.merge(processed_groups, left_index=True, right_index=True, how='left')
        #
        # df = pd.read_csv(p.folder + "manual_annotations/chung.csv", index_col="index", comment="#")
        #
        # adata_loom.obs = adata_loom.obs.merge(df, left_on="SEACell", right_index=True, how='left')
        #
        # adata_loom = adata_loom[~adata_loom.obs["SEACell"].isna()]
        #
        # scv.pp.filter_and_normalize(adata_loom, min_shared_counts=20, n_top_genes=2000)
        #
        # sc.tl.pca(adata_loom, use_highly_variable=True)
        # sc.pp.neighbors(adata_loom, n_neighbors=30, n_pcs=50)
        # sc.tl.umap(adata_loom)
        #
        # scv.pp.moments(adata_loom, n_neighbors=30)
        # scv.tl.velocity(adata_loom)
        # scv.tl.velocity_graph(adata_loom)
        # scv.pl.velocity_embedding_stream(adata_loom, basis='umap', color="manual_anno_L0", palette=custom_palette,
        #                                  save=f"original_{sid}.png", dpi=300)

        # Meta
        sc.tl.pca(adata_d7, use_highly_variable=True)
        sc.pp.neighbors(adata_d7, n_neighbors=30, n_pcs=50)
        sc.tl.umap(adata_d7)

        scv.pp.moments(adata_d7, n_neighbors=30)
        scv.tl.velocity(adata_d7)
        scv.tl.velocity_graph(adata_d7)
        scv.pl.velocity_embedding_stream(adata_d7, basis='umap', color=color_column,palette=custom_palette,
                                         save=f"meta_{sid}.png", dpi=300)

        print("Done")

    joblib.dump(ids, "ids.p")
