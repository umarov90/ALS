import multiprocessing
multiprocessing.set_start_method("fork")
import sys
import scvelo as scv
from params import Params
import anndata as ad
import joblib
import pandas as pd
import scanpy as sc
from utils import common as cm
import numpy as np
from scipy.sparse import csr_matrix
import matplotlib
matplotlib.use('macosx')


params = Params()
color_column = "manual_anno_L0"
custom_palette = joblib.load(params.folder + color_column + "_palette.p")
adata = ad.read_h5ad(params.file_path)
gene_order = adata.var_names.tolist()
with open(cm.ensure_file(f"{params.folder}pseudobulk/gene_order.tsv"), 'w') as file:
    file.write("\n".join(gene_order))

for lvl in range(3):
    print(f"Level {lvl}")
    column_name_level = "manual_anno_L" + str(lvl)
    cluster_list = adata.obs[column_name_level].unique().categories.tolist()
    adata_lvl = adata[adata.obs[column_name_level].isin(cluster_list)].copy()
    for cluster in cluster_list:
        adata_cluster = adata_lvl[adata_lvl.obs[column_name_level] == cluster].copy()
        print(f"Cluster {cluster}. Size {len(adata_cluster)}.")
        for day in adata_cluster.obs["day"].unique():
            adata_d = adata_cluster[adata_cluster.obs["day"] == day].copy()
            pseudobulk_list = []
            meta = ["Individual\tGender\tStatus\tDay\tmanual_anno_L0\tmanual_anno_L1\tmanual_anno_L2"] # \tTreatment\tPool

            # for pool in adata_cluster.obs["pool"].unique():
            #     adata_p = adata_cluster[adata_cluster.obs["pool"] == pool].copy()
            adata_p = adata_d

            for individual in adata_p.obs["individual"].unique():
                adata_i = adata_p[adata_p.obs["individual"] == individual].copy()

                # for treatment in adata_i.obs["treatment"].unique():
                #     adata_t = adata_i[adata_i.obs["treatment"] == treatment].copy()
                adata_t = adata_i

                celltype_expression = adata_t.X.sum(axis=0)
                pseudobulk_list.append(celltype_expression)
                meta.append(individual + "\t" + adata_t.obs["gender"].tolist()[0]
                            + "\t" + adata_t.obs["status"].tolist()[0] + "\t" + day + "\t"
                            # + treatment + "\t" + adata_t.obs["pool"].tolist()[0]  + "\t"
                            + adata_t.obs["manual_anno_L0"].tolist()[0] + "\t"
                            + adata_t.obs["manual_anno_L1"].tolist()[0] + "\t"
                            + adata_t.obs["manual_anno_L2"].tolist()[0])
            pseudobulk = np.squeeze(np.asarray(pseudobulk_list)).T
            print(f"Shape {pseudobulk.shape}, len {len(meta)}")
            np.savetxt(cm.ensure_file(f"{params.folder}pseudobulk/lvl_{lvl}/{cluster}_{day}.tsv"), pseudobulk, delimiter='\t', fmt='%d')
            with open(cm.ensure_file(f"{params.folder}pseudobulk/meta_lvl_{lvl}/{cluster}_{day}.tsv"), 'w') as file:
                file.write("\n".join(meta))




