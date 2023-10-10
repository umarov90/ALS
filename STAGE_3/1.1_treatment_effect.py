import anndata as ad
import numpy as np
from params import Params
from scipy import stats
import matplotlib
matplotlib.use('Agg')


cluster_col = "manual_anno_L0"
params = Params()
adata = ad.read_h5ad(params.folder + "ad_files/als_filtered_anno.h5ad")
adata.uns['log1p']["base"] = None
adata = adata[adata.obs["day"].isin(["D14"])].copy()
# adata = adata[adata.obs["gender"] == "Female"].copy()
individual_dict = {}
individuals = adata.obs['individual'].unique()
healthy_individuals = adata[adata.obs["status"] == "healthy"].obs['individual'].unique().tolist()
als_individuals = adata[adata.obs["status"] == "ALS"].obs['individual'].unique().tolist()
h_list_r = []
h_list_c = []
a_list_r = []
a_list_c = []
for individual in individuals:
    print(individual)
    adata_i = adata[adata.obs['individual'] == individual]
    adata_r = adata_i[adata_i.obs['treatment'] == "Ropi"]
    adata_c = adata_i[adata_i.obs['treatment'] != "Ropi"]
    if len(adata_r) > 0 and len(adata_c) > 0:
        a = len(adata_r[adata_r.obs[cluster_col] == "mMN"]) / len(adata_r)
        b = len(adata_c[adata_c.obs[cluster_col] == "mMN"]) / len(adata_c)
        if a > 0.1 and b > 0.1:
            individual_dict[individual] = a / b
            if individual in healthy_individuals:
                h_list_r.append(a)
                h_list_c.append(b)
            else:
                a_list_r.append(a)
                a_list_c.append(b)

ratio = []
for i in range(len(a_list_r)):
    ratio.append(a_list_r[i] / a_list_r[i])

# Perform the t-test
t_statistic, p_value = stats.ttest_ind(a_list_r, a_list_c)
print("ALS T-statistic:", t_statistic)
print("ALS P-value:", p_value)

t_statistic, p_value = stats.ttest_ind(h_list_r, h_list_c)
print("Healthy T-statistic:", t_statistic)
print("Healthy P-value:", p_value)
