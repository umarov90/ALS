import multiprocessing as mp
import os
from pathlib import Path

import joblib
import pylab as p
import scanpy as sc
import shutil
import numpy as np
from params import Params
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
params = Params()
import gc
import scvelo as scv
import scanpy as sc
import anndata as ad
import pandas as pd
from scipy.sparse import csr_matrix
# matplotlib.use('macosx')
scv.set_figure_params()


def load_layers(amats, mats, ldata_obs_names, name_to_row, id):
    for i, obs_name in enumerate(ldata_obs_names):
        if i % 500 == 0:
            print(i, end=" ")
        for key in ["spliced", "unspliced"]:
            if obs_name in name_to_row.keys():
                amats[key][name_to_row[obs_name], :] += mats[key][i, :]
    joblib.dump(amats, params.folder + "layer_parts/" + id + ".p")


def process_gene_rank(gene_rank, adata, n=-1):
    lowly_expressed = adata.var.index[np.ravel(adata.X.mean(axis=0)) < 0.05]
    exclude_genes = pd.read_csv(params.folder + "chrY_genes.csv")["gene"].tolist()
    exclude_genes.append("XIST")
    exclude_genes.append("MYC")
    gene_rank = gene_rank[~gene_rank['names'].isin(lowly_expressed)]
    gene_rank = gene_rank[~gene_rank['names'].isin(exclude_genes)]

    gene_rank["slogp"] = -np.log(gene_rank["pvals"].astype("float"))
    gene_rank['slogp'] *= np.sign(gene_rank['logfoldchanges'])
    gene_rank = gene_rank[gene_rank['pvals_adj'] <= 0.05]
    gene_rank.drop('logfoldchanges', axis=1, inplace=True)
    gene_rank.drop('pvals_adj', axis=1, inplace=True)
    gene_rank.drop('pvals', axis=1, inplace=True)
    # gene_rank = gene_rank[gene_rank['slogp'] > 0]
    gene_rank.sort_values(by=['slogp'], inplace=True, ascending=False)
    # calculate_qc_metrics will calculate number of cells per gene
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    # filter for genes expressed in at least 30 cells.
    gene_rank = gene_rank[gene_rank['names'].isin(adata.var_names[adata.var.n_cells_by_counts > 30])]
    gene_rank = gene_rank.reset_index(drop=True)
    if n > 0:
        return gene_rank['slogp'].abs().nlargest(n).index.tolist()
    else:
        return gene_rank


def process_gene_rank2(gene_rank, adata_cluster, cluster, n=-1):
    gene_rank = gene_rank[gene_rank['pvals_adj'] <= 0.05]

    sc.pp.calculate_qc_metrics(adata_cluster, percent_top=None, log1p=False, inplace=True)
    gene_rank = gene_rank[gene_rank['names'].isin(adata_cluster.var_names[adata_cluster.var.n_cells_by_counts > len(adata_cluster) / 2])]

    lowly_expressed = adata_cluster.var.index[np.ravel(adata_cluster.X.mean(axis=0)) < 0.05]
    exclude_genes = pd.read_csv(params.folder + "chrY_genes.csv")["gene"].tolist()
    exclude_genes.append("XIST")
    exclude_genes.append("MYC")

    gene_rank = gene_rank[~gene_rank['names'].isin(lowly_expressed)]
    gene_rank = gene_rank[~gene_rank['names'].isin(exclude_genes)]
    gene_rank = gene_rank[~gene_rank['names'].str.startswith("MT")]
    gene_rank["expression"] = np.ravel(adata_cluster[:, gene_rank["names"].tolist()].X.mean(axis=0))

    # gene_rank["scores"] = -np.log(gene_rank["pvals"].astype("float"))
    # gene_rank['scores'] *= np.sign(gene_rank['logfoldchanges'])
    # gene_rank['scores'] = gene_rank['logfoldchanges']
    gene_rank.reset_index(inplace=True, drop=True)
    if n > 0:
        # gene_rank = gene_rank.nlargest(n, 'scores')
        # gene_rank["clusters"] = cluster
        gene_rank = gene_rank.iloc[gene_rank['scores'].abs().nlargest(n).index]['names']
    return gene_rank



def ensure(path):
    Path(path).mkdir(parents=True, exist_ok=True)
    return path

def ensure_file(path):
    Path(os.path.dirname(path)).mkdir(parents=True, exist_ok=True)
    return path

def safe_save(adata, name):
    temp_name = name.replace(".h5ad", "_temp.h5ad")
    adata.write(temp_name)
    shutil.move(temp_name, name)


def copy_obsm(adata, adata_cluster, param1, param2):
    obs_names_list = adata_cluster.obs_names.tolist()
    n = 0
    for i, obs_name in enumerate(adata.obs_names):
        if obs_name in adata_cluster.obs_names:
            adata.obsm[param1][i] = adata_cluster.obsm[param2][obs_names_list.index(obs_name)]
            n += 1
    print(f"Big one {len(adata)} Small one {len(adata_cluster)} Copied {n}")


def volcano_plot(adata, group_key, title="", save=""):
    FDR = 0.01
    LOG_FOLD_CHANGE = 1.0
    result = sc.get.rank_genes_groups_df(adata, group=group_key, key="rank")

    exclude_genes = pd.read_csv(params.folder + "chrY_genes.csv")["gene"].tolist()
    exclude_genes.append("XIST")
    exclude_genes.append("MYC")
    result = result[~result['names'].isin(exclude_genes)]

    result["-logQ"] = -np.log(result["pvals"].astype("float"))

    max_non_inf_score = result.loc[result['-logQ'] != np.inf, 'scores'].max()
    max_non_inf_logQ = result.loc[result['-logQ'] != np.inf, '-logQ'].max()
    result['diff_score'] = np.abs(result['scores'] / max_non_inf_score)
    result.loc[result['-logQ'] == np.inf, '-logQ'] = \
        max_non_inf_logQ * result.loc[result['-logQ'] == np.inf, 'diff_score']

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

    sorted_data = lowqval_de.sort_values(by="-logQ", ascending=False)

    top_positive = sorted_data[sorted_data['logfoldchanges'] > 0].head(4)
    top_negative = sorted_data[sorted_data['logfoldchanges'] < 0].head(4)
    top_data = pd.concat([top_positive, top_negative])

    for i, (x, y, text) in enumerate(zip(top_data["logfoldchanges"], top_data["-logQ"], top_data["names"])):
        va = 'bottom' if i % 2 == 0 else 'top'
        y_adj = y + 0.5 if va == 'bottom' else y - 0.5
        plt.text(x, y_adj, text, fontsize=8, ha='left', va=va, rotation=0)

    ax.set_xlabel("log2 FC")
    ax.set_ylabel("-log Q-value")
    ax.set_xlim(-5, 5)
    plt.title(title)
    if len(save) > 0:
        fig.savefig(save, dpi=300)
    plt.close(fig)


def prepare_adata_layers(adata):
    # desired_order = adata.var_names.tolist()
    # print("Adata loaded")
    df = pd.read_csv(f'{params.folder}meta.tsv', sep='\t', index_col=False)
    # ps = []
    # for index, row in df.iterrows():
    #     ldata = scv.read(params.folder + f"/velocyto/{row['sample_id']}.loom")
    #     ldata.obs_names = ldata.obs_names.str.replace(":", "_").str.rstrip("x")
    #     ldata.var_names_make_unique()
    #     ldata = ldata[:, desired_order]
    #     # Copy layers from ldata to adata in parallel
    #     ldata_obs_names = ldata.obs_names.tolist()
    #     mats = ldata.layers.copy()
    #     mats.pop("matrix")
    #     mats.pop("ambiguous")
    #     del ldata
    #     gc.collect()
    #
    #     adata_p = adata[adata.obs["sample_id"] == row['sample_id']]
    #     if "old_obs_names" in adata_p.obs:
    #         name_to_row = {}
    #         for i, obs_name in enumerate(adata_p.obs_names):
    #             sublist = adata_p.obs.loc[obs_name, 'old_obs_names'].split("\t")
    #             for n in sublist:
    #                 name_to_row[n] = i
    #     else:
    #         name_to_row = {name: i for i, name in enumerate(adata_p.obs_names)}
    #
    #     amats = {}
    #     for key in ["spliced", "unspliced"]:
    #         amats[key] = csr_matrix((adata_p.shape[0], adata_p.shape[1]), dtype=np.float64)
    #     del adata_p
    #     gc.collect()
    #     load_proc = mp.Process(target=load_layers,
    #                            args=(amats, mats, ldata_obs_names, name_to_row,row['sample_id'],))
    #     load_proc.start()
    #     ps.append(load_proc)
    #     print(f"File {index + 1}")
    #     if len(ps) > 9:
    #         for load_proc in ps:
    #             load_proc.join()
    #         ps = []
    # for load_proc in ps:
    #     load_proc.join()
    for key in ["spliced", "unspliced"]:
        adata.layers[key] = csr_matrix((adata.shape[0], adata.shape[1]), dtype=np.float64)
    for index, row in df.iterrows():
        print(row['sample_id'])
        layer_part = joblib.load(params.folder + "layer_parts/" + row['sample_id'] + ".p")
        for key in ["spliced", "unspliced"]:
            adata.layers[key][adata.obs["sample_id"] == row['sample_id']] = layer_part[key]
    return adata
