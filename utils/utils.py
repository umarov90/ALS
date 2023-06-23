import scanpy as sc
import numpy as np
import matplotlib
matplotlib.use('Agg')


def process_gene_rank(gene_rank, adata, n=-1):
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


def process_gene_rank2(gene_rank, adata, n=-1):
    gene_rank = gene_rank[gene_rank['pvals_adj'] <= 0.05]
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    gene_rank = gene_rank[gene_rank['names'].isin(adata.var_names[adata.var.n_cells_by_counts > 30])]
    if n > 0:
        return gene_rank.nlargest(n, 'logfoldchanges')['names'].tolist()
    else:
        return gene_rank