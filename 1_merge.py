import re
import pandas as pd
import anndata as ad
import scanpy as sc
from params import Params

p = Params()
chunk_size = 200000
df = pd.read_csv(f'{p.folder}meta.tsv', sep='\t', index_col=False)
adatas = []
for index, row in df.iterrows():
    print(f'Processing file {index + 1}')
    adata = sc.read_10x_h5(f'{p.folder}cellranger_output/{row["sample_id"]}/outs/filtered_feature_bc_matrix.h5')
    adata.var_names_make_unique()
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adata = adata[(adata.obs['n_genes_by_counts'] >= 300), :]

    table = pd.read_csv(f'{p.folder}demuxlet/{row["sample_id"]}.best', sep='\t')
    sng_barcodes = table[table['DROPLET.TYPE'].isin(['SNG'])]['BARCODE'].tolist()
    adata = adata[adata.obs_names.isin(sng_barcodes)]
    adata.obs = adata.obs.merge(table[['BARCODE', 'SNG.BEST.GUESS']].set_index('BARCODE'),
                                left_index=True, right_index=True)
    adata.obs.rename(columns={'SNG.BEST.GUESS': 'individual'}, inplace=True)
    info = re.split('_|-', row["sample_name"])
    adata.obs['pool'] = info[0]
    adata.obs['day'] = info[1]
    if len(info) == 3:
        adata.obs['treatment'] = "na"
        adata.obs['GEM'] = info[2]
    elif len(info) == 4:
        adata.obs['treatment'] = info[2]
        adata.obs['GEM'] = info[3]
    else:
        adata.obs['treatment'] = "na"
        adata.obs['GEM'] = "na"

    adata.obs['sample_id'] = row["sample_id"]
    adata.obs_names = row['sample_id'] + '_' + adata.obs_names.str.strip('-1')
    adatas.append(adata)
adata = ad.concat(adatas)
print("Concat")
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[(adata.obs['total_counts'] > 500) &
              (adata.obs['n_genes_by_counts'] > 300) &
              (adata.obs['pct_counts_mt'] < 10), :]
adata = adata[adata.obs['sample_id'].isin(
    adata.obs.value_counts('sample_id')[adata.obs.value_counts('sample_id') > 1000].index
), :]
print("Filtered")
# adata_n = ad.read_h5ad(p.file_path)
# adata_n.obs["treatment"] = adata.obs["treatment"]
# adata_n.obs["GEM"] = adata.obs["GEM"]
# adata_n.write(p.file_path)
print("Starting variable genes")
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    subset=False,
    flavor="seurat_v3",
    batch_key="pool", span=0.5
)
print("Finished variable genes")
adata.var.loc[adata.var['mt'], 'highly_variable'] = False
adata.var.loc[adata.var_names.str.match('^IG[HIKL]'), 'highly_variable'] = False
adata.raw = adata.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
print("Normalized")
sc.pp.log1p(adata)
adata = adata[:, adata.var.highly_variable]
print("PCA starting")
sc.tl.pca(adata, chunked=True, chunk_size=chunk_size)
print("PCA finished")
print("Saving")
adata.write(p.file_path)
print("All done")
