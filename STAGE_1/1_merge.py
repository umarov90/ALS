import re
import pandas as pd
import anndata as ad
import scanpy as sc
import utils.common as cm
from params import Params

p = Params()
meta_df = pd.read_csv(f'{p.folder}meta.tsv', sep='\t', index_col=False)
adatas = []
for index, row in meta_df.iterrows():
    print(f'Processing file {index + 1}')
    adata = sc.read_10x_h5(f'{p.folder}cellranger_output/{row["sample_id"]}/outs/filtered_feature_bc_matrix.h5')
    adata.var_names_make_unique()
    table = pd.read_csv(f'{p.folder}demuxlet/{row["sample_id"]}.best', sep='\t')
    adata.obs = adata.obs.merge(table[['BARCODE', 'SNG.BEST.GUESS', 'DROPLET.TYPE']].set_index('BARCODE'),
                                left_index=True, right_index=True, how='left')
    adata.obs.rename(columns={'SNG.BEST.GUESS': 'individual', 'DROPLET.TYPE': 'droplet_type'}, inplace=True)
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
adata.obs["sample"] = adata.obs["day"].astype(str) + "_" + adata.obs["treatment"].astype(str)
adata.obs['individual'] = adata.obs['individual'].astype(float).fillna(-1).astype(int).astype(str)
print(f"Number of individuals: {len(list(adata.obs['individual'].unique()))}")
cm.safe_save(adata, p.file_path)
