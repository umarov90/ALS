import SEACells
import anndata as ad
import pandas as pd
from params import Params
import utils.common as cm
params = Params()

adata = ad.read_h5ad(params.folder + "ad_files/als_filtered.h5ad")
groups = pd.read_csv(params.folder + "meta_grouping.csv")

# Create a seacell id to old obs names mapping df
dfs = []
for index, row in groups.iterrows():
    obs_names = row['old_obs_names'].split('\t')
    new_df = pd.DataFrame({'obs_names': obs_names})
    new_df["SEACell"] = row["SEACell"]
    dfs.append(new_df)

processed_groups = pd.concat(dfs)
processed_groups.set_index("obs_names", inplace=True)

adata.obs = adata.obs.merge(processed_groups, left_index=True, right_index=True, how='left')
# Convert the adata to metacell level by summarize_by_SEACell
SEACell_ad = SEACells.core.summarize_by_SEACell(adata, SEACells_label='SEACell', summarize_layer='raw')

# Restore the obs columns by averaging, sum or picking first value
summary = adata.obs
for col in summary.columns:
    if summary[col].dtype == 'float64':
        if "total" in col:
            summary[col] = summary.groupby('SEACell')[col].transform('sum')
        else:
            summary[col] = summary.groupby('SEACell')[col].transform('mean')
    else:
        summary[col] = summary.groupby('SEACell')[col].transform('first')

summary = summary.drop_duplicates('SEACell', keep='first')
summary = summary.set_index('SEACell')
SEACell_ad.obs = SEACell_ad.obs.merge(summary, left_index=True, right_index=True, how='left')
# Add old_obs_names column to obs
groups = groups.set_index("SEACell")
SEACell_ad.obs = SEACell_ad.obs.merge(groups, left_index=True, right_index=True, how='left')
cm.safe_save(SEACell_ad, params.file_path)
