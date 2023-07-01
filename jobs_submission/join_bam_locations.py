import pandas as pd
from params import Params

p = Params()
df = pd.read_csv(f'{p.folder}meta.tsv', sep='\t', index_col=False)
unique_sample_ids = df['sample_id'].unique()

hca = pd.read_csv("../utils/hca.tsv", sep="\t", index_col=False)
hca["server"] = "HCA"
d7 = pd.read_csv("../utils/d7.tsv", sep="\t", index_col=False)
d7["server"] = "d7"
both = pd.concat([hca, d7])
both.drop_duplicates(subset='sample', inplace=True)
both = both[both['sample'].isin(unique_sample_ids)]
both.to_csv("both.tsv", sep="\t", index=False)