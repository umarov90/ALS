# Finding wrong individual labels based on pool
import gc
import anndata as ad
import utils.common as cm
from params import Params

p = Params()
adata = ad.read_h5ad(p.file_path)
print(f"Number of individuals: {len(list(adata.obs['individual'].unique()))}")
print(len(adata))
pool_size = 30
pool_num = 5
adata.obs['status'] = 'NA'
pools = {}
for i in range(pool_num):
    for j in range(pool_size):
        geno_id = str(i * pool_size + (j + 1))
        print(geno_id, end=" ")
        # 22 to 30 are healthy individuals
        if j + 1 > 22:
            adata.obs.loc[adata.obs['individual'] == geno_id, 'status'] = 'healthy'
            print("control")
        else:
            adata.obs.loc[adata.obs['individual'] == geno_id, 'status'] = 'ALS'
            print("ALS")
        pools.setdefault(str(i+1), []).append(geno_id)

print(adata.obs['status'].value_counts())

adata.obs['wrong_pool'] = True
adatas = []
unique_pools = list(adata.obs['pool'].unique())
for pool_id in unique_pools:
    print(pool_id)
    valid_individuals = pools[pool_id[-1]] + pools[pool_id[-2]]
    adata_p = adata[adata.obs["pool"] == pool_id].copy()
    print(adata_p.obs['wrong_pool'].value_counts())
    # (adata_p.obs['pool'] != pool_id) |
    adata_p.obs['wrong_pool'][adata_p.obs['individual'].isin(valid_individuals)] = False
    print(adata_p.obs['wrong_pool'].value_counts())
    adatas.append(adata_p)

del adata
gc.collect()

adata = ad.concat(adatas)
print(adata.obs['wrong_pool'].value_counts())
cm.safe_save(adata, p.file_path)


