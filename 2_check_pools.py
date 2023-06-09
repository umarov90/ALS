import gc
import anndata as ad
from params import Params

p = Params()
adata = ad.read_h5ad(p.file_path)
print(len(adata))
pool_size = 30
pool_num = 5
adata.obs['status'] = 'ALS'
pools = {}
for i in range(pool_num):
    for j in range(pool_size):
        geno_id = i * pool_size + (j + 1)
        print(geno_id, end=" ")
        if j + 1 > 22:
            adata.obs.loc[adata.obs['individual'] == geno_id, 'status'] = 'healthy'
            print("control", end=" ")
        print()
        pools.setdefault(str(i+1), []).append(geno_id)

chunk_size = 100000
chunks = [adata[i:i+chunk_size].copy() for i in range(0, len(adata), chunk_size)]
unique_pools = list(adata.obs['pool'].unique())
del adata
gc.collect()
for pool_id in unique_pools:
    valid_individuals = pools[pool_id[-1]] + pools[pool_id[-2]]
    for i, chunk in enumerate(chunks):
        start_n = len(chunks[i])
        chunks[i] = chunks[i][(chunks[i].obs['pool'] != pool_id) | (chunks[i].obs['individual'].isin(valid_individuals))].copy()
        gc.collect()
        print(f"{pool_id} removed {start_n - len(chunks[i])}")

adata = ad.concat(chunks)
print(len(adata))
adata.write(p.file_path)


