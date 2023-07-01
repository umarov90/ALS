import math
import numpy as np
import SEACells
import anndata as ad
import scanpy as sc
from utils.params import Params
import matplotlib
import seaborn as sns
sns.set_style('ticks')
matplotlib.rcParams['figure.figsize'] = [8, 8]
matplotlib.rcParams['figure.dpi'] = 300
params = Params()


def get_meta_ad(adata_i, n_SEACells):
    model = SEACells.core.SEACells(adata_i,
                                   build_kernel_on='X_pca',
                                   n_SEACells=n_SEACells,
                                   n_waypoint_eigs=10,
                                   convergence_epsilon=1e-5
                                   )
    model.construct_kernel_matrix()
    model.initialize_archetypes()
    model.fit(min_iter=10, max_iter=500)
    SEACell_ad = SEACells.core.summarize_by_SEACell(adata_i, SEACells_label='SEACell', summarize_layer='raw')
    for c in ['individual', 'pool', 'day', 'treatment', 'GEM', 'sample_id', "sample", 'status', "gender"]:
        SEACell_ad.obs[c] = adata_i.obs[c].iloc[0]
    SEACell_ad.obs["old_obs_names"] = ""
    for name in SEACell_ad.obs_names:
        old_names = adata_i.obs[adata_i.obs["SEACell"] == name].index.tolist()
        SEACell_ad.obs.loc[SEACell_ad.obs_names == name, 'old_obs_names'] = "\t".join(old_names)
    SEACell_ad.obs_names = prefix + SEACell_ad.obs_names
    return SEACell_ad


if __name__ == '__main__':
    SEACell_ads = []
    adata = ad.read_h5ad(params.folder + "als.h5ad")
    # sc.pp.normalize_total(adata, target_sum=1e4)
    # sc.pp.log1p(adata)
    print(len(adata))
    unique_pools = np.unique(adata.obs['pool'])
    for pool in unique_pools:
        adata_p = adata[adata.obs["pool"] == pool].copy()
        unique_samples = np.unique(adata_p.obs['sample_id'])
        for sample in unique_samples:
            adata_s = adata_p[adata_p.obs["sample_id"] == sample].copy()
            unique_individuals = np.unique(adata_s.obs['individual'])
            sc.tl.pca(adata_s)
            for individual in unique_individuals:
                adata_i = adata_s[adata_s.obs["individual"] == individual].copy()
                adata_i.obs["old_obs_names"] = adata_i.obs_names

                n_SEACells = int(math.sqrt(len(adata_i)))
                if n_SEACells <= 1:
                    SEACell_ads.append(adata_i.raw.to_adata())
                    continue

                prefix = f"{pool}_{sample}_{individual}_"
                print("=======================================================================")
                print("=======================================================================")
                print("=======================================================================")
                print(prefix)
                print(n_SEACells)
                SEACell_ad = None

                for i in range(n_SEACells + 1):
                    print(f"Trying n_SEACells {n_SEACells + i}")
                    try:
                        SEACell_ad = get_meta_ad(adata_i, n_SEACells + i)
                        break
                    except:
                        pass
                    if i != 0 and n_SEACells - i > 1:
                        print(f"Trying n_SEACells {n_SEACells - i}")
                        try:
                            SEACell_ad = get_meta_ad(adata_i, n_SEACells - i)
                            break
                        except:
                            pass
                if SEACell_ad is not None:
                    SEACell_ads.append(SEACell_ad)
                    with open("passed.txt", "a+") as file:
                        file.write(prefix + "\t" + str(len(adata_i)) + "\t" + str(int(math.sqrt(len(adata_i)))) + "\t" + str(len(SEACell_ad)) + "\n")
                else:
                    SEACell_ads.append(adata_i.raw.to_adata())
                    with open("failed.tsv", "a+") as file:
                        file.write(prefix + "\t" + str(len(adata_i)) + "\t" + str(int(math.sqrt(len(adata_i)))) + "\n")

    meta_ad = ad.concat(SEACell_ads)
    print(len(meta_ad))
    meta_ad.write(params.folder + "als_meta.h5ad")
