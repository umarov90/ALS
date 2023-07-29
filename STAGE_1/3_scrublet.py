import scanpy.external as sce
from params import Params
import anndata as ad
import utils.common as cm
import matplotlib
matplotlib.use('macosx')
import gc

p = Params()
adata = ad.read_h5ad(p.file_path)

adatas = []
unique_samples = list(adata.obs['sample_id'].unique())
for sample_id in unique_samples:
    adata_s = adata[adata.obs["sample_id"] == sample_id].copy()
    a = len(adata_s)
    # if sample_id in ["TFHS000382"]:
    dbl = adata_s.obs['droplet_type'].value_counts()['DBL']
    amb = adata_s.obs['droplet_type'].value_counts()['AMB']
    dbl_fraction = dbl / len(adata_s)
    amb_fraction = amb / len(adata_s)
    print(f"{sample_id} {dbl} {amb} {dbl_fraction} {amb_fraction}")
    sce.pp.scrublet(adata_s, expected_doublet_rate=dbl_fraction, stdev_doublet_rate=0.02) # , threshold=0.38

    # plt.hist(adata.obs['doublet_score'], bins=20, range=(0, 1.0))
    # plt.xlabel("Values")
    # plt.ylabel("Frequency")
    # plt.title("Histogram")
    # plt.show()

    del adata_s.uns['scrublet']
    adatas.append(adata_s)

del adata
gc.collect()

adata = ad.concat(adatas)
cm.safe_save(adata, p.file_path)
