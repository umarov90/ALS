import anndata as ad
import joblib
from params import Params
import utils.common as cm


p = Params()
adata = ad.read_h5ad(p.folder + "als_combined.h5ad")
ldata = adata.layers
joblib.dump(ldata, p.folder + "layers.p")
adata.layers.clear()
cm.safe_save(adata, p.file_path)