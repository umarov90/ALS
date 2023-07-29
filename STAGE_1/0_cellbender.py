import scrublet as scr
import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt
from params import Params
import anndata as ad

p = Params()
adata = ad.read_h5ad(p.file_path)
sce.pp.scrublet(adata)
adata.write(p.file_path)