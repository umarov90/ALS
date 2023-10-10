import scanpy as sc
from params import Params
import anndata as ad
import utils.common as cm
import matplotlib
matplotlib.use('macosx')

p = Params()
adata = ad.read_h5ad(p.folder + "ad_files/als_filtered.h5ad")
# sc.pp.scale(adata, max_value=10)
cell_cycle_genes = [x.strip() for x in open(p.folder + 'regev_lab_cell_cycle_genes.txt')]
print(len(cell_cycle_genes))

# Split into 2 lists
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]

cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
print(len(cell_cycle_genes))

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
sc.pl.violin(adata, ['S_score', 'G2M_score'], jitter=0.4, groupby='sample_id', rotation=45)
cm.safe_save(adata, p.folder + "als_filtered.h5ad")
