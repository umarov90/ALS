import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import anndata as ad
from datetime import datetime
from params import Params
import numpy as np

p = Params()
# adata1 = ad.read_h5ad(p.folder + p.file_name.replace(".h5ad", "_vis.h5ad"))
adata = ad.read_h5ad(p.folder + "als_filtered.h5ad")

# Compute total counts per cell
adata.obs['total_counts'] = adata.X.sum(axis=1)

# Sort cells based on total counts
sorted_indices = np.argsort(adata.obs['total_counts'])

# Create a scatter plot
fig, ax = plt.subplots()
sc.pl.scatter(adata, x=sorted_indices, y='total_counts', ax=ax)

# Customize the axis labels
ax.set_xlabel('Ranked Cells')
ax.set_ylabel('Total Counts')

# Show the plot
plt.show()
