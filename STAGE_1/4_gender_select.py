import anndata as ad
import numpy as np
from params import Params
import utils.common as cm
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('macosx')

p = Params()
adata = ad.read_h5ad(p.file_path)

unique_individuals = np.unique(adata.obs['individual'])
normalized_expression = {}
for ind in unique_individuals:
    print(ind)
    indices = adata.obs['individual'] == ind
    gene_expression = np.sum(adata[indices, 'RPS4Y1'].X)
    total_expression = np.sum(adata[indices, :].X)
    normalized_expression[ind] = gene_expression / total_expression

plt.hist(normalized_expression.values(), bins=20, range=(0, 0.0005))
plt.xlabel("Values")
plt.ylabel("Frequency")
plt.title("Histogram")
plt.show()

gender_mapping = {"Male": [], "Female": []}
for ind, expr in normalized_expression.items():
    if expr > 0.0001:
        gender_mapping["Male"].append(ind)
    else:
        gender_mapping["Female"].append(ind)

# sorted_expression_sum = sorted(normalized_expression.items(), key=lambda x: x[1])
# split_point = len(sorted_expression_sum) // 2
# gender_mapping = {"Female": [item[0] for item in sorted_expression_sum[:split_point]],
#                   "Male": [item[0] for item in sorted_expression_sum[split_point:]]}
print(gender_mapping)
adata.obs["gender"] = adata.obs["individual"].map(lambda x: "Male" if x in gender_mapping["Male"] else "Female")
cm.safe_save(adata, p.file_path)
