import scanpy as sc
import anndata as ad
import numpy as np
import SEACells
from sklearn.utils import shuffle
import scanpy.external as sce
import utils.common as cm
from sklearn.ensemble import RandomForestClassifier
from collections import Counter
from params import Params
import pandas as pd
from random import sample
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.metrics import make_scorer, f1_score
from sklearn.model_selection import cross_val_score
import matplotlib
matplotlib.use('Agg')


cluster_col = "manual_anno_L0"
params = Params()
adata = ad.read_h5ad(params.folder + "ad_files/als_meta.h5ad")
adata = adata[adata.obs["gender"] == "Female"].copy()
# adata = adata[adata.obs["treatment"] != "Ropi"].copy()

# healthy_individuals = adata[(adata.obs["status"] == "healthy")
#                             & (adata.obs["gender"] == "Male")].obs['individual'].unique().tolist()
# als_individuals = adata[(adata.obs["status"] == "ALS")
#                         & (adata.obs["gender"] == "Male")].obs['individual'].unique().tolist()[len(healthy_individuals):]
# adata = adata[~adata.obs["individual"].isin(als_individuals)].copy()

n_pcs = 10
individuals = adata.obs["individual"].unique().tolist()
individual_to_status = adata.obs.groupby('individual')["status"].first().to_dict()
data = {}
print(adata.obs[cluster_col].unique())
PC_genes = {}
for cluster in ['iPSC', 'pMN', 'NSC', 'oMN']:
    adata_cluster = adata[adata.obs[cluster_col] == cluster].copy()
    SEACell_ad = SEACells.core.summarize_by_SEACell(adata_cluster, SEACells_label='individual', summarize_layer='raw')
    # for name in SEACell_ad.obs_names:
    #     SEACell_ad.obs["gender"] = adata_cluster[adata_cluster.obs["individual"] == name].obs["gender"].tolist()[0]
    # cm.highly_var(SEACell_ad)
    # sc.pp.normalize_total(SEACell_ad, target_sum=1e4)
    sc.pp.log1p(SEACell_ad)
    sc.tl.pca(SEACell_ad, n_comps=n_pcs)

    # sce.pp.harmony_integrate(SEACell_ad, 'gender', max_iter_harmony=50, max_iter_kmeans=100, block_size=0.01, sigma=0.2,
    #                          epsilon_cluster=float('-inf'), epsilon_harmony=float('-inf'))

    feature_loadings = SEACell_ad.varm['PCs']
    for pc_index in range(feature_loadings.shape[1]):
        # Get the genes contributing to the current PC and their corresponding loadings
        gene_indices = feature_loadings[:, pc_index].argsort()[::-1]  # Sort in descending order
        genes = adata_cluster.var_names[gene_indices]
        PC_genes[cluster + "_" + str(pc_index)] = genes[:500]

    for individual in individuals:
        if individual in SEACell_ad.obs_names:
            data.setdefault(individual, []).extend(SEACell_ad[SEACell_ad.obs_names == individual].obsm['X_pca'].flatten())
        else:
            data.setdefault(individual, []).extend(np.zeros(n_pcs))


pd_data = pd.DataFrame.from_dict(data, orient='index', columns=list(PC_genes.keys()))
pd_data["status"] = pd_data.index.map(individual_to_status)
pd_data['status'] = pd_data['status'].replace({'ALS': 0, 'healthy': 1})
pd_data = pd_data.sample(frac=1) # , random_state=1


X = pd_data.drop('status', axis=1)
y = pd_data['status']
X, y = shuffle(X, y)
print(y.value_counts())
# Create a Logistic Regression model
classifier = LogisticRegression()

# Perform 3-fold cross-validation
cv_scores = cross_val_score(classifier, X, y, cv=3, scoring='accuracy')

# Print cross-validation scores
for fold, score in enumerate(cv_scores, start=1):
    print(f"Fold {fold}: Accuracy = {score:.2f}")

# Print mean and standard deviation of cross-validation scores
print(f"Mean Accuracy: {np.mean(cv_scores):.2f}")
print(f"Std Deviation: {np.std(cv_scores):.2f}")

classifier = LogisticRegression()
classifier.fit(X, y)

# Get the coefficients assigned to each feature
coefficients = classifier.coef_[0]

# Create a dictionary with feature names and their corresponding coefficients
feature_coeff_dict = dict(zip(X.columns, coefficients))

# Sort the features by their coefficients in descending order
sorted_features = sorted(feature_coeff_dict.items(), key=lambda x: -x[1], reverse=True) # abs(x[1])

all_genes = []
# Print the top N most important features
top_n = 10  # Change this number to the desired top features count
i = 1
for feature, coeff in sorted_features[:top_n]:
    print(f"Feature: {feature}, Coefficient: {coeff:.4f}, Genes {PC_genes[feature]}")
    all_genes.extend(PC_genes[feature][:1 + int(len(PC_genes[feature]) / i)]) # Less genes from un important features
    i += 1
print("\n\n\n")

# Use Counter to count occurrences of each element
element_counts = Counter(all_genes)
# Sort the elements by count in descending order
sorted_elements = sorted(element_counts.items(), key=lambda x: x[1], reverse=True)
# Print the counts for each element

for element, count in sorted_elements:
    # print(f"Element: {element}, Count: {count}")
    print(element)
