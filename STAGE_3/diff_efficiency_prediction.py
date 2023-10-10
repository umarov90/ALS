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
import xgboost as xgb
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('macosx')

cluster_col = "manual_anno_L0"
params = Params()
adata = ad.read_h5ad(params.folder + "ad_files/als_filtered_anno.h5ad")
individuals = adata.obs["individual"].unique().tolist()
X = []
y = []
adata_cluster_mMN = adata[(adata.obs[cluster_col] == "mMN") & (adata.obs["day"] == "D7")].copy()
adata_cluster_iPSC = adata[adata.obs[cluster_col] == "iPSC"].copy()
for individual in individuals:
    initial_num = len(adata_cluster_iPSC[adata_cluster_iPSC.obs['individual'] == individual])
    val = (adata_cluster_mMN.obs['individual'] == individual).sum() / initial_num
    try:
        X.append(adata_cluster_iPSC[adata_cluster_iPSC.obs['individual'] == individual].X.mean(axis=0))
        y.append(val < 2)
    except:
        print(f"Individual {individual} has no cells in iPSC.")

X = np.squeeze(np.asarray(X))

# plt.hist(y, bins=20)
# plt.show()

print(f"Positives {y.count(True)}")
print(f"Negatives {y.count(False)}")
# Define the XGBoost classifier
xgb_classifier = xgb.XGBClassifier(
    objective="binary:logistic",
    n_estimators=100,
    max_depth=3,
    learning_rate=0.1,
    random_state=1
)

# Define the StratifiedKFold cross-validator with 5 folds
cv = StratifiedKFold(n_splits=5, shuffle=True)

# Perform cross-validation and calculate AUC for each fold
auc_scores = cross_val_score(xgb_classifier, X, y, cv=cv, scoring="roc_auc")

# Print the AUC scores for each fold
for fold, auc in enumerate(auc_scores, start=1):
    print(f"Fold {fold}: AUC = {auc:.4f}")

# Calculate and print the mean AUC across all folds
mean_auc = np.mean(auc_scores)
print(f"Mean AUC: {mean_auc:.4f}")

# Fit the model to your data
xgb_classifier.fit(X, y)

# Get feature importances
feature_importances = xgb_classifier.feature_importances_

# Get the names of your features
feature_names = adata.var_names.tolist()

# Create a sorted list of (feature, importance) pairs
feature_importance_pairs = [(feature, importance) for feature, importance in zip(feature_names, feature_importances) if importance >= 0.02]
sorted_feature_importance_pairs = sorted(feature_importance_pairs, key=lambda x: x[1], reverse=True)

# Print or visualize the feature importances
for feature, importance in sorted_feature_importance_pairs:
    print(f"{feature}: {importance:.4f}")

# Plot feature importances
plt.figure(figsize=(5, 10))
plt.barh([pair[0] for pair in sorted_feature_importance_pairs], [pair[1] for pair in sorted_feature_importance_pairs])
plt.xlabel('Importance')
plt.title('Gene Importance')
plt.gca().invert_yaxis()  # Invert the y-axis for better visualization
plt.savefig(params.folder + "gene_importance.pdf")