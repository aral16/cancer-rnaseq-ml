import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve

X = pd.read_csv("data/processed/ml_top200_features.csv", index_col=0)
meta = pd.read_csv("data/processed/brca_metadata.csv")

# label encoding: normal=0, tumor=1
y = meta["condition"].map({"normal": 0, "tumor": 1}).to_numpy()
if np.any(pd.isna(y)):
    raise ValueError("Found condition values not in {normal, tumor}. Check brca_metadata.csv")

# IMPORTANT: counts are skewed; log transform helps stability
X_log = np.log1p(X.values)

pipe = Pipeline([
    ("scaler", StandardScaler()),
    ("clf", LogisticRegression(max_iter=5000, class_weight="balanced"))
])

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# Out-of-fold predicted probabilities (no leakage)
proba = cross_val_predict(pipe, X_log, y, cv=cv, method="predict_proba")[:, 1]

auc = roc_auc_score(y, proba)
print(f"5-fold CV ROC-AUC (out-of-fold): {auc:.4f}")

# ROC curve plot
fpr, tpr, _ = roc_curve(y, proba)
plt.figure()
plt.plot(fpr, tpr)
plt.plot([0, 1], [0, 1], linestyle="--")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title(f"BRCA Tumor vs Normal ROC (5-fold OOF) | AUC={auc:.3f}")
plt.tight_layout()
plt.savefig("results/ml_roc_curve.png", dpi=200)
print("Saved: results/ml_roc_curve.png")
