import numpy as np
import pandas as pd

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

X = pd.read_csv("data/processed/ml_er_features.csv", index_col=0)
meta = pd.read_csv("data/processed/brca_metadata_er.csv")
y = meta["ER_status"].map({"ER_neg": 0, "ER_pos": 1}).to_numpy()

X_log = np.log1p(X.values)

pipe = Pipeline([
    ("scaler", StandardScaler()),
    ("clf", LogisticRegression(max_iter=5000, class_weight="balanced"))
])

pipe.fit(X_log, y)

coefs = pipe.named_steps["clf"].coef_.ravel()
genes = X.columns

df = pd.DataFrame({"gene": genes, "coef": coefs})
df["abs_coef"] = df["coef"].abs()
df = df.sort_values("abs_coef", ascending=False)

df.head(25).to_csv("results/ml_er_top25_logreg_genes.csv", index=False)
df.to_csv("results/ml_er_logreg_all_coeffs.csv", index=False)

print("Saved: results/ml_er_top25_logreg_genes.csv")
print(df.head(15))
