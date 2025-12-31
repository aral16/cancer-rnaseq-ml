import numpy as np
import pandas as pd

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

X = pd.read_csv("data/processed/ml_top200_features.csv", index_col=0)
meta = pd.read_csv("data/processed/brca_metadata.csv")
y = meta["condition"].map({"normal": 0, "tumor": 1}).to_numpy()

X_log = np.log1p(X.values)

pipe = Pipeline([
    ("scaler", StandardScaler()),
    ("clf", LogisticRegression(max_iter=5000, class_weight="balanced"))
])

pipe.fit(X_log, y)

# Coefficients correspond to standardized features (still useful for ranking)
coefs = pipe.named_steps["clf"].coef_.ravel()
genes = X.columns

df = pd.DataFrame({"gene": genes, "coef": coefs})
df["abs_coef"] = df["coef"].abs()
df = df.sort_values("abs_coef", ascending=False)

df.head(25).to_csv("results/ml_top25_logreg_genes.csv", index=False)
df.to_csv("results/ml_logreg_all_coeffs.csv", index=False)

print("Saved: results/ml_top25_logreg_genes.csv")
print("Saved: results/ml_logreg_all_coeffs.csv")
print(df.head(10))
