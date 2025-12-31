import pandas as pd
import numpy as np

counts = pd.read_csv("data/processed/brca_counts.csv", index_col=0)
meta   = pd.read_csv("data/processed/brca_metadata.csv")

# Extract ESR1 expression
esr1 = counts.loc["ESR1"]
esr1_log = np.log2(esr1 + 1)

# Quantiles
q_low  = esr1_log.quantile(0.40)
q_high = esr1_log.quantile(0.60)

er_status = pd.Series(index=esr1_log.index, dtype="object")
er_status[esr1_log <= q_low]  = "ER_neg"
er_status[esr1_log >= q_high] = "ER_pos"

meta_er = meta.copy()
meta_er["ER_status"] = er_status.values

# Drop ambiguous samples
meta_er = meta_er.dropna(subset=["ER_status"])

meta_er.to_csv("data/processed/brca_metadata_er.csv", index=False)

print(meta_er["ER_status"].value_counts())
print("Saved: data/processed/brca_metadata_er.csv")
