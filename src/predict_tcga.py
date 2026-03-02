import numpy as np
import pandas as pd
from joblib import load

# load data
X_tcga = np.load(r"C:\Users\mba22ew\test\cleaned\X_tcga.npy")

# load model + scaler
model = load(r"C:\Users\mba22ew\test\models\model.joblib")
scaler = load(r"C:\Users\mba22ew\test\models\scaler.joblib")

# scale
X_tcga = scaler.transform(X_tcga)

# predict
probs = model.predict_proba(X_tcga)[:, 1]

# save
pd.DataFrame({
    "Sample": range(len(probs)),
    "Favourable_Probability": probs
}).to_csv(
    r"C:\Users\mba22ew\test\results\tcga_predictions.csv",
    index=False
)

print("TCGA prediction complete")