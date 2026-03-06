import pandas as pd
import numpy as np
import joblib

# -----------------------------
# LOAD MODEL
# -----------------------------
model = joblib.load(r"C:\Users\mba22ew\test\models\ridge_eln_model.joblib")
scaler = joblib.load(r"C:\Users\mba22ew\test\models\ridge_scaler.joblib")

print("Model + scaler loaded")

# -----------------------------
# LOAD TCGA EXPRESSION
# -----------------------------
X = pd.read_csv(
    r"C:\Users\mba22ew\test\cleaned\tcga_cleaned_expression.csv",
    index_col=0
)

# remove metadata row if present
X = X[~X.index.str.contains("Entrez", case=False, na=False)]

print("Expression matrix shape:", X.shape)

# -----------------------------
# SCALE USING TRAINING SCALER
# -----------------------------
X_scaled = scaler.transform(X.values)

# -----------------------------
# PREDICT
# -----------------------------
probs = model.predict_proba(X_scaled)[:,1]

results = pd.DataFrame({
    "Sample_ID": X.index,
    "ELN_probability": probs
})

results["ELN_predicted"] = (results["ELN_probability"] > 0.7).astype(int)

print("\nPrediction summary:")
print(results["ELN_predicted"].value_counts())

# -----------------------------
# SAVE RESULTS
# -----------------------------
results.to_csv(r"C:\Users\mba22ew\test\results\tcga_eln_predictions.csv", index=False)

print("\n✅ TCGA predictions saved")