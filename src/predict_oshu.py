import pandas as pd
import numpy as np
import joblib
import os

# Get the directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)

# LOAD MODEL
model = joblib.load(os.path.join(parent_dir, "models", "ridge_eln_model.joblib"))
scaler = joblib.load(os.path.join(parent_dir, "models", "ridge_scaler.joblib"))

print("Model + scaler loaded")

# LOAD OHSU VALIDATION EXPRESSION
X = pd.read_csv(
    os.path.join(parent_dir, "cleaned", "validata_cleaned_expression.csv"),
    index_col=0
)

print("Expression matrix shape:", X.shape)

print("Expression matrix shape:", X.shape)

# SCALE USING TRAINING SCALER
X_scaled = scaler.transform(X.values)

# HANDLE NaN ENTRIES (rare due to missing/dropped genes)
nan_rows = np.isnan(X_scaled).any(axis=1)
if nan_rows.any():
    print(f"Warning: dropping {nan_rows.sum()} samples with NaN after scaling")
    X = X.loc[~nan_rows]
    X_scaled = X_scaled[~nan_rows]

# PREDICT
probs = model.predict_proba(X_scaled)[:, 1]

results = pd.DataFrame({
    "Sample_ID": X.index,
    "favourable_fusion_probability": probs
})

results["favourable_fusion_predicted"] = (results["favourable_fusion_probability"] > 0.675).astype(int)

print("\nPrediction summary:")
print(results["favourable_fusion_predicted"].value_counts())

# SAVE RESULTS
results.to_csv(os.path.join(parent_dir, "results", "validata_favourable_fusion_predictions.csv"), index=False)

print("\n✅ VALIDATA predictions saved")
