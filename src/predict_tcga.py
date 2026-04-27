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


# LOAD TCGA EXPRESSION

X = pd.read_csv(
    os.path.join(parent_dir, "cleaned", "tcga_cleaned_expression.csv"),
    index_col=0
)

# remove metadata row if present
X = X[~X.index.str.contains("Entrez", case=False, na=False)]

print("Expression matrix shape:", X.shape)


# SCALE USING TRAINING SCALER

X_scaled = scaler.transform(X.values)


# PREDICT

probs = model.predict_proba(X_scaled)[:,1]

results = pd.DataFrame({
    "Sample_ID": X.index,
    "favourable_fusion_probability": probs
})

results["favourable_fusion_predicted"] = (results["favourable_fusion_probability"] > 0.675).astype(int)

print("\nPrediction summary:")
print(results["favourable_fusion_predicted"].value_counts())


# SAVE RESULTS

results.to_csv(os.path.join(parent_dir, "results", "tcga_favourable_fusion_predictions.csv"), index=False)

print("\n✅ TCGA predictions saved")