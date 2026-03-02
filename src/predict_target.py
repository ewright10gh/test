import numpy as np
import pandas as pd
from joblib import load
import os

# =========================
# PATHS
# =========================
X_TARGET_PATH = r"C:\Users\mba22ew\test\cleaned\X_target.npy"
MODEL_PATH  = r"C:\Users\mba22ew\test\models\ridge_eln_model.joblib"
SCALER_PATH = r"C:\Users\mba22ew\test\models\ridge_scaler.joblib"
OUTPUT_PATH = r"C:\Users\mba22ew\test\results\target_eln_predictions.csv"
EXPR_TARGET_PATH = r"C:\Users\mba22ew\test\cleaned\target_cleaned_expression.csv"  # DataFrame with gene names

# =========================
# LOAD DATA
# =========================
print("Loading TARGET expression matrix...")
X_target = np.load(X_TARGET_PATH)
expr_target = pd.read_csv(EXPR_TARGET_PATH, index_col=0)  # For gene-level analysis
genes = expr_target.columns
samples = expr_target.index

print("Loading Ridge model + scaler...")
model  = load(MODEL_PATH)
scaler = load(SCALER_PATH)

# =========================
# SCALE
# =========================
X_target_scaled = scaler.transform(X_target)

# =========================
# PREDICT
# =========================
print("Predicting ELN favourable status...")
probs = model.predict_proba(X_target_scaled)[:, 1]
preds = (probs >= 0.5).astype(int)

results = pd.DataFrame({
    "ELN_favourable_probability": probs,
    "ELN_predicted_class": preds
}, index=samples)

# =========================
# SAVE PREDICTIONS
# =========================
os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
results.to_csv(OUTPUT_PATH)
print("Predictions saved ✅")
print(results["ELN_predicted_class"].value_counts())

# =========================
# ANALYZE TOP PREDICTIVE GENES
# =========================
print("\nAnalyzing gene expression patterns for predicted classes...")

# Load top 20 genes from Ridge model
feature_importances_path = r"C:\Users\mba22ew\test\cleaned\ridge_feature_importance.csv"
top_genes = pd.read_csv(feature_importances_path, index_col=0).head(20).index.tolist()

# Include HOXA cluster explicitly if not in top 20
hoxa_genes = [g for g in genes if g.startswith("HOXA")]
genes_to_check = list(set(top_genes + hoxa_genes))

expr_subset = expr_target[genes_to_check]

# Split expression by predicted class
expr_fav = expr_subset.loc[results["ELN_predicted_class"] == 1]
expr_nonfav = expr_subset.loc[results["ELN_predicted_class"] == 0]

# Summary statistics
summary = pd.DataFrame({
    "Favorable_mean": expr_fav.mean(),
    "NonFavorable_mean": expr_nonfav.mean(),
    "Favorable_median": expr_fav.median(),
    "NonFavorable_median": expr_nonfav.median()
}).sort_values("Favorable_mean", ascending=False)

print("\nTop genes expression summary by predicted class:")
print(summary)

