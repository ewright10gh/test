# =========================
# TARGET ELN PREDICTION + ANALYSIS
# =========================

import os
import numpy as np
import pandas as pd
import joblib
import matplotlib.pyplot as plt
import seaborn as sns

# -------------------------
# PATHS
# -------------------------

BASE = r"C:\Users\mba22ew\test"

X_TARGET_PATH = os.path.join(BASE, "cleaned", "X_target.npy")
TARGET_EXPR_PATH = os.path.join(BASE, "cleaned", "target_cleaned_expression.csv")

MODEL_PATH = os.path.join(BASE, "models", "ridge_eln_model.joblib")
SCALER_PATH = os.path.join(BASE, "models", "ridge_scaler.joblib")

OUT_DIR = os.path.join(BASE, "results")
os.makedirs(OUT_DIR, exist_ok=True)

# -------------------------
# LOAD DATA
# -------------------------

print("Loading TARGET matrix...")
X_target = np.load(X_TARGET_PATH)

print("Loading TARGET expression for patient IDs...")
target_expr = pd.read_csv(TARGET_EXPR_PATH, index_col=0)

print("Loading model + scaler...")
model = joblib.load(MODEL_PATH)
scaler = joblib.load(SCALER_PATH)

# -------------------------
# SCALE + PREDICT
# -------------------------

X_target_scaled = scaler.transform(X_target)

print("Predicting ELN favourable status...")
prob = model.predict_proba(X_target_scaled)[:, 1]
pred = (prob >= 0.5).astype(int)

# -------------------------
# BUILD RESULTS TABLE
# -------------------------

results = pd.DataFrame({
    "ELN_favourable_probability": prob,
    "ELN_predicted_class": pred
}, index=target_expr.index)

results.to_csv(os.path.join(OUT_DIR, "target_eln_predictions.csv"))

print("\nPrediction counts:")
print(results["ELN_predicted_class"].value_counts())

# -------------------------
# CONFIDENCE STRATIFICATION
# -------------------------
high_conf_fav = results[results.ELN_favourable_probability >= 0.6]
high_conf_adv = results[results.ELN_favourable_probability <= 0.2]

print("High-confidence favourable:", len(high_conf_fav))
print("High-confidence adverse:", len(high_conf_adv))

high_conf_fav.to_csv(os.path.join(OUT_DIR, "target_high_conf_favourable.csv"))
high_conf_adv.to_csv(os.path.join(OUT_DIR, "target_high_conf_adverse.csv"))

