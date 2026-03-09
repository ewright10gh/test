
# PREDICT ELN FAVOURABLE STATUS IN TARGET


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from joblib import load


# PATHS


BASE = r"C:\Users\mba22ew\test"

X_TARGET_PATH = os.path.join(BASE, "cleaned", "X_target.npy")
TARGET_EXPR_PATH = os.path.join(BASE, "cleaned", "target_cleaned_expression.csv")

MODEL_PATH = os.path.join(BASE, "models", "ridge_eln_model.joblib")

OUTDIR = os.path.join(BASE, "results")
os.makedirs(OUTDIR, exist_ok=True)


# LOAD DATA


print("Loading TARGET matrix...")
X_target = np.load(X_TARGET_PATH)

print("Loading TARGET expression for sample IDs...")
target_expr = pd.read_csv(TARGET_EXPR_PATH, index_col=0)
sample_ids = target_expr.index

print("Loading trained model...")
model = load(MODEL_PATH)

print("TARGET shape:", X_target.shape)


# PREDICT


print("Predicting ELN favourable status...")

prob = model.predict_proba(X_target)[:, 1]
pred = (prob >= 0.5).astype(int)

results = pd.DataFrame({
    "ELN_favourable_probability": prob,
    "ELN_predicted_class": pred
}, index=sample_ids)


# HIGH-CONFIDENCE CALLS


high_conf_fav = results[results["ELN_favourable_probability"] >= 0.7]
high_conf_adv = results[results["ELN_favourable_probability"] <= 0.3]

print("\nPrediction counts:")
print(results["ELN_predicted_class"].value_counts())

print("\nHigh-confidence favourable:", len(high_conf_fav))
print("High-confidence adverse:", len(high_conf_adv))


# PROBABILITY DIAGNOSTICS


print("\nProbability summary:")
print(results["ELN_favourable_probability"].describe())


# SAVE TABLES


results.to_csv(os.path.join(OUTDIR, "TARGET_ELN_predictions.csv"))
high_conf_fav.to_csv(os.path.join(OUTDIR, "TARGET_high_conf_favourable.csv"))
high_conf_adv.to_csv(os.path.join(OUTDIR, "TARGET_high_conf_adverse.csv"))


# PLOTS (FOR PAPER)


plt.figure()
results["ELN_favourable_probability"].hist(bins=50)
plt.title("TARGET ELN favourable probability distribution")
plt.xlabel("Probability")
plt.ylabel("Number of patients")
plt.savefig(os.path.join(OUTDIR, "probability_histogram.png"))
plt.close()

plt.figure()
results["ELN_favourable_probability"].sort_values().reset_index(drop=True).plot()
plt.title("TARGET patients ranked by ELN favourable probability")
plt.xlabel("Patients")
plt.ylabel("Probability")
plt.savefig(os.path.join(OUTDIR, "probability_ranked_curve.png"))
plt.close()

print("\nAll outputs saved to:", OUTDIR)
print("Done ✅")


