
# PREDICT FAVOURABLE FUSION STATUS IN TARGET


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from joblib import load


# PATHS

script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)

X_TARGET_PATH = os.path.join(parent_dir, "cleaned", "X_target.npy")
TARGET_EXPR_PATH = os.path.join(parent_dir, "cleaned", "target_cleaned_expression.csv")

MODEL_PATH = os.path.join(parent_dir, "models", "ridge_eln_model.joblib")

OUTDIR = os.path.join(parent_dir, "results")
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


print("Predicting favourable fusion probability...")

prob = model.predict_proba(X_target)[:, 1]
pred = (prob >= 0.5).astype(int)

results = pd.DataFrame({
    "favourable_fusion_probability": prob,
    "favourable_fusion_predicted_class": pred
}, index=sample_ids)


# HIGH-CONFIDENCE CALLS


high_conf_fav = results[results["favourable_fusion_probability"] >= 0.675]
high_conf_nonfav = results[results["favourable_fusion_probability"] <= 0.3]

print("\nPrediction counts:")
print(results["favourable_fusion_predicted_class"].value_counts())

print("\nHigh-confidence favourable:", len(high_conf_fav))
print("High-confidence non-favourable:", len(high_conf_nonfav))


# PROBABILITY DIAGNOSTICS


print("\nProbability summary:")
print(results["favourable_fusion_probability"].describe())


# SAVE TABLES


results.to_csv(os.path.join(OUTDIR, "TARGET_favourable_fusion_predictions.csv"))
high_conf_fav.to_csv(os.path.join(OUTDIR, "TARGET_high_conf_favourable.csv"))
high_conf_nonfav.to_csv(os.path.join(OUTDIR, "TARGET_high_conf_nonfavourable.csv"))


# PLOTS (FOR PAPER)


plt.figure()
results["favourable_fusion_probability"].hist(bins=50)
plt.title("TARGET favourable fusion probability distribution")
plt.xlabel("Probability")
plt.ylabel("Number of patients")
plt.savefig(os.path.join(OUTDIR, "probability_histogram.png"))
plt.close()

plt.figure()
results["favourable_fusion_probability"].sort_values().reset_index(drop=True).plot()
plt.title("TARGET patients ranked by favourable fusion probability")
plt.xlabel("Patients")
plt.ylabel("Probability")
plt.savefig(os.path.join(OUTDIR, "probability_ranked_curve.png"))
plt.close()

print("\nAll outputs saved to:", OUTDIR)
print("Done ✅")


