import os
import pandas as pd
import numpy as np

BASE = r"C:\Users\mba22ew\test"

OHSU_EXPR = os.path.join(BASE, "cleaned", "ohsu_cleaned_expression.csv")
LABELS = os.path.join(BASE, "cleaned", "ohsu_favorable_labels.csv")

OUTDIR = os.path.join(BASE, "results")
os.makedirs(OUTDIR, exist_ok=True)

print("Loading OHSU expression...")
expr = pd.read_csv(OHSU_EXPR, index_col=0)

print("Loading ELN labels...")
labels = pd.read_csv(LABELS, index_col=0).iloc[:, 0]

# Ensure same order
expr = expr.loc[labels.index]

print("Samples:", expr.shape[0])
print("Genes:", expr.shape[1])

# Split groups
fav = expr[labels == 1]
nonfav = expr[labels == 0]

print("Favourable:", fav.shape[0])
print("Non-favourable:", nonfav.shape[0])

# Compute statistics
stats = pd.DataFrame({
    "Favorable_mean": fav.mean(),
    "NonFavorable_mean": nonfav.mean(),
    "Favorable_median": fav.median(),
    "NonFavorable_median": nonfav.median(),
})

stats["mean_diff"] = stats["Favorable_mean"] - stats["NonFavorable_mean"]

# Sort
top_fav = stats.sort_values("mean_diff", ascending=False).head(25)
top_adv = stats.sort_values("mean_diff", ascending=True).head(25)

# Save
top_fav.to_csv(os.path.join(OUTDIR, "OHSU_top_favourable_expression_genes.csv"))
top_adv.to_csv(os.path.join(OUTDIR, "OHSU_top_adverse_expression_genes.csv"))

print("\nTop favourable genes:")
print(top_fav.head(20))

print("\nTop adverse genes:")
print(top_adv.head(20))

print("\n✅ OHSU expression-based top genes saved")