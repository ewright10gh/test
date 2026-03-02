import os
import pandas as pd

BASE = r"C:\Users\mba22ew\test"

PRED_PATH = os.path.join(BASE, "results", "TARGET_ELN_predictions.csv")
EXPR_PATH = os.path.join(BASE, "cleaned", "target_cleaned_expression.csv")
OUT_PATH = os.path.join(BASE, "results", "TARGET_top_genes_expression_table.csv")

print("Loading data...")

pred = pd.read_csv(PRED_PATH, index_col=0)
expr = pd.read_csv(EXPR_PATH, index_col=0)

expr = expr.loc[pred.index]

fav = expr[pred["ELN_predicted_class"] == 1]
adv = expr[pred["ELN_predicted_class"] == 0]

print("Fav:", fav.shape)
print("Adv:", adv.shape)

# =============================
# CALCULATE STATS
# =============================

table = pd.DataFrame(index=expr.columns)

table["Favorable_mean"] = fav.mean()
table["NonFavorable_mean"] = adv.mean()

table["Favorable_median"] = fav.median()
table["NonFavorable_median"] = adv.median()

# Difference for ranking
table["mean_diff"] = table["Favorable_mean"] - table["NonFavorable_mean"]

# Sort by difference
table = table.sort_values("mean_diff", ascending=False)

# Save full table
table.to_csv(OUT_PATH)

print("\nTop favourable genes:")
print(table.head(20))

print("\nTop adverse genes:")
print(table.tail(20))