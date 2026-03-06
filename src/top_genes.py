# ============================================
# COMPARE TOP GENES ACROSS OHSU, TCGA, TARGET
# ============================================

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from joblib import load

BASE = r"C:\Users\mba22ew\test"

OHSU_EXPR = os.path.join(BASE, "cleaned", "ohsu_cleaned_expression.csv")
TCGA_EXPR = os.path.join(BASE, "cleaned", "tcga_cleaned_expression.csv")
TARGET_EXPR = os.path.join(BASE, "cleaned", "target_cleaned_expression.csv")

OHSU_LABELS = os.path.join(BASE, "cleaned", "ohsu_favorable_labels.csv")
TCGA_PRED = os.path.join(BASE, "results", "TCGA_ELN_predictions.csv")
TARGET_PRED = os.path.join(BASE, "results", "TARGET_ELN_predictions.csv")

MODEL_PATH = os.path.join(BASE, "models", "ridge_eln_model.joblib")

OUTDIR = os.path.join(BASE, "results")
os.makedirs(OUTDIR, exist_ok=True)

# ============================================
# LOAD DATA
# ============================================

print("Loading expression matrices...")

ohsu = pd.read_csv(OHSU_EXPR, index_col=0)
tcga = pd.read_csv(TCGA_EXPR, index_col=0)
target = pd.read_csv(TARGET_EXPR, index_col=0)

print("OHSU:", ohsu.shape)
print("TCGA:", tcga.shape)
print("TARGET:", target.shape)

# ============================================
# LOAD LABELS / PREDICTIONS
# ============================================

ohsu_labels = pd.read_csv(OHSU_LABELS, index_col=0).iloc[:,0]
tcga_pred = pd.read_csv(TCGA_PRED, index_col=0)
target_pred = pd.read_csv(TARGET_PRED, index_col=0)

# Align sample order
ohsu = ohsu.loc[ohsu_labels.index]
tcga = tcga.loc[tcga_pred.index]
target = target.loc[target_pred.index]

# ============================================
# FUNCTION TO COMPUTE TOP GENES
# ============================================

def compute_stats(expr, labels):

    fav = expr[labels == 1]
    nonfav = expr[labels == 0]

    stats = pd.DataFrame({
        "Favorable_mean": fav.mean(),
        "NonFavorable_mean": nonfav.mean(),
        "Favorable_median": fav.median(),
        "NonFavorable_median": nonfav.median()
    })

    stats["mean_diff"] = stats["Favorable_mean"] - stats["NonFavorable_mean"]

    return stats


# ============================================
# OHSU
# ============================================

print("\nComputing OHSU stats...")

ohsu_stats = compute_stats(ohsu, ohsu_labels)

top_ohsu = ohsu_stats.sort_values("mean_diff", ascending=False).head(50)

top_ohsu.to_csv(os.path.join(OUTDIR,"OHSU_top_genes.csv"))

# ============================================
# TCGA
# ============================================

print("Computing TCGA stats...")

tcga_stats = compute_stats(tcga, tcga_pred["ELN_predicted"])

top_tcga = tcga_stats.sort_values("mean_diff", ascending=False).head(50)

top_tcga.to_csv(os.path.join(OUTDIR,"TCGA_top_genes.csv"))

# ============================================
# TARGET
# ============================================

print("Computing TARGET stats...")

target_stats = compute_stats(target, target_pred["ELN_predicted_class"])

top_target = target_stats.sort_values("mean_diff", ascending=False).head(50)

top_target.to_csv(os.path.join(OUTDIR,"TARGET_top_genes.csv"))

# ============================================
# GENE OVERLAP
# ============================================

set_ohsu = set(top_ohsu.index)
set_tcga = set(top_tcga.index)
set_target = set(top_target.index)

overlap_all = set_ohsu & set_tcga & set_target
overlap_ohsu_tcga = set_ohsu & set_tcga
overlap_ohsu_target = set_ohsu & set_target

print("\nGene overlap:")
print("OHSU ∩ TCGA:", len(overlap_ohsu_tcga))
print("OHSU ∩ TARGET:", len(overlap_ohsu_target))
print("All three:", len(overlap_all))
print("Shared genes:", overlap_all)

# Save overlaps
pd.Series(list(overlap_all)).to_csv(
    os.path.join(OUTDIR,"shared_genes_all_datasets.csv"),
    index=False
)

# ============================================
# CHECK MYH11 RANKING
# ============================================

def check_gene(gene, stats):

    if gene in stats.index:
        rank = stats["mean_diff"].rank(ascending=False)[gene]
        print(f"{gene} rank:", int(rank))
        print(stats.loc[gene])
    else:
        print(gene, "not found")


print("\nChecking MYH11:")
check_gene("MYH11", ohsu_stats)
check_gene("MYH11", tcga_stats)
check_gene("MYH11", target_stats)

# ============================================
# COMPARE WITH RIDGE MODEL COEFFICIENTS
# ============================================

print("\nComparing model coefficients...")

model = load(MODEL_PATH)

coefs = pd.Series(model.coef_[0], index=ohsu.columns)

coef_df = pd.DataFrame({
    "ridge_coef": coefs,
    "expression_diff": ohsu_stats["mean_diff"]
})

coef_df = coef_df.dropna()

# ============================================
# PLOT COEFFICIENT VS EXPRESSION DIFFERENCE
# ============================================

plt.figure(figsize=(6,6))

plt.scatter(
    coef_df["expression_diff"],
    coef_df["ridge_coef"],
    alpha=0.4
)

plt.xlabel("Expression difference (favourable - adverse)")
plt.ylabel("Ridge coefficient")

plt.title("Model weight vs biological expression difference")

plt.savefig(os.path.join(OUTDIR,"ridge_vs_expression.png"))

plt.close()

# ============================================
# SAVE FULL TABLE
# ============================================

coef_df.to_csv(os.path.join(OUTDIR,"model_gene_weights_vs_expression.csv"))

print("\n✅ Gene comparison complete")
print("Results saved to:", OUTDIR)

# ============================================
# SIMPLE INTERPRETABLE FIGURES
# ============================================

import seaborn as sns
import matplotlib.pyplot as plt

FIGDIR = os.path.join(BASE, "results", "figures")
os.makedirs(FIGDIR, exist_ok=True)

# --------------------------------------------
# 1. DATASET GENE HEATMAP
# --------------------------------------------

print("Generating gene heatmap...")

# Combine mean expression differences
heatmap_df = pd.DataFrame({
    "OHSU": ohsu_stats["mean_diff"],
    "TCGA": tcga_stats["mean_diff"],
    "TARGET": target_stats["mean_diff"]
})

# Use only genes shared across datasets
heatmap_df = heatmap_df.loc[heatmap_df.index.intersection(
    ohsu_stats.index
).intersection(tcga_stats.index).intersection(target_stats.index)]

# Take strongest genes overall
heatmap_df["abs_mean"] = heatmap_df.abs().mean(axis=1)
heatmap_df = heatmap_df.sort_values("abs_mean", ascending=False).head(18)
heatmap_df = heatmap_df.drop(columns="abs_mean")

plt.figure(figsize=(8,6))

sns.heatmap(
    heatmap_df,
    cmap="vlag",
    center=0,
    annot=True,
    fmt=".2f"
)

plt.title("Top AML gene expression differences across datasets")
plt.ylabel("Gene")
plt.xlabel("Dataset")

plt.tight_layout()

plt.savefig(os.path.join(FIGDIR, "gene_dataset_heatmap.png"), dpi=300)
plt.close()


# --------------------------------------------
# 2. EFFECT SIZE COMPARISON
# --------------------------------------------

print("Generating dataset comparison scatter plot...")

common = ohsu_stats.index.intersection(tcga_stats.index)

plt.figure(figsize=(6,6))

plt.scatter(
    ohsu_stats.loc[common, "mean_diff"],
    tcga_stats.loc[common, "mean_diff"],
    alpha=0.6
)

plt.axhline(0)
plt.axvline(0)

plt.xlabel("OHSU gene effect size")
plt.ylabel("TCGA gene effect size")

plt.title("Gene signature consistency between datasets")

plt.tight_layout()

plt.savefig(os.path.join(FIGDIR, "gene_effect_scatter.png"), dpi=300)
plt.close()


# --------------------------------------------
# 3. TOP GENES BARPLOT
# --------------------------------------------

print("Generating top gene barplot...")

top_genes = ohsu_stats.sort_values("mean_diff", ascending=False).head(15)

plt.figure(figsize=(8,5))

top_genes["mean_diff"].sort_values().plot.barh()

plt.xlabel("Expression difference (favourable - adverse)")
plt.ylabel("Gene")

plt.title("Top genes associated with favourable AML (OHSU)")

plt.tight_layout()

plt.savefig(os.path.join(FIGDIR, "ohsu_top_genes_barplot.png"), dpi=300)
plt.close()


# --------------------------------------------
# 4. SHARED GENE EXPRESSION COMPARISON
# --------------------------------------------

print("Generating shared gene comparison...")

shared = list(overlap_all)

if len(shared) > 0:

    shared_df = pd.DataFrame({
        "OHSU": ohsu_stats.loc[shared, "mean_diff"],
        "TCGA": tcga_stats.loc[shared, "mean_diff"],
        "TARGET": target_stats.loc[shared, "mean_diff"]
    })

    shared_df.plot(kind="bar", figsize=(8,5))

    plt.ylabel("Expression difference")
    plt.title("Shared favourable AML genes across datasets")

    plt.tight_layout()

    plt.savefig(os.path.join(FIGDIR, "shared_gene_barplot.png"), dpi=300)
    plt.close()


print("✅ New simplified figures saved")