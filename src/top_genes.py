
# COMPARE TOP GENES ACROSS OHSU, TCGA, TARGET


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

OUTDIR = os.path.join(BASE, "results")
os.makedirs(OUTDIR, exist_ok=True)


# LOAD DATA

print("Loading expression matrices...")

ohsu = pd.read_csv(OHSU_EXPR, index_col=0)
tcga = pd.read_csv(TCGA_EXPR, index_col=0)
target = pd.read_csv(TARGET_EXPR, index_col=0)

print("OHSU:", ohsu.shape)
print("TCGA:", tcga.shape)
print("TARGET:", target.shape)


# LOAD LABELS / PREDICTIONS

ohsu_labels = pd.read_csv(OHSU_LABELS, index_col=0).iloc[:,0]
tcga_pred = pd.read_csv(TCGA_PRED, index_col=0)
target_pred = pd.read_csv(TARGET_PRED, index_col=0)

# Align sample order
ohsu = ohsu.loc[ohsu_labels.index]
tcga = tcga.loc[tcga_pred.index]
target = target.loc[target_pred.index]


# FUNCTION TO COMPUTE TOP GENES

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


# OHSU

print("\nComputing OHSU stats...")

ohsu_stats = compute_stats(ohsu, ohsu_labels)

top_ohsu = ohsu_stats.sort_values("mean_diff", ascending=False).head(50)

# add gene column
top_ohsu = top_ohsu.reset_index().rename(columns={"index": "gene"})

top_ohsu.to_csv(
    os.path.join(OUTDIR, "OHSU_top_genes.csv"),
    index=False
)


# TCGA

print("Computing TCGA stats...")

tcga_stats = compute_stats(tcga, tcga_pred["ELN_predicted"])

top_tcga = tcga_stats.sort_values("mean_diff", ascending=False).head(50)

top_tcga = top_tcga.reset_index().rename(columns={"index": "gene"})

top_tcga.to_csv(
    os.path.join(OUTDIR, "TCGA_top_genes.csv"),
    index=False
)


# TARGET

print("Computing TARGET stats...")

target_stats = compute_stats(target, target_pred["ELN_predicted_class"])

top_target = target_stats.sort_values("mean_diff", ascending=False).head(50)

top_target = top_target.reset_index().rename(columns={"index": "gene"})

top_target.to_csv(
    os.path.join(OUTDIR, "TARGET_top_genes.csv"),
    index=False
)


# GENE OVERLAP

set_ohsu = set(top_ohsu["gene"])
set_tcga = set(top_tcga["gene"])
set_target = set(top_target["gene"])

overlap_all = set_ohsu & set_tcga & set_target
overlap_ohsu_tcga = set_ohsu & set_tcga
overlap_ohsu_target = set_ohsu & set_target

print("\nGene overlap:")
print("OHSU ∩ TCGA:", len(overlap_ohsu_tcga))
print("OHSU ∩ TARGET:", len(overlap_ohsu_target))
print("All three:", len(overlap_all))
print("Shared genes:", overlap_all)


# SAVE SHARED GENES WITH HEADER

pd.DataFrame({"gene": list(overlap_all)}).to_csv(
    os.path.join(OUTDIR, "shared_genes_all_datasets.csv"),
    index=False
)

# CHECK MYH11 RANKING


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


# COMPARE WITH RIDGE MODEL COEFFICIENTS


print("\nComparing model coefficients...")

model = load(MODEL_PATH)

coefs = pd.Series(model.coef_[0], index=ohsu.columns)

coef_df = pd.DataFrame({
    "ridge_coef": coefs,
    "expression_diff": ohsu_stats["mean_diff"]
})

coef_df = coef_df.dropna()


# PLOT COEFFICIENT VS EXPRESSION DIFFERENCE


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


# SAVE FULL TABLE


coef_df.to_csv(os.path.join(OUTDIR,"model_gene_weights_vs_expression.csv"))

print("\n✅ Gene comparison complete")
print("Results saved to:", OUTDIR)


# SIMPLE INTERPRETABLE FIGURES


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
print("\n--- Generating fusion subtype heatmap ---")

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# -------------------------------
# Load expression
# -------------------------------

expr = pd.read_csv(
    r"C:\Users\mba22ew\test\cleaned\ohsu_cleaned_expression.csv",
    index_col=0
)

# genes as rows → transpose if needed
if expr.shape[0] > expr.shape[1]:
    expr = expr.T

print("Expression shape:", expr.shape)

# -------------------------------
# Load clinical
# -------------------------------

import pandas as pd

clin = pd.read_csv(
    r"C:\Users\mba22ew\test\aml_ohsu_2022_clinical_data.tsv",
    sep="\t"
)

print(clin.columns)
clin = clin.set_index("Sample ID")

print("Clinical shape:", clin.shape)

# -------------------------------
# Map fusion types
# -------------------------------

def map_fusion(x):

    x = str(x).upper()

    if "PML" in x or "15;17" in x:
        return "PML_RARA"

    if "RUNX1" in x or "8;21" in x:
        return "RUNX1_RUNX1T1"

    if "CBFB" in x or "INV(16)" in x:
        return "CBFB_MYH11"

    return "OTHER"


clin["fusion_program"] = clin["Cancer Type Detailed"].apply(map_fusion)

# -------------------------------
# Align samples
# -------------------------------

common = expr.index.intersection(clin.index)

expr = expr.loc[common]
clin = clin.loc[common]

print("Samples aligned:", len(common))

# keep only fusion samples
fusion_mask = clin["fusion_program"] != "OTHER"

expr_fusion = expr.loc[fusion_mask]
labels = clin.loc[fusion_mask, "fusion_program"]

print("\nFusion counts:")
print(labels.value_counts())

# -------------------------------
# Load top genes from cohorts
# -------------------------------

ohsu_genes = pd.read_csv(
    r"C:\Users\mba22ew\test\results\OHSU_top_genes.csv"
)["gene"]

tcga_genes = pd.read_csv(
    r"C:\Users\mba22ew\test\results\TCGA_top_genes.csv"
)["gene"]

target_genes = pd.read_csv(
    r"C:\Users\mba22ew\test\results\TARGET_top_genes.csv"
)["gene"]

# take top 50 from each
ohsu_genes = set(ohsu_genes.head(50))
tcga_genes = set(tcga_genes.head(50))
target_genes = set(target_genes.head(50))

# -------------------------------
# Find shared genes
# -------------------------------

shared_genes = list(ohsu_genes & tcga_genes & target_genes)

print("\nShared predictive genes:", len(shared_genes))

# fallback if overlap small
if len(shared_genes) < 10:
    shared_genes = list(ohsu_genes | tcga_genes | target_genes)[:40]
    print("Using union genes:", len(shared_genes))

# -------------------------------
# Build heatmap matrix
# -------------------------------

heat = expr_fusion[shared_genes]

# Z-score normalisation
heat = (heat - heat.mean()) / heat.std()

# order samples by subtype
order = labels.sort_values().index
heat = heat.loc[order]

# -------------------------------
# Color bar for fusions
# -------------------------------

fusion_colors = {
    "PML_RARA": "#e41a1c",
    "RUNX1_RUNX1T1": "#377eb8",
    "CBFB_MYH11": "#4daf4a"
}

col_colors = labels.map(fusion_colors)

# -------------------------------
# Plot heatmap
# -------------------------------

g = sns.clustermap(
    heat.T,
    cmap="vlag",
    col_cluster=False,
    xticklabels=False,
    yticklabels=True,
    col_colors=col_colors,
    figsize=(10,8)
)

plt.savefig(
    r"C:\Users\mba22ew\test\results\fusion_subtype_heatmap.png",
    dpi=400
)

print("\nFusion subtype heatmap saved:")
print(r"C:\Users\mba22ew\test\results\fusion_subtype_heatmap.png")
