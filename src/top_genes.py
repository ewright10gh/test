
# COMPARE TOP GENES ACROSS OHSU, TCGA, TARGET


import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-GUI backend
import matplotlib.pyplot as plt
from joblib import load


BASE = r"C:\Users\mba22ew\test"

OHSU_EXPR = os.path.join(BASE, "cleaned", "ohsu_cleaned_expression.csv")
TCGA_EXPR = os.path.join(BASE, "cleaned", "tcga_cleaned_expression.csv")
TARGET_EXPR = os.path.join(BASE, "cleaned", "target_cleaned_expression.csv")

OHSU_LABELS = os.path.join(BASE, "cleaned", "ohsu_favorable_labels.csv")
TCGA_PRED = os.path.join(BASE, "results", "TCGA_favourable_fusion_predictions.csv")
TARGET_PRED = os.path.join(BASE, "results", "TARGET_favourable_fusion_predictions.csv")

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

tcga_stats = compute_stats(tcga, tcga_pred["favourable_fusion_predicted"])

top_tcga = tcga_stats.sort_values("mean_diff", ascending=False).head(50)

top_tcga = top_tcga.reset_index().rename(columns={"index": "gene"})

top_tcga.to_csv(
    os.path.join(OUTDIR, "TCGA_top_genes.csv"),
    index=False
)


# TARGET

print("Computing TARGET stats...")

target_stats = compute_stats(target, target_pred["favourable_fusion_predicted_class"])

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
MODEL_PATH = "models/ridge_eln_model.joblib"

print("\nComparing model coefficients...")

model = load(MODEL_PATH)

coefs = pd.Series(model.coef_[0], index=ohsu.columns)

coef_df = pd.DataFrame({
    "ridge_coef": coefs,
    "expression_diff": ohsu_stats["mean_diff"]
})

coef_df = coef_df.dropna()


# PLOT COEFFICIENT VS EXPRESSION DIFFERENCE

plt.figure(figsize=(7,7))

plt.scatter(
    coef_df["expression_diff"],
    coef_df["ridge_coef"],
    alpha=0.5
)


plt.xlabel("Expression difference (favourable - adverse)")
plt.ylabel("Ridge coefficient")

plt.title("Model weight vs biological expression difference")

plt.tight_layout()
plt.subplots_adjust(left=0.18)

plt.savefig(os.path.join(OUTDIR,"ridge_vs_expression.png"))

plt.close()


# SAVE FULL TABLE


coef_df.to_csv(os.path.join(OUTDIR,"model_gene_weights_vs_expression.csv"))

print("\n✅ Gene comparison complete")
print("Results saved to:", OUTDIR)


# ==========================================================
# PUBLICATION-QUALITY FIGURE GENERATION
# ==========================================================

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ----------------------------------------------------------
# GLOBAL STYLE
# ----------------------------------------------------------

sns.set_style("white")
sns.set_context("paper", font_scale=1.4)

plt.rcParams["figure.dpi"] = 120
plt.rcParams["savefig.dpi"] = 400
plt.rcParams["font.family"] = "Arial"

FIGDIR = os.path.join(BASE, "results", "figures")
os.makedirs(FIGDIR, exist_ok=True)

# ==========================================================
# 1. DATASET GENE HEATMAP
# ==========================================================

print("Generating dataset gene heatmap...")

heatmap_df = pd.DataFrame({
    "OHSU": ohsu_stats["mean_diff"],
    "TCGA": tcga_stats["mean_diff"],
    "TARGET": target_stats["mean_diff"]
})

heatmap_df = heatmap_df.loc[
    heatmap_df.index
    .intersection(ohsu_stats.index)
    .intersection(tcga_stats.index)
    .intersection(target_stats.index)
]

heatmap_df["abs_mean"] = heatmap_df.abs().mean(axis=1)

heatmap_df = (
    heatmap_df
    .sort_values("abs_mean", ascending=False)
    .head(20)
    .drop(columns="abs_mean")
)

plt.figure(figsize=(7,8))

sns.heatmap(
    heatmap_df,
    cmap="vlag",
    center=0,
    linewidths=0.5,
    linecolor="lightgrey",
    annot=True,
    fmt=".2f",
    cbar_kws={"label":"Expression difference"}
)

plt.title("Top AML transcriptional differences across datasets", pad=15)
plt.xlabel("Dataset")
plt.ylabel("Gene")

plt.yticks(rotation=0)

plt.tight_layout(rect=[0,0,1,0.97])

plt.savefig(os.path.join(FIGDIR,"gene_dataset_heatmap.png"))
plt.close()

# Store top genes from dataset heatmap for later use
top_heatmap_genes = heatmap_df.index.tolist()

# ==========================================================
# 2. EFFECT SIZE COMPARISON SCATTER
# ==========================================================

print("Generating dataset comparison scatter plot...")

common = ohsu_stats.index.intersection(tcga_stats.index)

plt.figure(figsize=(6,6))

plt.scatter(
    ohsu_stats.loc[common,"mean_diff"],
    tcga_stats.loc[common,"mean_diff"],
    alpha=0.6,
    s=40
)

plt.axhline(0, linestyle="--", linewidth=1)
plt.axvline(0, linestyle="--", linewidth=1)

plt.xlabel("OHSU gene effect size")
plt.ylabel("TCGA gene effect size")

plt.title("Consistency of gene signatures between cohorts", pad=12)

plt.tight_layout()

plt.savefig(os.path.join(FIGDIR,"gene_effect_scatter.png"))
plt.close()

# ==========================================================
# 3. TOP GENES BARPLOT
# ==========================================================

print("Generating top gene barplot...")

top_genes = ohsu_stats.sort_values("mean_diff", ascending=False).head(15)

plt.figure(figsize=(7,6))

top_genes["mean_diff"].sort_values().plot.barh(color="#4C72B0")

plt.xlabel("Expression difference (favourable − adverse)")
plt.ylabel("Gene")

plt.title("Top genes associated with favourable AML (OHSU)", pad=12)

plt.tight_layout()

plt.savefig(os.path.join(FIGDIR,"ohsu_top_genes_barplot.png"))
plt.close()

# ==========================================================
# 4. SHARED GENE EXPRESSION COMPARISON
# ==========================================================

print("Generating shared gene comparison...")

shared = list(overlap_all)

if len(shared) > 0:

    shared_df = pd.DataFrame({
        "OHSU": ohsu_stats.loc[shared,"mean_diff"],
        "TCGA": tcga_stats.loc[shared,"mean_diff"],
        "TARGET": target_stats.loc[shared,"mean_diff"]
    })

    shared_df.plot(
        kind="bar",
        figsize=(8,6),
        colormap="Set2"
    )

    plt.ylabel("Expression difference")
    plt.xlabel("Gene")

    plt.title("Genes shared across AML cohorts", pad=12)

    plt.xticks(rotation=45, ha="right")

    plt.tight_layout()

    plt.savefig(os.path.join(FIGDIR,"shared_gene_barplot.png"))
    plt.close()

print("Gene comparison figures saved.")



# ==========================================================
# 5. FUSION SUBTYPE HEATMAP (OHSU)
# ==========================================================

print("\nGenerating fusion subtype heatmap...")

# ----------------------------------------------------------
# Load expression
# ----------------------------------------------------------

expr = pd.read_csv(
    r"C:\Users\mba22ew\test\cleaned\ohsu_cleaned_expression.csv",
    index_col=0
)

if expr.shape[0] > expr.shape[1]:
    expr = expr.T

print("Expression shape:", expr.shape)

# ----------------------------------------------------------
# Load clinical
# ----------------------------------------------------------

clin = pd.read_csv(
    r"C:\Users\mba22ew\test\aml_ohsu_2022_clinical_data.tsv",
    sep="\t"
)

clin = clin.set_index("Sample ID")

print("Clinical shape:", clin.shape)

# ----------------------------------------------------------
# Map fusion types
# ----------------------------------------------------------

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

# ----------------------------------------------------------
# Align samples
# ----------------------------------------------------------

common = expr.index.intersection(clin.index)

expr = expr.loc[common]
clin = clin.loc[common]

fusion_mask = clin["fusion_program"] != "OTHER"

expr_fusion = expr.loc[fusion_mask]
labels = clin.loc[fusion_mask,"fusion_program"]

print("\nFusion counts:")
print(labels.value_counts())

# ----------------------------------------------------------
# Load top genes
# ----------------------------------------------------------

ohsu_genes = pd.read_csv(
    r"C:\Users\mba22ew\test\results\OHSU_top_genes.csv"
)["gene"]

tcga_genes = pd.read_csv(
    r"C:\Users\mba22ew\test\results\TCGA_top_genes.csv"
)["gene"]

target_genes = pd.read_csv(
    r"C:\Users\mba22ew\test\results\TARGET_top_genes.csv"
)["gene"]

ohsu_genes = set(ohsu_genes.head(50))
tcga_genes = set(tcga_genes.head(50))
target_genes = set(target_genes.head(50))

shared_genes = list(ohsu_genes & tcga_genes & target_genes)

if len(shared_genes) < 10:
    shared_genes = list(ohsu_genes | tcga_genes | target_genes)[:40]

print("Genes used:", len(shared_genes))

# ==========================================================
# FUSION TRANSCRIPTIONAL PROGRAM HEATMAP (IMPROVED)
# ==========================================================

print("Generating fusion transcriptional program heatmap...")

# compute mean expression per fusion subtype
fusion_means = expr_fusion.groupby(labels).mean()

# transpose so genes are rows
fusion_means = fusion_means.T

# select favorable genes from the dataset heatmap (genes with positive expression difference in OHSU)
favorable_top_genes = [gene for gene in top_heatmap_genes if ohsu_stats.loc[gene, "mean_diff"] > 0]

if len(favorable_top_genes) >= 5:
    selected_genes = favorable_top_genes
else:
    # Fallback: select top informative genes based on variance
    gene_variance = fusion_means.var(axis=1)
    selected_genes = gene_variance.sort_values(ascending=False).head(20).index.tolist()

fusion_means = fusion_means.loc[selected_genes]

# z-score genes
fusion_means = (fusion_means - fusion_means.mean(axis=1).values.reshape(-1,1)) / fusion_means.std(axis=1).values.reshape(-1,1)

plt.figure(figsize=(6,8))

sns.heatmap(
    fusion_means,
    cmap="vlag",
    center=0,
    linewidths=0.5,
    linecolor="lightgrey",
    cbar_kws={"label":"Z-scored expression"}
)

plt.title("Distinct transcriptional programs driven by AML fusion subtypes")

plt.xlabel("Fusion subtype")
plt.ylabel("Gene")

plt.xticks(rotation=0)

plt.tight_layout()

plt.savefig(os.path.join(FIGDIR,"fusion_program_heatmap_clean.png"))
plt.close()

print("Fusion program heatmap saved.")

# ==========================================================
# FUSION FREQUENCY ACROSS DATASETS
# ==========================================================

print("Generating fusion frequency figure...")

# OHSU counts
ohsu_counts = labels.value_counts()

# TCGA predicted fusions
tcga_fusion = pd.read_csv(
    os.path.join(BASE,"results","fusion_inference","TCGA_fusion_inference.csv")
)

tcga_counts = tcga_fusion["Predicted_fusion_program"].value_counts()

# TARGET predicted fusions
target_fusion = pd.read_csv(
    os.path.join(BASE,"results","fusion_inference","TARGET_fusion_inference.csv")
)

target_counts = target_fusion["Predicted_fusion_program"].value_counts()

# ----------------------------------------------------------
# Combine counts
# ----------------------------------------------------------

fusion_freq = pd.DataFrame({
    "OHSU": ohsu_counts,
    "TCGA": tcga_counts,
    "TARGET": target_counts
})

# convert to percentages
fusion_freq_percent = fusion_freq.div(
    fusion_freq.sum(axis=0),
    axis=1
) * 100

# ----------------------------------------------------------
# Plot
# ----------------------------------------------------------

fusion_freq_percent.T.plot(
    kind="bar",
    figsize=(7,5),
    colormap="Set1"
)

plt.ylabel("Frequency (%)")
plt.xlabel("Dataset")

plt.title(
    "Frequency of major AML fusion subtypes across cohorts",
    pad=15
)

plt.xticks(rotation=0)

plt.legend(title="Fusion subtype")

plt.tight_layout()

plt.savefig(
    os.path.join(FIGDIR,"fusion_frequency_across_datasets.png")
)

plt.close()

print("Fusion frequency figure saved.")