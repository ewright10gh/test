
# COMPARE TOP GENES ACROSS OHSU, TCGA, TARGET


import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-GUI backend
import matplotlib.pyplot as plt
from joblib import load
from scipy.stats import mannwhitneyu, pearsonr

# Get the directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)

# 2022 OHSU (training cohort)
OHSU_EXPR = os.path.join(parent_dir, "cleaned", "ohsu_cleaned_expression.csv")
OHSU_LABELS = os.path.join(parent_dir, "cleaned", "ohsu_favorable_labels.csv")

# 2018 OHSU from VALIDATA (external validation cohort)
OHSU_2018_EXPR = os.path.join(parent_dir, "cleaned", "validata_cleaned_expression.csv")
OHSU_2018_PRED = os.path.join(parent_dir, "results", "validata_favourable_fusion_predictions.csv")

TCGA_EXPR = os.path.join(parent_dir, "cleaned", "tcga_cleaned_expression.csv")
TARGET_EXPR = os.path.join(parent_dir, "cleaned", "target_cleaned_expression.csv")

TCGA_PRED = os.path.join(parent_dir, "results", "TCGA_favourable_fusion_predictions.csv")
TARGET_PRED = os.path.join(parent_dir, "results", "TARGET_favourable_fusion_predictions.csv")

OUTDIR = os.path.join(parent_dir, "results")
os.makedirs(OUTDIR, exist_ok=True)


# LOAD DATA

print("Loading expression matrices...")

ohsu_2022 = pd.read_csv(OHSU_EXPR, index_col=0)
ohsu_2018 = pd.read_csv(OHSU_2018_EXPR, index_col=0)
tcga = pd.read_csv(TCGA_EXPR, index_col=0)
target = pd.read_csv(TARGET_EXPR, index_col=0)

print("OHSU 2022:", ohsu_2022.shape)
print("OHSU 2018 (VALIDATA):", ohsu_2018.shape)
print("TCGA:", tcga.shape)
print("TARGET:", target.shape)


# LOAD LABELS / PREDICTIONS

ohsu_2022_labels = pd.read_csv(OHSU_LABELS, index_col=0).iloc[:,0]
ohsu_2018_pred = pd.read_csv(OHSU_2018_PRED, index_col=0)
tcga_pred = pd.read_csv(TCGA_PRED, index_col=0)
target_pred = pd.read_csv(TARGET_PRED, index_col=0)

# Ensure valid binary values for 2018 prediction
ohsu_2018_pred["favourable_fusion_predicted"] = pd.to_numeric(
    ohsu_2018_pred["favourable_fusion_predicted"], errors="coerce"
).fillna(0).astype(int)

# Align sample order
ohsu_2022 = ohsu_2022.loc[ohsu_2022_labels.index]
ohsu_2018 = ohsu_2018.loc[ohsu_2018_pred.index]
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


# OHSU (2022 training data)

print("\nComputing OHSU 2022 stats...")

ohsu_2022_stats = compute_stats(ohsu_2022, ohsu_2022_labels)

top_ohsu_2022 = ohsu_2022_stats.sort_values("mean_diff", ascending=False).head(50)

# add gene column
top_ohsu_2022 = top_ohsu_2022.reset_index().rename(columns={"index": "gene"})

top_ohsu_2022.to_csv(
    os.path.join(OUTDIR, "OHSU_2022_top_genes.csv"),
    index=False
)

# OHSU (2018 VALIDATA)

print("Computing OHSU 2018 (VALIDATA) stats...")

ohsu_2018_stats = compute_stats(ohsu_2018, ohsu_2018_pred["favourable_fusion_predicted"])

top_ohsu_2018 = ohsu_2018_stats.sort_values("mean_diff", ascending=False).head(50)

top_ohsu_2018 = top_ohsu_2018.reset_index().rename(columns={"index": "gene"})

top_ohsu_2018.to_csv(
    os.path.join(OUTDIR, "OHSU_2018_top_genes.csv"),
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

from itertools import combinations

# Store all sets in a dictionary
gene_sets = {
    "OHSU 2022": set(top_ohsu_2022["gene"]),
    "OHSU 2018": set(top_ohsu_2018["gene"]),
    "TCGA": set(top_tcga["gene"]),
    "TARGET": set(top_target["gene"])
}

print("\nPairwise gene overlaps:")

# Loop through all 2-cohort combinations
for (name1, set1), (name2, set2) in combinations(gene_sets.items(), 2):
    overlap = set1 & set2
    print(f"{name1} ∩ {name2}: {len(overlap)}")

# Keep your 4-way overlap as well
overlap_all = set.intersection(*gene_sets.values())

print("\nAll four datasets:", len(overlap_all))
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
check_gene("MYH11", ohsu_2022_stats)
check_gene("MYH11", ohsu_2018_stats)
check_gene("MYH11", tcga_stats)
check_gene("MYH11", target_stats)


# COMPARE WITH RIDGE MODEL COEFFICIENTS
MODEL_PATH = os.path.join(parent_dir, "models", "ridge_eln_model.joblib")

print("\nComparing model coefficients...")

model = load(MODEL_PATH)

coefs = pd.Series(model.coef_[0], index=ohsu_2022.columns)

coef_df = pd.DataFrame({
    "ridge_coef": coefs,
    "expression_diff": ohsu_2022_stats["mean_diff"]
})

coef_df = coef_df.dropna()


# PLOT COEFFICIENT VS EXPRESSION DIFFERENCE

plt.figure(figsize=(7,7))

plt.scatter(
    coef_df["expression_diff"],
    coef_df["ridge_coef"],
    alpha=0.5
)


plt.xlabel("Expression difference Fav − Non-fav (Log2FC)")
plt.ylabel("Ridge coefficient (standardized log-odds effect size)")

plt.title("Model coefficient vs expression difference")

plt.tight_layout()
plt.subplots_adjust(left=0.18)

plt.savefig(os.path.join(OUTDIR,"ridge_vs_expression.png"))

plt.close()


# SAVE FULL TABLE


coef_df.to_csv(os.path.join(OUTDIR,"model_gene_weights_vs_expression.csv"))

print("\n✅ Gene comparison complete")
print("Results saved to:", OUTDIR)



# PUBLICATION-QUALITY FIGURE GENERATION


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

FIGDIR = os.path.join(parent_dir, "results", "figures")
os.makedirs(FIGDIR, exist_ok=True)


# 1. DATASET GENE HEATMAP


print("Generating dataset gene heatmap...")

heatmap_df = pd.DataFrame({
    "OHSU 2022": ohsu_2022_stats["mean_diff"],
    "OHSU 2018": ohsu_2018_stats["mean_diff"],
    "TCGA": tcga_stats["mean_diff"],
    "TARGET": target_stats["mean_diff"]
})

heatmap_df = heatmap_df.loc[
    heatmap_df.index
    .intersection(ohsu_2022_stats.index)
    .intersection(ohsu_2018_stats.index)
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
    cbar_kws={"label":"Expression difference (Log2FC)"}
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


# 2. EFFECT SIZE COMPARISON SCATTER


print("Generating dataset comparison scatter plot...")

common = ohsu_2022_stats.index.intersection(tcga_stats.index)

plt.figure(figsize=(6,6))

plt.scatter(
    ohsu_2022_stats.loc[common,"mean_diff"],
    tcga_stats.loc[common,"mean_diff"],
    alpha=0.6,
    s=40
)

plt.axhline(0, linestyle="--", linewidth=1)
plt.axvline(0, linestyle="--", linewidth=1)

plt.xlabel("OHSU Expression Difference Fav − Non-fav(Log2FC)")
plt.ylabel("TCGA Expression Difference Fav − Non-fav(Log2FC)")

plt.title("Consistency of gene signatures between cohorts", pad=12)

plt.tight_layout()

plt.savefig(os.path.join(FIGDIR,"gene_effect_scatter.png"))
plt.close()

# Compute Pearson correlation for gene effect consistency
pearson_r, pearson_p = pearsonr(
    ohsu_2022_stats.loc[common,"mean_diff"],
    tcga_stats.loc[common,"mean_diff"]
)
print(f"Pearson correlation (OHSU 2022 vs TCGA mean_diff across all common genes): {pearson_r:.4f}, p-value: {pearson_p:.4e}")

# 3. TOP GENES BARPLOT


print("Generating top gene barplot...")

top_genes = ohsu_2022_stats.sort_values("mean_diff", ascending=False).head(15)

plt.figure(figsize=(7,6))

top_genes["mean_diff"].sort_values().plot.barh(color="#ce3581")

plt.xlabel("Expression difference Fav − Non-fav (Log2FC)")
plt.ylabel("Gene")

plt.title("Top genes associated with favourable AML (OHSU)", pad=12)

plt.tight_layout()

plt.savefig(os.path.join(FIGDIR,"ohsu_top_genes_barplot.png"))
plt.close()


# 4. SHARED GENE EXPRESSION COMPARISON


print("Generating shared gene comparison...")

shared = list(overlap_all)

if len(shared) > 0:

    shared_df = pd.DataFrame({
        "OHSU 2022": ohsu_2022_stats.loc[shared,"mean_diff"],
        "OHSU 2018": ohsu_2018_stats.loc[shared,"mean_diff"],
        "TCGA": tcga_stats.loc[shared,"mean_diff"],
        "TARGET": target_stats.loc[shared,"mean_diff"]
    })

    # Compute Mann-Whitney U tests separately for shared genes
    mannwhitney_rows = []
    for gene in shared_df.index:
        fav_values = ohsu_2022.loc[ohsu_2022_labels == 1, gene].dropna()
        nonfav_values = ohsu_2022.loc[ohsu_2022_labels == 0, gene].dropna()
        if len(fav_values) > 0 and len(nonfav_values) > 0:
            u_stat, p_val = mannwhitneyu(fav_values, nonfav_values, alternative="two-sided")
        else:
            u_stat, p_val = np.nan, np.nan
        mannwhitney_rows.append({
            "gene": gene,
            "U_statistic": u_stat,
            "p_value": p_val
        })

    mannwhitney_df = pd.DataFrame(mannwhitney_rows)
    mannwhitney_df.to_csv(
        os.path.join(OUTDIR, "shared_genes_mannwhitney_results.csv"),
        index=False
    )

    # Plot only the mean_diff values
    # Sort genes by average expression difference (descending) for ordering
    shared_df['avg_diff'] = shared_df.mean(axis=1)
    shared_df = shared_df.sort_values('avg_diff', ascending=False).drop(columns='avg_diff')

    purple_pink_blue = sns.color_palette(["#D6C20C", "#ce3581", "#EC7608", "#1f77b4"])
    shared_df.plot(
        kind="bar",
        figsize=(8,6),
        color=purple_pink_blue
    )

    plt.ylabel("Expression difference Fav − Non-fav (Log2FC)")
    plt.xlabel("Gene")

    plt.title("Genes shared across AML cohorts", pad=12)

    plt.xticks(rotation=45, ha="right")

    plt.tight_layout()

    plt.savefig(os.path.join(FIGDIR,"shared_gene_barplot.png"))
    plt.close()

print("Gene comparison figures saved.")




# 5. FUSION SUBTYPE HEATMAP (OHSU)


print("\nGenerating fusion subtype heatmap...")

# ----------------------------------------------------------
# Load expression
# ----------------------------------------------------------

expr = pd.read_csv(
    os.path.join(parent_dir, "cleaned", "ohsu_cleaned_expression.csv"),
    index_col=0
)

if expr.shape[0] > expr.shape[1]:
    expr = expr.T

print("Expression shape:", expr.shape)

# ----------------------------------------------------------
# Load clinical
# ----------------------------------------------------------

clin = pd.read_csv(
    os.path.join(parent_dir, "aml_ohsu_2022_clinical_data.tsv"),
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

    if "RUNX1T1" in x or "8;21" in x:
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
# Load top genes (with error handling for missing files)
# ----------------------------------------------------------

try:
    ohsu_genes = pd.read_csv(
        os.path.join(parent_dir, "results", "OHSU_top_genes.csv")
    )["gene"].head(50).tolist()
except FileNotFoundError:
    print("Warning: OHSU_top_genes.csv not found. Skipping fusion heatmap generation or using fallback genes from heatmap_df.")
    ohsu_genes = top_heatmap_genes[:50] if 'top_heatmap_genes' in locals() else []

try:
    tcga_genes = pd.read_csv(
        os.path.join(parent_dir, "results", "TCGA_top_genes.csv")
    )["gene"].head(50).tolist()
except FileNotFoundError:
    print("Warning: TCGA_top_genes.csv not found. Using OHSU genes as fallback for TCGA genes in fusion heatmap.")
    tcga_genes = ohsu_genes.copy() if ohsu_genes else []

try:
    target_genes = pd.read_csv(
        os.path.join(parent_dir, "results", "TARGET_top_genes.csv")
    )["gene"].head(50).tolist()
except FileNotFoundError:
    print("Warning: TARGET_top_genes.csv not found. Using OHSU genes as fallback for TARGET genes in fusion heatmap.")
    target_genes = ohsu_genes.copy() if ohsu_genes else []

ohsu_genes = set(ohsu_genes) if ohsu_genes else set()
tcga_genes = set(tcga_genes) if tcga_genes else set(ohsu_genes)  # fallback to ohsu if empty
target_genes = set(target_genes) if target_genes else set(ohsu_genes)  # fallback to ohsu if empty

shared_genes = list(ohsu_genes & tcga_genes & target_genes)

if len(shared_genes) < 10:
    shared_genes = list(ohsu_genes | tcga_genes | target_genes)[:40]

print("Genes used:", len(shared_genes))


# FUSION TRANSCRIPTIONAL PROGRAM HEATMAP (IMPROVED)


print("Generating fusion transcriptional program heatmap...")

# compute mean expression per fusion subtype
fusion_means = expr_fusion.groupby(labels).mean()

# transpose so genes are rows
fusion_means = fusion_means.T

# select favorable genes from the dataset heatmap (genes with positive expression difference in OHSU)
favorable_top_genes = [gene for gene in top_heatmap_genes if ohsu_2022_stats.loc[gene, "mean_diff"] > 0]

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
    cbar_kws={"label":"Relative expression (z-score, SD units)"}
)

plt.title("Fusion subtype–specific transcriptional programs in AML")

plt.xlabel("Fusion subtype")
plt.ylabel("Gene")

plt.xticks(rotation=0)

plt.tight_layout()

plt.savefig(os.path.join(FIGDIR,"fusion_program_heatmap_clean.png"))
plt.close()

print("Fusion program heatmap saved.")


# FUSION FREQUENCY ACROSS DATASETS


print("Generating fusion frequency figure...")

# OHSU counts (2022)
ohsu_counts = labels.value_counts()

# OHSU 2018 VALIDATA predicted fusions
try:
    validata_fusion = pd.read_csv(
        os.path.join(parent_dir,"results","fusion_inference","VALIDATA_fusion_inference.csv")
    )
    validata_counts = validata_fusion["Predicted_fusion_program"].value_counts()
except FileNotFoundError:
    print("Warning: VALIDATA_fusion_inference.csv not found. Skipping VALIDATA in fusion frequency plot.")
    validata_counts = pd.Series(dtype=int)

# TCGA predicted fusions
try:
    tcga_fusion = pd.read_csv(
        os.path.join(parent_dir,"results","fusion_inference","TCGA_fusion_inference.csv")
    )
    tcga_counts = tcga_fusion["Predicted_fusion_program"].value_counts()
except FileNotFoundError:
    print("Warning: TCGA_fusion_inference.csv not found. Skipping TCGA in fusion frequency plot.")
    tcga_counts = pd.Series(dtype=int)

# TARGET predicted fusions
try:
    target_fusion = pd.read_csv(
        os.path.join(parent_dir,"results","fusion_inference","TARGET_fusion_inference.csv")
    )
    target_counts = target_fusion["Predicted_fusion_program"].value_counts()
except FileNotFoundError:
    print("Warning: TARGET_fusion_inference.csv not found. Skipping TARGET in fusion frequency plot.")
    target_counts = pd.Series(dtype=int)

# ----------------------------------------------------------
# Combine counts
# ----------------------------------------------------------

fusion_freq = pd.DataFrame({
    "OHSU 2022": ohsu_counts,
    "OHSU 2018": validata_counts,
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
    color=purple_pink_blue
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