import os
import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

# Get the directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)


# FILE PATHS


OHSU_EXPR = os.path.join(parent_dir, "cleaned", "ohsu_cleaned_expression.csv")
OHSU_CLIN = os.path.join(parent_dir, "aml_ohsu_2022_clinical_data.tsv")

TCGA_EXPR = os.path.join(parent_dir, "cleaned", "tcga_cleaned_expression.csv")
TCGA_PRED = os.path.join(parent_dir, "results", "TCGA_favourable_fusion_predictions.csv")

TARGET_EXPR = os.path.join(parent_dir, "cleaned", "target_cleaned_expression.csv")
TARGET_PRED = os.path.join(parent_dir, "results", "TARGET_favourable_fusion_predictions.csv")

VALIDATA_EXPR = os.path.join(parent_dir,"cleaned","validata_cleaned_expression.csv")
VALIDATA_PRED = os.path.join(parent_dir,"results","validata_favourable_fusion_predictions.csv")

OUTDIR = os.path.join(parent_dir,"results","fusion_inference")
os.makedirs(OUTDIR,exist_ok=True)


# LOAD DATA


print("Loading OHSU expression...")
ohsu_expr = pd.read_csv(OHSU_EXPR,index_col=0)

print("Loading OHSU clinical...")
ohsu_clin = pd.read_csv(OHSU_CLIN,sep="\t",low_memory=False)
ohsu_clin = ohsu_clin.set_index("Sample ID")

print("Loading TCGA expression...")
tcga_expr = pd.read_csv(TCGA_EXPR,index_col=0)

print("Loading TCGA predictions...")
tcga_pred = pd.read_csv(TCGA_PRED,index_col=0)

print("Loading TARGET expression...")
target_expr = pd.read_csv(TARGET_EXPR,index_col=0)

print("Loading TARGET predictions...")
target_pred = pd.read_csv(TARGET_PRED,index_col=0)

print("Loading VALIDATA expression...")
validata_expr = pd.read_csv(VALIDATA_EXPR,index_col=0)

print("Loading VALIDATA predictions...")
validata_pred = pd.read_csv(VALIDATA_PRED,index_col=0)


# ALIGN OHSU


ohsu_clin = ohsu_clin.loc[ohsu_expr.index]

fusion_series = ohsu_clin["Cancer Type Detailed"].astype(str).str.upper()

fusion_group = pd.Series("OTHER",index=ohsu_expr.index)

fusion_group[fusion_series.str.contains("PML-RARA",na=False)] = "PML_RARA"
fusion_group[fusion_series.str.contains("RUNX1-RUNX1T1",na=False)] = "RUNX1_RUNX1T1"
fusion_group[fusion_series.str.contains("CBFB-MYH11",na=False)] = "CBFB_MYH11"

print("\nOHSU Fusion counts:")
print(fusion_group.value_counts())


# COMPUTE FUSION CENTROIDS


centroids = {}

for fusion in ["PML_RARA", "RUNX1_RUNX1T1", "CBFB_MYH11"]:

    subset = ohsu_expr[fusion_group == fusion]

    if len(subset) == 0:
        continue

    centroids[fusion] = subset.mean()

centroids = pd.DataFrame(centroids)

print("\nCentroid matrix shape:", centroids.shape)


# =========================================================
# 2. FUSION INFERENCE FUNCTION (COSINE SIMILARITY)
# =========================================================

def infer_fusion(expr, pred, dataset_name):

    print(f"\nProcessing {dataset_name}")

    # -----------------------------
    # Ensure Sample ID index
    # -----------------------------
    if "Sample_ID" in pred.columns:
        pred = pred.set_index("Sample_ID")

    # -----------------------------
    # Find probability column
    # -----------------------------
    prob_col = None
    for c in pred.columns:
        if "prob" in c.lower():
            prob_col = c
            break

    if prob_col is None:
        raise ValueError("No probability column found")

    print("Using probability column:", prob_col)

    # -----------------------------
    # Align samples
    # -----------------------------
    common_samples = expr.index.intersection(pred.index)

    expr = expr.loc[common_samples]
    pred = pred.loc[common_samples]

    print("Aligned samples:", len(common_samples))

    # -----------------------------
    # Select favourable samples
    # -----------------------------
    fav_mask = pred[prob_col] > 0.675
    fav_expr = expr.loc[fav_mask]

    print(f"{dataset_name} favourable samples:", fav_expr.shape[0])

    if fav_expr.shape[0] == 0:
        print("No favourable samples found")
        return None

    # -----------------------------
    # Align genes
    # -----------------------------
    common_genes = fav_expr.columns.intersection(centroids.index)

    fav_expr = fav_expr[common_genes]
    cent = centroids.loc[common_genes]

    print("Genes used:", len(common_genes))

    # -----------------------------
    # COSINE SIMILARITY
    # -----------------------------
    sim = cosine_similarity(fav_expr, cent.T)

    sim_df = pd.DataFrame(
        sim,
        index=fav_expr.index,
        columns=cent.columns
    )

    # Assign subtype
    sim_df["Predicted_fusion_program"] = sim_df.idxmax(axis=1)

    # Save
    outfile = os.path.join(OUTDIR, f"{dataset_name}_fusion_inference.csv")
    sim_df.to_csv(outfile)

    print("Saved:", outfile)

    print("\nFusion distribution:")
    print(sim_df["Predicted_fusion_program"].value_counts())

    return sim_df


# =========================================================
# 3. EVALUATION FUNCTION (CONFUSION MATRIX + REPORT)
# =========================================================

from sklearn.metrics import classification_report, confusion_matrix
import seaborn as sns
import matplotlib.pyplot as plt

def evaluate_fusion_predictions(sim_df, true_labels, dataset_name):

    print(f"\nEvaluating {dataset_name}")

    # Align
    common = sim_df.index.intersection(true_labels.index)

    sim_df = sim_df.loc[common]
    true_labels = true_labels.loc[common]

    y_true = true_labels
    y_pred = sim_df["Predicted_fusion_program"]

    # Classification report
    print(f"\n=== {dataset_name} Fusion Classification Report ===")
    print(classification_report(y_true, y_pred))

    # Confusion matrix
    labels = ["PML_RARA", "RUNX1_RUNX1T1", "CBFB_MYH11"]

    cm = confusion_matrix(y_true, y_pred, labels=labels)

    plt.figure(figsize=(5, 4))
    sns.heatmap(
        cm,
        annot=True,
        fmt="d",
        cmap="Blues",
        xticklabels=labels,
        yticklabels=labels
    )

    plt.xlabel("Predicted")
    plt.ylabel("True")
    plt.title(f"{dataset_name} Fusion Subtype Confusion Matrix")

    plt.tight_layout()
    plt.savefig(f"{dataset_name}_fusion_confusion_matrix.png", dpi=300)
    plt.show()


# =========================================================
# 4. OPTIONAL: CORRELATION-BASED INFERENCE
# =========================================================

def infer_fusion_correlation(expr, pred, dataset_name):

    print(f"\nProcessing {dataset_name} (correlation)")

    if "Sample_ID" in pred.columns:
        pred = pred.set_index("Sample_ID")

    prob_col = [c for c in pred.columns if "prob" in c.lower()][0]

    common_samples = expr.index.intersection(pred.index)

    expr = expr.loc[common_samples]
    pred = pred.loc[common_samples]

    fav_expr = expr.loc[pred[prob_col] > 0.675]

    common_genes = fav_expr.columns.intersection(centroids.index)

    fav_expr = fav_expr[common_genes]
    cent = centroids.loc[common_genes]

    # Correlation matrix
    corr = np.corrcoef(fav_expr.values, cent.T.values)[:len(fav_expr), len(fav_expr):]

    corr_df = pd.DataFrame(
        corr,
        index=fav_expr.index,
        columns=cent.columns
    )

    corr_df["Predicted_fusion_program"] = corr_df.idxmax(axis=1)

    return corr_df


# =========================================================
# 5. OPTIONAL: KNN BASELINE
# =========================================================

from sklearn.neighbors import KNeighborsClassifier

def run_knn_baseline(fav_expr):

    knn = KNeighborsClassifier(n_neighbors=1)

    X_train = ohsu_expr.loc[
        fusion_group.isin(["PML_RARA","RUNX1_RUNX1T1","CBFB_MYH11"])
    ]

    y_train = fusion_group.loc[X_train.index]

    knn.fit(X_train, y_train)

    y_pred = knn.predict(fav_expr)

    return pd.Series(y_pred, index=fav_expr.index)


# =========================================================
# 6. RUN ANALYSIS
# =========================================================

tcga_results = infer_fusion(tcga_expr, tcga_pred, "TCGA")
validata_results = infer_fusion(validata_expr, validata_pred, "VALIDATA")
target_results = infer_fusion(target_expr, target_pred, "TARGET")


# =========================================================
# 7. EVALUATE (ONLY WHERE TRUE LABELS EXIST)
# =========================================================

def assign_fusion(label):
    label = str(label).upper()
    if "PML-RARA" in label or "15;17" in label:
        return "PML_RARA"
    if "RUNX1-RUNX1T1" in label or "8;21" in label:
        return "RUNX1_RUNX1T1"
    if "CBFB-MYH11" in label or "INV(16)" in label:
        return "CBFB_MYH11"
    return np.nan

# Load TCGA clinical data from the workspace
tcga_clinical = pd.read_csv(os.path.join(parent_dir, "laml_tcga_pub_clinical_data.tsv"), sep="\t", low_memory=False)

tcga_clinical = tcga_clinical.set_index("Sample ID")

tcga_clinical["fusion"] = tcga_clinical["Cytogenetics"].apply(assign_fusion)

# Keep only the 3 fusions
tcga_true_fusions = tcga_clinical["fusion"].dropna()

print(tcga_true_fusions.value_counts())

validata_clinical = pd.read_csv(os.path.join(parent_dir, "aml_ohsu_2018_clinical_data.tsv"), sep="\t", low_memory=False)

if "sampleId" in validata_clinical.columns:
    validata_clinical = validata_clinical.set_index("sampleId")
else:
    validata_clinical = validata_clinical.set_index("Sample ID")

validata_clinical["fusion"] = validata_clinical["CANCER_TYPE_DETAILED"].apply(assign_fusion)

validata_true_fusions = validata_clinical["fusion"].dropna()

print(validata_true_fusions.value_counts())

print(tcga_true_fusions.index.isin(tcga_results.index).mean())

if tcga_results is not None:
    evaluate_fusion_predictions(tcga_results, tcga_true_fusions, "TCGA")
else:
    print("TCGA fusion inference results missing; skipping evaluation.")

if validata_results is not None:
    evaluate_fusion_predictions(validata_results, validata_true_fusions, "OHSU 2018")
else:
    print("VALIDATA fusion inference results missing; skipping evaluation.")

if tcga_results is not None:
    tcga_eval = tcga_results.loc[tcga_true_fusions.index.intersection(tcga_results.index)]

tcga_true = tcga_true_fusions.loc[tcga_eval.index]

print("\nDone.")